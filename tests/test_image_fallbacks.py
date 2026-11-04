"""Offline image-discovery, product-selection, and fallback-chain regressions."""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.parse import parse_qs, urlsplit

import numpy as np
import requests
from astropy.io import fits
from astropy.io.votable import from_table, writeto
from astropy.table import Table

import get_finding_chart as chart


# Encode synthetic service metadata with FIELD IDs different from FIELD names.
def metadata_response(rows, status="OK"):
    table = Table(rows=rows)
    document = from_table(table)
    for index, field in enumerate(document.get_first_table().fields):
        field.ID = f"field{index}"
    output = io.BytesIO()
    writeto(document, output)
    xml = output.getvalue().replace(
        b"</RESOURCE>", f'<INFO name="QUERY_STATUS" value="{status}"/></RESOURCE>'.encode()
    )
    response = requests.Response()
    response.status_code = 200
    response._content = xml
    return response


# Produce an actual celestial image for validation and final-field coverage checks.
def write_image(path, *, value=7., ra=150., side=100):
    wcs = chart.make_pa_wcs(ra, 20., side, 1., 0.)
    fits.PrimaryHDU(np.full((side, side), value), header=wcs.to_header()).writeto(path, overwrite=True)


# Check product filtering, cutout construction, and final-field coverage.
class ImageDiscoveryTests(unittest.TestCase):
    # Reject malformed tables and VO errors even if HTTP returned success.
    def test_sia_errors(self):
        for status in ("ERROR", "OVERFLOW"):
            reply = metadata_response([{"access_url": "https://example.invalid/image.fits"}], status)
            with self.subTest(status=status), patch.object(chart.requests, "get", return_value=reply):
                with self.assertRaisesRegex(RuntimeError, status):
                    chart.query_image_sia("https://example.invalid/sia", {})
        reply = metadata_response([{"wrong_column": "value"}])
        with patch.object(chart.requests, "get", return_value=reply):
            with self.assertRaisesRegex(RuntimeError, "access_url"):
                chart.query_image_sia("https://example.invalid/sia", {})

    # Retain only same-band stacked science images, including the DES Y filter.
    def test_noirlab_product_selection(self):
        rows = [
            {"access_url": "https://example.invalid/g.fits", "obs_bandpass": "g DECam", "prodtype": "image", "proctype": "Stack"},
            {"access_url": "https://example.invalid/mask.fits", "obs_bandpass": "Y DECam", "prodtype": "dqmask", "proctype": "Stack"},
            {"access_url": "https://example.invalid/raw.fits", "obs_bandpass": "Y DECam", "prodtype": "image", "proctype": "InstCal"},
            {"access_url": "https://example.invalid/y.fits", "obs_bandpass": "Y DECam", "prodtype": "image", "proctype": "Stack"},
        ]
        with patch.object(chart.requests, "get", return_value=metadata_response(rows)) as get:
            with patch.object(chart, "download_image_candidates") as download:
                chart.get_noirlab_image("des", 150., 60., 6., "y", Path("unused.fits"), 4., 1., 0.)
        self.assertEqual(download.call_args.args[0], ["https://example.invalid/y.fits"])
        self.assertEqual(get.call_args.args[0], chart.NOIRLAB_SIA_URLS["des"])
        self.assertAlmostEqual(float(get.call_args.kwargs["params"]["SIZE"].split(",")[0]), .2)

    # IRSA discovers the All-Sky collection and converts Ks metadata to K cutouts.
    def test_irsa_collection_band_and_angular_cutout(self):
        rows = [
            {"access_url": "https://example.invalid/ibe/data/twomass/image.fits.gz?existing=yes", "energy_bandpassname": "Ks"},
            {"access_url": "https://example.invalid/ibe/data/twomass/j.fits.gz", "energy_bandpassname": "J"},
        ]
        with patch.object(chart.requests, "get", return_value=metadata_response(rows)) as get:
            with patch.object(chart, "download_image_candidates") as download:
                chart.get_2mass_irsa(150., 20., 6., "K", Path("unused.fits"), 4., 1., 0.)
        params = get.call_args.kwargs["params"]
        self.assertEqual(params["COLLECTION"], "twomass_allsky")
        self.assertEqual(params["BAND"], "2.2e-6")
        urls = download.call_args.args[0]
        self.assertEqual(len(urls), 1)
        cutout = parse_qs(urlsplit(urls[0]).query)
        self.assertEqual(cutout["size"], ["6.0000000000arcmin"])
        self.assertEqual(cutout["center"], ["150.0000000000,20.0000000000deg"])
        self.assertEqual(cutout["gzip"], ["false"])
        self.assertEqual(cutout["existing"], ["yes"])

    # Broken and off-center images are skipped until a real usable image is found.
    def test_candidate_validation_and_center_coverage(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "output.fits"

            # Model a failed download, a wrong field, and a valid field in order.
            def download(url, path):
                if url == "broken":
                    path.write_bytes(b"not FITS")
                else:
                    write_image(path, ra=160. if url == "elsewhere" else 150.)

            with patch.object(chart, "download_file", side_effect=download) as get:
                with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
                    chart.download_image_candidates(
                        ["broken", "elsewhere", "valid", "valid"], output, 150., 20., 1., 1., 0., "synthetic"
                    )
            self.assertEqual(get.call_count, 3)
            data, _ = chart.read_fits_image(output)
            self.assertTrue(np.all(data == 7.))

    # A later complete image takes precedence over a usable partial candidate.
    def test_complete_candidate_preferred_over_partial(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "output.fits"

            # Partial images cover the center but leave missing chart edges.
            def download(url, path):
                write_image(path, side=40 if url == "partial" else 100, value=1 if url == "partial" else 2)

            with patch.object(chart, "download_file", side_effect=download):
                with contextlib.redirect_stdout(io.StringIO()):
                    chart.download_image_candidates(
                        ["partial", "complete"], output, 150., 20., 1., 1., 0., "synthetic"
                    )
            self.assertTrue(np.all(chart.read_fits_image(output)[0] == 2.))

    # Warn when a native image only partly covers the requested chart.
    def test_partial_image_warning(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "output.fits"
            warnings = io.StringIO()

            # Supply one small, centered image with incomplete final-field coverage.
            def download(url, path):
                write_image(path, side=40)

            with patch.object(chart, "download_file", side_effect=download):
                with contextlib.redirect_stderr(warnings), contextlib.redirect_stdout(io.StringIO()):
                    chart.download_image_candidates(["partial"], output, 150., 20., 1., 1., 0., "synthetic")
            self.assertIn("incomplete final-field coverage", warnings.getvalue())
            self.assertTrue(output.is_file())

    # Failed discovery never overwrites a pre-existing output with unusable data.
    def test_no_matches_preserve_output(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "output.fits"
            output.write_bytes(b"existing")
            with self.assertRaisesRegex(RuntimeError, "No matching"):
                chart.download_image_candidates([], output, 150., 20., 1., 1., 0., "synthetic")
            self.assertEqual(output.read_bytes(), b"existing")


# Verify fallback order, backend reporting, and preservation of custom layers.
class ImageFallbackRoutingTests(unittest.TestCase):
    # Run retrieval using the same geometry for all route checks.
    def retrieve(self, survey, path, layer="ls-dr10"):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return chart.retrieve_image(survey, 150., 20., 100, 1.5, 1., 1., 0., "r" if survey != "2mass" else "J", path, layer)

    # A healthy primary remains preferred for Legacy, DES, and SDSS.
    def test_healthy_primaries_skip_noirlab(self):
        for survey in ("legacy", "des", "sdss"):
            with self.subTest(survey=survey), tempfile.TemporaryDirectory() as folder:
                output = Path(folder) / "output.fits"

                # Produce valid primary-service images in the existing call signature.
                def primary(*args):
                    write_image(args[4] if survey == "legacy" else args[5])

                with patch.object(chart, "get_legacy", side_effect=primary), patch.object(chart, "get_hips2fits", side_effect=primary):
                    with patch.object(chart, "get_noirlab_image") as alternative:
                        result = self.retrieve(survey, output)
                alternative.assert_not_called()
                self.assertNotIn("NOIRLab", result.backend)

    # Service failures dispatch to the corresponding NOIRLab survey collection.
    def test_noirlab_fallback_routes(self):
        for survey in ("legacy", "des", "sdss"):
            with self.subTest(survey=survey), tempfile.TemporaryDirectory() as folder:
                output = Path(folder) / "output.fits"

                # Preserve the alternate image as native FITS for --fits retention.
                def alternative(*args):
                    write_image(args[5])

                with patch.object(chart, "get_legacy", side_effect=requests.Timeout()), patch.object(chart, "get_hips2fits", side_effect=requests.Timeout()):
                    with patch.object(chart, "get_noirlab_image", side_effect=alternative) as alternate:
                        result = self.retrieve(survey, output)
                self.assertEqual(alternate.call_args.args[0], survey)
                self.assertIn("NOIRLab", result.backend)
                self.assertFalse(result.from_hips2fits)
                self.assertTrue(output.is_file())

    # Preserve explicit model, residual, and non-DR10 layer semantics on failure.
    def test_custom_legacy_layer_does_not_switch_products(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(chart, "get_legacy", side_effect=requests.Timeout()):
                with patch.object(chart, "get_noirlab_image") as alternative:
                    with self.assertRaisesRegex(RuntimeError, "only the ls-dr10"):
                        self.retrieve("legacy", Path(folder) / "output.fits", "ls-dr10-model")
            alternative.assert_not_called()

    # Two primary failures reach IRSA without querying any stellar catalog.
    def test_2mass_reaches_irsa_third(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "output.fits"
            order = []

            # Record the previous services failing in priority order.
            def skyview(*args):
                order.append("SkyView")
                raise requests.Timeout("SkyView down")

            # Record the existing CDS alternative failing before IRSA is attempted.
            def hips(*args):
                order.append("CDS")
                raise requests.Timeout("CDS down")

            # Supply an actual native image from the final alternative.
            def irsa(*args):
                order.append("IRSA")
                write_image(args[4])

            with patch.object(chart, "get_2mass_skyview", side_effect=skyview), patch.object(chart, "get_hips2fits", side_effect=hips):
                with patch.object(chart, "get_2mass_irsa", side_effect=irsa):
                    result = self.retrieve("2mass", output)
            self.assertEqual(order, ["SkyView", "CDS", "IRSA"])
            self.assertIn("IRSA", result.backend)
            self.assertFalse(result.from_hips2fits)

    # Keep both existing preferred 2MASS paths without unnecessary IRSA requests.
    def test_2mass_healthy_services_skip_irsa(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "output.fits"

            # Write a fully covered primary SkyView field.
            def skyview(*args):
                write_image(args[5])

            with patch.object(chart, "get_2mass_skyview", side_effect=skyview):
                with patch.object(chart, "get_2mass_irsa") as irsa, patch.object(chart, "get_hips2fits") as hips:
                    result = self.retrieve("2mass", output)
            irsa.assert_not_called()
            hips.assert_not_called()
            self.assertIn("SkyView", result.backend)

    # A working CDS image still wins before trying the newly added IRSA route.
    def test_2mass_healthy_cds_skips_irsa(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "output.fits"

            # Supply valid FITS in the existing HiPS function's output argument.
            def hips(*args):
                write_image(args[5])

            with patch.object(chart, "get_2mass_skyview", side_effect=requests.Timeout()):
                with patch.object(chart, "get_hips2fits", side_effect=hips), patch.object(chart, "get_2mass_irsa") as irsa:
                    result = self.retrieve("2mass", output)
            irsa.assert_not_called()
            self.assertTrue(result.from_hips2fits)

    # Exercise discovery, selection, FITS validation, rendering, and native FITS retention.
    def test_full_des_fallback_chart_and_fits(self):
        rows = [{
            "access_url": "https://example.invalid/science.fits", "obs_bandpass": "r DECam",
            "prodtype": "image", "proctype": "Stack",
        }]
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder) / "chart.png"
            retained = Path(folder) / "native.fits"
            argv = [
                "get_finding_chart.py", "-c", "150", "20", "-s", "des", "-f", "1",
                "-o", str(output), "--fits", str(retained),
            ]

            # Replace only remote traffic; keep all chart-processing code active.
            def download(url, path):
                write_image(path)

            with patch.object(sys, "argv", argv), patch.object(chart, "get_hips2fits", side_effect=requests.Timeout()):
                with patch.object(chart.requests, "get", return_value=metadata_response(rows)), patch.object(chart, "download_file", side_effect=download):
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        self.assertEqual(chart.main(), 0)
            self.assertTrue(output.read_bytes().startswith(b"\x89PNG"))
            self.assertTrue(np.all(chart.read_fits_image(retained)[0] == 7.))

    # Report all failure causes when every 2MASS service is unavailable.
    def test_all_2mass_services_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            with patch.object(chart, "get_2mass_skyview", side_effect=requests.Timeout("sky-error")), patch.object(chart, "get_hips2fits", side_effect=requests.Timeout("cds-error")):
                with patch.object(chart, "get_2mass_irsa", side_effect=requests.Timeout("irsa-error")):
                    with self.assertRaisesRegex(RuntimeError, "sky-error.*cds-error.*irsa-error"):
                        self.retrieve("2mass", Path(folder) / "output.fits")


# Allow direct execution with the standard unittest runner.
if __name__ == "__main__":
    unittest.main()
