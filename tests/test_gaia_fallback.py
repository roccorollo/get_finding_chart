"""Offline regressions for Gaia TAP failover and equivalent CSV/VOTable data."""

import contextlib
import io
import sys
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.votable import from_table, writeto
from astropy.table import MaskedColumn, Table

import get_finding_chart as chart


# Keep synthetic Gaia identifiers beyond the exact integer range of float64.
SOURCE_ID = 1234567890123456789
HEADER = "source_id,ra,dec,pmra,pmdec,ref_epoch,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag\n"
CSV = HEADER + f"{SOURCE_ID},150.001,20,100,-20,2016,18,18.5,17.5\n"


# Construct an HTTP response without contacting a remote TAP server.
def response(body, status=200):
    result = requests.Response()
    result.status_code = status
    result._content = body.encode() if isinstance(body, str) else body
    result.encoding = "utf-8"
    return result


# Build realistic TABLEDATA with differing FIELD IDs and optional masked cells.
def votable(empty=False, masked=False):
    table = Table(
        rows=[(SOURCE_ID, 150.001, 20., 100., -20., 2016., 18., 18.5, 17.5)],
        names=HEADER.strip().split(","),
    )
    if masked:
        for name in ("pmra", "pmdec", "phot_bp_mean_mag"):
            table.replace_column(name, MaskedColumn(table[name], mask=[True]))
    if empty:
        table = table[:0]
    document = from_table(table)
    for index, field in enumerate(document.get_first_table().fields):
        field.ID = f"col{index}"
    output = io.BytesIO()
    writeto(document, output)
    return output.getvalue()


# Verify service routing without changing the established acquisition calculation.
class GaiaFallbackTests(unittest.TestCase):
    # Capture messages and requests while querying a synthetic sky position.
    def query(self, replies):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(chart.requests, "post", side_effect=replies) as post:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                result = chart.query_gaia_candidates(
                    SkyCoord(150., 20., unit="deg"), chart.FaintSelection(1, 16., 19.), 2026.
                )
        return result, post, stdout.getvalue(), stderr.getvalue()

    # Prefer ESA without making an unnecessary mirror request.
    def test_esa_success(self):
        candidates, post, output, warnings = self.query([response(CSV)])
        self.assertEqual(post.call_count, 1)
        self.assertEqual(candidates[0].source_id, str(SOURCE_ID))
        self.assertTrue(candidates[0].propagated)
        self.assertIn("Gaia DR3 backend: ESA", output)
        self.assertEqual(warnings, "")

    # Exercise timeouts, connection errors, HTTP errors, and malformed CSV/XML.
    def test_failure_types_fall_back_to_equivalent_votable(self):
        expected = self.query([response(CSV)])[0][0]
        failures = [
            requests.Timeout("read timed out"), requests.ConnectionError("offline"),
            response("unavailable", 503), response("<html>maintenance</html>"),
            response("bad,columns\n1,2\n"), response(HEADER + "123,150\n"),
            response("<VOTABLE>broken"),
        ]
        for failure in failures:
            with self.subTest(failure=str(failure)):
                candidates, post, output, warnings = self.query([failure, response(votable())])
                actual = candidates[0]
                self.assertEqual(actual.source_id, expected.source_id)
                self.assertEqual(actual.rp_mag, expected.rp_mag)
                self.assertEqual(actual.bp_mag, expected.bp_mag)
                self.assertEqual(actual.propagated, expected.propagated)
                self.assertAlmostEqual(actual.coord.separation(expected.coord).arcsec, 0.)
                self.assertEqual(post.call_count, 2)
                self.assertEqual(post.call_args.args[0], chart.MAST_GAIA_TAP_SYNC_URL)
                query = post.call_args.kwargs["data"]["QUERY"]
                self.assertIn("FROM dbo.gaia_source", query)
                self.assertIn("phot_rp_mean_mag > 16.00000000", query)
                self.assertIn("0.0250000000", query)
                self.assertIn("MAST/STScI", output)
                self.assertIn("trying MAST", warnings)

    # Empty responses indicate no matches rather than a service outage.
    def test_empty_csv_and_votable(self):
        for body in (HEADER, votable(empty=True)):
            with self.subTest(format=body[:20]):
                candidates, post, _, _ = self.query([response(body)])
                self.assertEqual(candidates, [])
                self.assertEqual(post.call_count, 1)

    # Missing proper motions retain the catalog position and masked photometry.
    def test_masked_votable_values(self):
        candidates, _, _, _ = self.query([requests.Timeout(), response(votable(masked=True))])
        self.assertFalse(candidates[0].propagated)
        self.assertIsNone(candidates[0].bp_mag)
        self.assertAlmostEqual(candidates[0].coord.ra.deg, 150.001)

    # Surface both backend errors when neither can return a valid result.
    def test_both_fail(self):
        with self.assertRaisesRegex(RuntimeError, "Both Gaia DR3 services failed.*ESA:.*MAST/STScI:"):
            self.query([requests.Timeout("ESA timeout"), response("offline", 503)])

    # Reject TAP errors or overflow even if accompanied by a syntactically valid table.
    def test_votable_query_status(self):
        for status in ("ERROR", "OVERFLOW"):
            xml = votable().replace(
                b"</RESOURCE>",
                f'<INFO name="QUERY_STATUS" value="{status}">incomplete</INFO></RESOURCE>'.encode(),
            )
            with self.subTest(status=status):
                with self.assertRaisesRegex(RuntimeError, status):
                    self.query([requests.Timeout(), response(xml)])

    # Render a real chart from synthetic FITS after a simulated ESA outage.
    def test_local_faint_chart_with_mast_fallback(self):
        with tempfile.TemporaryDirectory() as directory:
            image = Path(directory) / "field.fits"
            output = Path(directory) / "chart.png"
            wcs = chart.make_pa_wcs(150., 20., 120, 1., 0.)
            fits.PrimaryHDU(
                np.arange(120 * 120, dtype="float32").reshape(120, 120),
                header=wcs.to_header(),
            ).writeto(image)

            # Use the actual CLI, Gaia parser, propagation, and chart renderer.
            argv = [
                "get_finding_chart.py", "-c", "150", "20", "-i", str(image),
                "--faint", "-f", "1", "-o", str(output),
            ]
            stdout = io.StringIO()
            with patch.object(sys, "argv", argv), patch.object(
                chart.requests, "post", side_effect=[requests.Timeout(), response(votable())]
            ), contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(chart.main(), 0)
            self.assertTrue(output.read_bytes().startswith(b"\x89PNG"))
            self.assertIn("Gaia DR3 backend: MAST/STScI", stdout.getvalue())
            self.assertIn(str(SOURCE_ID), stdout.getvalue())

    # Keep sorting and the propagated-radius exclusion shared between backends.
    def test_radius_and_candidate_order(self):
        body = HEADER + (
            f"{SOURCE_ID},150.002,20,0,0,2016,18,18.5,17.5\n"
            "123,150.001,20,0,0,2016,18,18.5,17.5\n"
            "456,150.026,20,1000,0,2016,18,18.5,17.5\n"
        )
        candidates, _, _, _ = self.query([response(body)])
        self.assertEqual([star.source_id for star in candidates], ["123", str(SOURCE_ID)])


# Allow direct execution with the standard unittest runner.
if __name__ == "__main__":
    unittest.main()
