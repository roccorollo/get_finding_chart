"""Checks for local FITS input and the revised short options."""

import contextlib
import io
import os
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import numpy as np
from astropy.io import fits

import get_finding_chart as chart


class LocalImageTests(unittest.TestCase):
    def test_image_option_renders_local_fits_without_survey_download(self):
        with TemporaryDirectory() as directory:
            workdir = Path(directory)
            image_path = workdir / "input.fits"
            y, x = np.mgrid[:100, :100]
            data = (x + y + np.sin(x / 3)).astype("float32")
            wcs = chart.make_pa_wcs(256.398, -23.456, 100, 1.0, 0.0)
            fits.PrimaryHDU(data, header=wcs.to_header()).writeto(image_path)

            previous_dir = Path.cwd()
            output = io.StringIO()
            try:
                os.chdir(workdir)
                argv = [
                    "get_finding_chart.py", "-c", "256.398", "-23.456",
                    "-f", "0.5", "-i", str(image_path), "-n", "local",
                ]
                with patch.object(sys, "argv", argv), patch.object(
                    chart, "retrieve_image", side_effect=AssertionError("network used")
                ), contextlib.redirect_stdout(output):
                    self.assertEqual(chart.main(), 0)
            finally:
                os.chdir(previous_dir)

            result = workdir / "local_image_256p398_-23p456_0p5arcmin.png"
            self.assertTrue(result.is_file())
            self.assertGreater(result.stat().st_size, 0)
            self.assertIn("Input FITS image:", output.getvalue())
            self.assertNotIn("Survey download:", output.getvalue())

    def test_image_and_downloaded_fits_output_are_incompatible(self):
        argv = [
            "get_finding_chart.py", "-c", "256.398", "-23.456",
            "-i", "input.fits", "--fits", "copy.fits",
        ]
        with patch.object(sys, "argv", argv), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as failure:
                chart.main()
        self.assertEqual(failure.exception.code, 2)

    def test_reassigned_flags_have_no_parser_collision(self):
        parser = chart.make_parser()
        args = parser.parse_args(["-c", "256.398 -23.456", "-n", "target", "-i", "input.fits"])
        self.assertEqual(args.name, "target")
        self.assertEqual(args.image, "input.fits")
        self.assertNotIn("--id", parser.format_help())

    # Confirm graphics launches only after a local chart was saved.
    def test_graphics_opens_successful_local_chart(self):
        with TemporaryDirectory() as directory:
            workdir = Path(directory)
            image_path = workdir / "input.fits"
            wcs = chart.make_pa_wcs(256.398, -23.456, 80, 1.0, 0.0)
            fits.PrimaryHDU(
                np.arange(80 * 80, dtype="float32").reshape(80, 80),
                header=wcs.to_header(),
            ).writeto(image_path)

            previous_dir = Path.cwd()
            try:
                os.chdir(workdir)
                argv = [
                    "get_finding_chart.py", "-c", "256.398", "-23.456",
                    "-f", "0.5", "-i", str(image_path), "-n", "local", "-g",
                ]
                with patch.object(sys, "argv", argv), patch.object(
                    chart, "open_chart"
                ) as opener, contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(chart.main(), 0)
                result = workdir / "local_image_256p398_-23p456_0p5arcmin.png"
                self.assertTrue(result.is_file())
                opener.assert_called_once_with(Path(result.name))
            finally:
                os.chdir(previous_dir)

    # Failed rendering must not launch a viewer for an incomplete chart.
    def test_graphics_does_not_open_failed_chart(self):
        with TemporaryDirectory() as directory:
            image_path = Path(directory) / "empty.fits"
            wcs = chart.make_pa_wcs(256.398, -23.456, 80, 1.0, 0.0)
            fits.PrimaryHDU(
                np.full((80, 80), np.nan, dtype="float32"),
                header=wcs.to_header(),
            ).writeto(image_path)
            argv = [
                "get_finding_chart.py", "-c", "256.398", "-23.456",
                "-f", "0.5", "-i", str(image_path), "-g",
            ]
            with patch.object(sys, "argv", argv), patch.object(
                chart, "open_chart"
            ) as opener, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(chart.main(), 1)
            opener.assert_not_called()


if __name__ == "__main__":
    unittest.main()
