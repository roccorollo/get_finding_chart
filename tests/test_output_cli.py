"""Offline checks for finding-chart output paths and replacement rules."""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import get_finding_chart as chart


# Replace a survey download with a small staged file for CLI routing checks.
def fake_retrieval(*args):
    args[9].write_bytes(b"survey FITS")
    return chart.RetrievalInfo("mock archive")


# Replace rendering with a file write that records the requested chart format.
def fake_render(*args):
    args[1].write_bytes(b"chart " + args[1].suffix.lower().encode())


# Check output format inference, collision refusal, and optional FITS publication.
class OutputCliTests(unittest.TestCase):
    # Run the chart CLI with offline download and rendering substitutes.
    def run_chart(self, argv):
        with patch.object(chart, "retrieve_image", side_effect=fake_retrieval) as retrieval:
            with patch.object(chart, "make_finding_chart", side_effect=fake_render):
                with patch.object(sys, "argv", ["get_finding_chart.py", "-c", "150", "20",
                                                "-s", "legacy", *argv]):
                    with contextlib.redirect_stdout(io.StringIO()):
                        return chart.main(), retrieval

    # Infer PDF from the explicit filename and publish a native FITS when asked.
    def test_output_suffix_and_fits(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "chart.pdf"
            fits_output = Path(directory) / "field.fits"
            result, _ = self.run_chart(["-o", str(output), "--fits", str(fits_output)])
            self.assertEqual(result, 0)
            self.assertEqual(output.read_bytes(), b"chart .pdf")
            self.assertEqual(fits_output.read_bytes(), b"survey FITS")

    # Refuse both output collisions before requesting data, then allow replacement.
    def test_existing_outputs_need_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "chart.png"
            fits_output = Path(directory) / "field.fits"
            output.write_bytes(b"old chart")
            fits_output.write_bytes(b"old FITS")
            argv = ["-o", str(output), "--fits", str(fits_output)]
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as result:
                    self.run_chart(argv)
            self.assertEqual(result.exception.code, 2)
            self.assertEqual(output.read_bytes(), b"old chart")
            self.assertEqual(fits_output.read_bytes(), b"old FITS")
            result, _ = self.run_chart([*argv, "--overwrite"])
            self.assertEqual(result, 0)
            self.assertEqual(output.read_bytes(), b"chart .png")
            self.assertEqual(fits_output.read_bytes(), b"survey FITS")

    # Reject a chart format that conflicts with the explicit filename suffix.
    def test_extension_conflict(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "chart.pdf"
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as result:
                    self.run_chart(["-o", str(output), "-e", "png"])
            self.assertEqual(result.exception.code, 2)
            self.assertFalse(output.exists())

    # A rendering error must leave an existing chart untouched when replacing.
    def test_failed_render_preserves_chart(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "chart.png"
            output.write_bytes(b"old chart")
            with patch.object(chart, "retrieve_image", side_effect=fake_retrieval):
                with patch.object(chart, "make_finding_chart", side_effect=OSError("render failed")):
                    argv = ["get_finding_chart.py", "-c", "150", "20", "-s", "legacy",
                            "-o", str(output), "--overwrite"]
                    with patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()):
                        with contextlib.redirect_stderr(io.StringIO()):
                            self.assertEqual(chart.main(), 1)
            self.assertEqual(output.read_bytes(), b"old chart")


# Allow this module to run directly through unittest.
if __name__ == "__main__":
    unittest.main()
