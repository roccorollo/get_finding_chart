"""Offline regressions for source order, chart centering, and target visibility."""

import contextlib
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from astropy import units as u
from astropy.coordinates import SkyCoord

import get_finding_chart as chart


# Check input order and the two-line heading under a constrained chart width.
class SourceOrderTests(unittest.TestCase):
    # Preserve the target label even when intermediate sources are truncated.
    def test_heading_keeps_final_target(self):
        coords = [SkyCoord(150 + i * 0.001, 20, unit=u.deg) for i in range(10)]
        with patch.object(chart, "text_fits_axis", side_effect=lambda fig, ax, text: len(text) <= 75):
            lines = chart.build_info_lines(None, None, coords, ["x"] * 10, 1.0, 240.0)

        self.assertEqual(len(lines), 2)
        self.assertIn("S:", lines[0])
        self.assertIn("T: RA", lines[1])
        self.assertIn("[...]", lines[1])

    # Center a normal downloaded chart on S while retaining T for offsets.
    def test_normal_mode_centers_on_first_source(self):
        observed, output = self.run_chart(faint=False)
        self.assertEqual(observed["retrieval_center"], (150.0, 20.0))
        self.assertEqual(observed["render_center"], (150.0, 20.0))
        self.assertEqual(observed["coords"], [(150.0, 20.0), (150.004, 20.004)])
        self.assertIn("PA(S->T)=", output)

    # Let the selected Gaia star take precedence over S in faint mode.
    def test_faint_mode_centers_on_selected_gaia_star(self):
        observed, output = self.run_chart(faint=True)
        self.assertEqual(observed["retrieval_center"], (150.0041, 20.0041))
        self.assertEqual(observed["render_center"], (150.0041, 20.0041))
        self.assertEqual(observed["gaia_target"], (150.004, 20.004))
        self.assertEqual(observed["markers"][-1], "slit")
        self.assertIn("PA(S->T)=", output)

    # Run the CLI with local substitutes for survey and Gaia services.
    def run_chart(self, *, faint):
        observed = {}
        target = SkyCoord(150.004, 20.004, unit=u.deg)
        star = SkyCoord(150.0041, 20.0041, unit=u.deg)
        candidate = chart.GaiaCandidate(
            "123", 18.0, 18.1, 18.2, star, True, star.separation(target).arcsec
        )

        # Capture the center passed to the survey retrieval and chart renderer.
        def retrieve(*args):
            observed["retrieval_center"] = args[1:3]
            args[9].write_bytes(b"synthetic FITS")
            return chart.RetrievalInfo("synthetic")

        # Record the chart coordinates and publish a small stand-in image.
        def render(*args):
            observed["render_center"] = args[3]
            observed["coords"] = args[2]
            observed["markers"] = args[6]
            args[1].write_bytes(b"synthetic chart")

        # Ensure faint acquisition queries the final target T.
        def query(coord, selection, epoch):
            observed["gaia_target"] = (coord.ra.deg, coord.dec.deg)
            return [candidate]

        with tempfile.TemporaryDirectory() as directory:
            output_path = Path(directory) / "chart.png"
            argv = [
                "get_finding_chart.py", "-c", "150", "20", "150.004", "20.004",
                "-s", "legacy", "-o", str(output_path),
            ]
            if faint:
                argv.append("--faint")
            output = io.StringIO()
            with patch.object(sys, "argv", argv), patch.object(
                chart, "retrieve_image", side_effect=retrieve
            ), patch.object(chart, "make_finding_chart", side_effect=render), patch.object(
                chart, "query_gaia_candidates", side_effect=query
            ), contextlib.redirect_stdout(output):
                self.assertEqual(chart.main(), 0)
            self.assertEqual(output_path.read_bytes(), b"synthetic chart")

        return observed, output.getvalue()


# Allow direct execution through Python's unittest runner.
if __name__ == "__main__":
    unittest.main()
