"""Checks for UHS J-band survey selection and HiPS dispatch."""

import unittest
from pathlib import Path
from unittest.mock import patch

import get_finding_chart as chart


# Verify the UHS option uses the registered J-band image survey.
class UhsSurveyTests(unittest.TestCase):
    # The CLI accepts UHS and uses J as the default band.
    def test_uhs_is_a_cli_survey_with_j_default(self):
        args = chart.make_parser().parse_args(["-c", "10 20", "-s", "uhs"])
        self.assertEqual(args.survey, "uhs")
        self.assertEqual(chart.SURVEYS[args.survey].default_band, "J")

    # The resolver accepts either case of J and rejects unavailable bands.
    def test_uhs_band_resolution(self):
        expected = ("j", "wfau.roe.ac.uk/P/UHSDR1/J")
        self.assertEqual(chart.resolve_primary_hips_band("uhs", "J"), expected)
        self.assertEqual(chart.resolve_primary_hips_band("uhs", "j"), expected)
        with self.assertRaisesRegex(ValueError, "band must be one of: J"):
            chart.resolve_primary_hips_band("uhs", "H")

    # Image retrieval passes the UHS HiPS ID and marks its FITS as a HiPS product.
    def test_uhs_retrieval_uses_hips2fits(self):
        output = Path("uhs.fits")
        with patch.object(chart, "get_hips2fits") as download:
            result = chart.retrieve_image(
                "uhs", 10.0, 20.0, 100, 1.0, 1.0, 0.403, 0.0, "J", output, "ls-dr10"
            )
        download.assert_called_once_with(
            10.0, 20.0, 100, 1.0, "wfau.roe.ac.uk/P/UHSDR1/J", output
        )
        self.assertTrue(result.from_hips2fits)
        self.assertIn("UKIRT Hemisphere Survey DR1", result.backend)
