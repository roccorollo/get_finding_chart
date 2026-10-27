"""Regression checks for deprecated FITS WCS matrix keywords."""

import warnings
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from get_finding_chart import normalize_deprecated_pc_keywords, read_fits_image


class WcsKeywordTests(unittest.TestCase):
    def test_legacy_pc_keywords_are_renamed_without_changing_values(self):
        header = fits.Header()
        header["PC001001"] = (-1.0, "legacy matrix entry")
        header["PC001002"] = 0.25
        normalized = normalize_deprecated_pc_keywords(header)

        self.assertEqual(normalized["PC1_1"], -1.0)
        self.assertEqual(normalized.comments["PC1_1"], "legacy matrix entry")
        self.assertEqual(normalized["PC1_2"], 0.25)
        self.assertNotIn("PC001001", normalized)
        self.assertIn("PC001001", header)  # The downloaded header is untouched.

        header["PC1_1"] = 2.0
        normalized = normalize_deprecated_pc_keywords(header)
        self.assertEqual(normalized["PC1_1"], 2.0)
        self.assertNotIn("PC001001", normalized)

    def test_legacy_fits_has_same_wcs_without_pc_warnings(self):
        header = fits.Header()
        header["CTYPE1"] = "RA---TAN"
        header["CTYPE2"] = "DEC--TAN"
        header["CRVAL1"] = 256.398
        header["CRVAL2"] = -23.456
        header["CRPIX1"] = 10.5
        header["CRPIX2"] = 10.5
        header["CDELT1"] = 1 / 3600
        header["CDELT2"] = 1 / 3600
        header["PC001001"] = -1.0
        header["PC001002"] = 0.0
        header["PC002001"] = 0.0
        header["PC002002"] = 1.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            original_wcs = WCS(header).celestial

        with TemporaryDirectory() as directory:
            path = Path(directory) / "legacy_pc.fits"
            fits.PrimaryHDU(np.ones((20, 20)), header=header).writeto(path)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                data, normalized_wcs = read_fits_image(path)

        self.assertEqual(data.shape, (20, 20))
        self.assertFalse(
            any("PCi_ja keyword is deprecated" in str(item.message) for item in caught)
        )
        x = np.array([0, 10, 19])
        y = np.array([0, 10, 19])
        for before, after in zip(
            original_wcs.pixel_to_world_values(x, y),
            normalized_wcs.pixel_to_world_values(x, y),
        ):
            np.testing.assert_allclose(after, before, rtol=0, atol=1e-10)


if __name__ == "__main__":
    unittest.main()
