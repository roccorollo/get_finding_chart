"""Regression checks for acquisition moves toward the main target."""

import contextlib
import io
import unittest

from astropy import units as u
from astropy.coordinates import SkyCoord

from get_finding_chart import (
    FaintSelection,
    GaiaCandidate,
    relative_geometry_values,
    report_gaia_candidates,
    report_relative_geometry,
    source_to_target_offset_text,
)


TARGET = SkyCoord(
    "17:05:35.520", "-23:27:21.60", unit=(u.hourangle, u.deg), frame="icrs"
)
S = SkyCoord(
    "17:05:34.0", "-23:27:31.60", unit=(u.hourangle, u.deg), frame="icrs"
)


class OffsetDirectionTests(unittest.TestCase):
    def test_two_source_move_points_east_and_north_to_target(self):
        separation, east, north, pa, opposite_pa = relative_geometry_values(
            S, TARGET
        )
        self.assertAlmostEqual(separation, 23.1833, places=3)
        self.assertAlmostEqual(east, 20.9159, places=3)
        self.assertAlmostEqual(north, 9.9995, places=3)
        self.assertAlmostEqual(pa, 64.4484, places=3)
        self.assertAlmostEqual(opposite_pa, 244.4484, places=3)
        self.assertEqual(
            source_to_target_offset_text(S, TARGET, "S"),
            "S -> T offsets\nDeltaRA=+20.92 arcsec\nDeltaDec=+10.00 arcsec",
        )

        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            report_relative_geometry(
                [(S.ra.deg, S.dec.deg), (TARGET.ra.deg, TARGET.dec.deg)]
            )
        self.assertIn("DeltaRA=+20.92 arcsec; DeltaDec=+10.00 arcsec", output.getvalue())
        self.assertIn("PA(S->T)=64.45 deg", output.getvalue())

    def test_gaia_acquisition_move_points_east_and_south_to_target(self):
        # Reproduce the G2 displacement shown on the existing chart.
        g2 = TARGET.spherical_offsets_by(-15.79 * u.arcsec, 4.87 * u.arcsec)
        self.assertEqual(
            source_to_target_offset_text(g2, TARGET, "G2"),
            "G2 -> T offsets\nDeltaRA=+15.79 arcsec\nDeltaDec=-4.87 arcsec",
        )
        candidates = [
            GaiaCandidate(
                source_id="1",
                g_mag=15.0,
                bp_mag=15.5,
                rp_mag=14.5,
                coord=TARGET.spherical_offsets_by(-5 * u.arcsec, 0 * u.arcsec),
                propagated=True,
                separation_arcsec=5.0,
            ),
            GaiaCandidate(
                source_id="2",
                g_mag=16.0,
                bp_mag=16.5,
                rp_mag=15.5,
                coord=g2,
                propagated=True,
                separation_arcsec=TARGET.separation(g2).arcsec,
            ),
        ]
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            report_gaia_candidates(TARGET, candidates, FaintSelection(2, 14, 17), 2026.0)
        g2_report = output.getvalue().split("G2 [selected]:", 1)[1]
        self.assertIn("DeltaRA=+15.79 arcsec; DeltaDec=-4.87 arcsec", g2_report)
        self.assertIn("PA(G2->T)=107.14 deg", g2_report)


if __name__ == "__main__":
    unittest.main()
