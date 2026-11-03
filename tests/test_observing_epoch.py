"""Offline checks for observing dates, current-time defaults, and explicit epochs."""

import argparse
import contextlib
import io
import unittest
from unittest.mock import patch

from astropy.time import Time

import get_finding_chart as chart


# Check date validation and ensure the default clock is resolved at invocation time.
class ObservingEpochTests(unittest.TestCase):
    # A fresh call must follow the clock, even across a year boundary.
    def test_default_tracks_current_time(self):
        times = [Time("2027-12-31T23:59:59"), Time("2028-01-01T00:00:01")]
        with patch.object(chart.Time, "now", side_effect=times) as clock:
            epochs = [chart.resolve_gaia_epoch(None, None) for _ in times]
        self.assertEqual(clock.call_count, 2)
        self.assertEqual(epochs, [instant.jyear for instant in times])
        self.assertLess(epochs[0], epochs[1])

    # Midnight precedes the J2000 noon origin by half of one Julian day.
    def test_calendar_date_is_midnight_utc(self):
        epoch = chart.resolve_gaia_epoch(None, "2000-01-01")
        self.assertAlmostEqual(epoch, 2000.0 - 0.5 / 365.25, places=10)

    # Accept leap days only in leap years and require unambiguous calendar dates.
    def test_date_validation(self):
        self.assertIsInstance(chart.resolve_gaia_epoch(None, "2028-02-29"), float)
        for value in ("2027-02-29", "2028-02-30", "2028-13-01", "2028-2-01",
                      "01/02/2028", "2028-02-01T12:00:00", "today", ""):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    chart.resolve_gaia_epoch(None, value)

    # Preserve numeric epoch overrides and their established validity constraints.
    def test_explicit_epoch_validation(self):
        self.assertEqual(chart.resolve_gaia_epoch(2000.0, None), 2000.0)
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    chart.resolve_gaia_epoch(value, None)

    # Reject conflicting controls before chart retrieval or rendering can begin.
    def test_mutually_exclusive_cli_controls(self):
        errors = io.StringIO()
        with contextlib.redirect_stderr(errors):
            with self.assertRaises(SystemExit) as failure:
                chart.make_parser().parse_args(chart.normalize_coordinate_argv([
                    "-c", "150", "20", "--date", "2028-02-29", "--epoch", "2028.5"
                ]))
        self.assertEqual(failure.exception.code, 2)
        self.assertIn("not allowed with argument", errors.getvalue())


# Allow direct execution with the standard unittest runner.
if __name__ == "__main__":
    unittest.main()
