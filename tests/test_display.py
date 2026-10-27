"""Checks for automatic sky contrast and desktop viewer dispatch."""

import contextlib
import io
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

import get_finding_chart as chart


# Exercise visible sky brightness, explicit cuts, and platform viewer commands.
class DisplayTests(unittest.TestCase):
    # Auto mode should reach a lighter sky while keeping bright stars black.
    def test_auto_scaling_lightens_sky(self):
        rng = np.random.default_rng(42)
        pixels = rng.normal(100.0, 10.0, 10000)
        pixels[:20] = 400.0
        low, high = chart.compute_intensity_limits(
            pixels, "auto", 1.0, 99.5, None, None
        )
        sky = np.median(pixels[20:])
        normalized = (sky - low) / (high - low)
        brightness = 1.0 - float(
            chart.get_stretch("asinh")(np.array([normalized]))[0]
        )
        self.assertGreaterEqual(brightness, 0.74)
        self.assertLessEqual(brightness, 0.78)
        self.assertLess(high, pixels[0])

    # Absolute user cuts and percentile mode keep their requested limits.
    def test_manual_and_percentile_cuts_are_preserved(self):
        pixels = np.arange(100.0)
        self.assertEqual(
            chart.compute_intensity_limits(pixels, "auto", 1, 99, 10, 80),
            (10.0, 80.0),
        )
        low, high = chart.compute_intensity_limits(
            pixels, "percentile", 10, 90, None, None
        )
        self.assertAlmostEqual(low, 9.9)
        self.assertAlmostEqual(high, 89.1)

    # Linux and macOS use their respective default desktop file handlers.
    def test_graphics_uses_platform_default_opener(self):
        output = Path("chart.png")
        for platform, command in (("linux", "xdg-open"), ("darwin", "open")):
            with self.subTest(platform=platform), patch.object(
                sys, "platform", platform
            ), patch.object(chart.shutil, "which", return_value=f"/usr/bin/{command}"), patch.object(
                chart.subprocess, "Popen"
            ) as launch:
                chart.open_chart(output)
                self.assertEqual(launch.call_args.args[0], [command, str(output.resolve())])

    # Missing desktop helpers should warn without losing the saved chart.
    def test_graphics_without_opener_warns(self):
        errors = io.StringIO()
        with patch.object(sys, "platform", "linux"), patch.object(
            chart.shutil, "which", return_value=None
        ), contextlib.redirect_stderr(errors):
            chart.open_chart(Path("chart.png"))
        self.assertIn("No supported system opener", errors.getvalue())


if __name__ == "__main__":
    unittest.main()
