#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Create annotated astronomical finding charts from public sky surveys.

Supports Pan-STARRS1, DESI Legacy Surveys, SkyMapper DR4, SDSS DR9,
GALEX GR6/7, DES DR2, VISTA VIKING, VISTA Hemisphere Survey (VHS),
DSS2, 2MASS, and AllWISE. Includes an optional Gaia DR3 faint-source
acquisition mode and inverted grayscale display.
See README.md for usage examples.
"""

from __future__ import annotations

__author__ = "Andrea Rossi"
__credits__ = "Development assistance: OpenAI ChatGPT"
__license__ = "MIT"
__version__ = "2.0.5"

import argparse
import csv
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from io import BytesIO, StringIO
from pathlib import Path
from urllib.parse import urljoin

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import Angle, Latitude, Longitude, SkyCoord
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.stats import sigma_clipped_stats
from astropy.visualization import (
    AsinhStretch,
    ImageNormalize,
    LinearStretch,
    LogStretch,
    SqrtStretch,
    ZScaleInterval,
)
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
try:
    from astroquery.hips2fits import hips2fits
except ImportError:  # Allow older astroquery installs to use the fallback backends.
    hips2fits = None

from astroquery.skyview import SkyView
from matplotlib.colors import is_color_like
from matplotlib.patches import Circle, Polygon, Rectangle
from scipy.ndimage import binary_propagation, map_coordinates


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

DEFAULT_TIMEOUT = 120
DOWNLOAD_CHUNK_SIZE = 1024 * 1024
REPROJECT_BLOCK_SIZE = 256
MAX_SCALING_SAMPLE = 1_000_000
MAX_SOURCES = 10
DEFAULT_FOV_ARCMIN = 4.0
DEFAULT_SLIT_LENGTH_ARCSEC = 240.0
DEFAULT_SLIT_WIDTH_ARCSEC = 1.0
MIN_SLIT_WIDTH_ARCSEC = 0.1
MAX_SLIT_WIDTH_ARCSEC = 60.0
ALLOWED_MARKERS = ("cross", "circle", "diamond", "x", "slit")

# Base size of compact markers. User --marker-size values are multiplicative
# factors relative to this default appearance. Slit length/width remain physical.
MARKER_SCALE = 0.50

# Okabe-Ito color-blind-friendly palette. The palette is cycled only when
# more sources are requested than distinct default colors are available.
COLORBLIND_PALETTE = (
    "#D55E00",  # vermillion
    "#0072B2",  # blue
    "#009E73",  # bluish green
    "#CC79A7",  # reddish purple
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
)
ACQUISITION_COLOR = COLORBLIND_PALETTE[1]
COMPASS_COLOR = "#006400"  # dark green, intentionally independent of marker palette
DEFAULT_MARKER_SIZE = 1.0

SKYMAPPER_SIAP_URL = "https://api.skymapper.nci.org.au/public/siap/dr4/query"
SKYMAPPER_BASE_URL = "https://api.skymapper.nci.org.au/"
SKYMAPPER_MAX_CUTOUT_ARCMIN = 10.0
GAIA_TAP_SYNC_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
ESO_TAP_SYNC_URL = "https://archive.eso.org/tap_obs/sync"
ESO_SODA_SYNC_URL = "https://dataportal.eso.org/dataPortal/soda/sync"

GAIA_DR3_TABLE = "gaiadr3.gaia_source"
GAIA_SEARCH_RADIUS_ARCSEC = 90.0
DEFAULT_GAIA_RP_MAX = 19.0
DEFAULT_GAIA_EPOCH = 2026.0
MAX_GAIA_TERMINAL_SOURCES = 5

# Native PS1/2MASS cutouts are preferred. If more than this fraction of the
# final displayed field consists of edge-connected non-finite pixels, the code
# retries with the mosaicked CDS HiPS2FITS service.  Isolated NaNs away from the borders do
# not trigger a fallback.
MAX_EDGE_MISSING_FRACTION = 0.01
MAX_CENTER_MISSING_FRACTION = 0.50
CENTER_CHECK_FRACTION = 0.05

HIPS_PS1 = {
    band: f"CDS/P/PanSTARRS/DR1/{band}" for band in ("g", "r", "i", "z", "y")
}
HIPS_2MASS = {band: f"CDS/P/2MASS/{band}" for band in ("J", "H", "K")}

# Surveys retrieved directly through CDS HiPS2FITS. Keys are normalized
# lower-case CLI band names; values are the corresponding CDS HiPS IDs.
HIPS_PRIMARY = {
    "sdss": {
        band: f"CDS/P/SDSS9/{band}" for band in ("u", "g", "r", "i", "z")
    },
    "galex": {
        "fuv": "CDS/P/GALEXGR6_7/FUV",
        "nuv": "CDS/P/GALEXGR6_7/NUV",
    },
    "des": {
        "g": "CDS/P/DES-DR2/g",
        "r": "CDS/P/DES-DR2/r",
        "i": "CDS/P/DES-DR2/i",
        "z": "CDS/P/DES-DR2/z",
        "y": "CDS/P/DES-DR2/Y",
    },
    "viking": {
        "z": "CDS/P/VISTA/VIKING/z",
        "y": "CDS/P/VISTA/VIKING/Y",
        "j": "CDS/P/VISTA/VIKING/J",
        "h": "CDS/P/VISTA/VIKING/H",
        "k": "CDS/P/VISTA/VIKING/K",
    },
    "allwise": {
        "w1": "CDS/P/allWISE/W1",
        "w2": "CDS/P/allWISE/W2",
        "w3": "CDS/P/allWISE/W3",
        "w4": "CDS/P/allWISE/W4",
    },
}

HIPS_BAND_ALIASES = {
    "viking": {"ks": "k", "k_s": "k"},
    "allwise": {
        "1": "w1",
        "3.4": "w1",
        "2": "w2",
        "4.6": "w2",
        "3": "w3",
        "12": "w3",
        "4": "w4",
        "22": "w4",
    },
}

HIPS_BAND_HELP = {
    "sdss": "u g r i z",
    "galex": "FUV NUV",
    "des": "g r i z Y",
    "viking": "z Y J H K",
    "allwise": "w1 w2 w3 w4",
}


@dataclass(frozen=True)
class SurveySpec:
    label: str
    tag: str
    default_band: str
    pixel_scale: float  # arcsec/pixel used for requested output sampling


@dataclass(frozen=True)
class RetrievalInfo:
    backend: str
    from_hips2fits: bool = False


@dataclass(frozen=True)
class CoverageInfo:
    finite_fraction: float
    edge_missing_fraction: float
    center_missing_fraction: float
    acceptable: bool


@dataclass(frozen=True)
class FaintSelection:
    rank: int
    mag_min: float | None
    mag_max: float


@dataclass(frozen=True)
class GaiaCandidate:
    source_id: str
    g_mag: float | None
    bp_mag: float | None
    rp_mag: float
    coord: SkyCoord
    propagated: bool
    separation_arcsec: float


SURVEYS = {
    "ps1": SurveySpec("Pan-STARRS1", "ps1", "r", 0.25),
    "legacy": SurveySpec("DESI Legacy Survey DR10", "ls10", "r", 0.262),
    "skymapper": SurveySpec("SkyMapper DR4", "smdr4", "r", 0.50),
    "sdss": SurveySpec("SDSS DR9", "sdss9", "r", 0.396),
    "galex": SurveySpec("GALEX GR6/7", "galex", "nuv", 1.5),
    "des": SurveySpec("DES DR2", "desdr2", "r", 0.263),
    "viking": SurveySpec("VISTA VIKING", "viking", "K", 0.339),
    "vhs": SurveySpec("VISTA Hemisphere Survey", "vhs", "K", 0.339),
    "dss2": SurveySpec("DSS2", "dss2", "red", 1.0),
    "2mass": SurveySpec("2MASS", "2mass", "J", 1.0),
    "allwise": SurveySpec("AllWISE", "allwise", "w1", 1.375),
}

SURVEY_ALIASES = {
    "ls": "legacy",
    "sm": "skymapper",
    "vik": "viking",
    "aw": "allwise",
    "panstarrs": "ps1",
}

# Keep the --help survey list in a deliberate human-friendly order:
# short names first, followed by the long/canonical names.
SURVEY_CHOICES = (
    "ps1", "ls", "sm", "sdss", "galex", "des", "vik", "vhs", "dss2", "2mass", "aw",
    "panstarrs", "legacy", "skymapper", "sdss", "galex", "des", "viking", "vhs", "dss2", "2mass", "allwise",
)


# -----------------------------------------------------------------------------
# Command-line parsing helpers
# -----------------------------------------------------------------------------


def looks_like_option(token: str) -> bool:
    """Return True for CLI options, but not for negative coordinate values."""
    return token.startswith("--") or bool(re.match(r"^-[A-Za-z]", token))


def normalize_coordinate_argv(argv: list[str]) -> list[str]:
    """Combine all coordinate tokens after -c/--coordinates into one value.

    This prevents argparse from interpreting negative sexagesimal declinations
    such as ``-23:27:21.60`` as command-line options.
    """
    args = list(argv)

    for i, token in enumerate(args):
        if token not in {"-c", "--coordinates"}:
            continue

        j = i + 1
        values: list[str] = []
        while j < len(args) and not looks_like_option(args[j]):
            values.append(args[j])
            j += 1

        if values:
            args[i + 1:j] = [" ".join(values)]
        break

    return args


def parse_coordinate_pair(ra_text: str, dec_text: str) -> tuple[float, float]:
    """Parse one RA/Dec pair and return decimal degrees."""
    try:
        ra = float(ra_text)
        dec = float(dec_text)
    except ValueError:
        try:
            coord = SkyCoord(
                ra_text,
                dec_text,
                unit=(u.hourangle, u.deg),
                frame="icrs",
            )
        except Exception as exc:
            raise argparse.ArgumentTypeError(
                f"Could not parse coordinates: {ra_text} {dec_text}"
            ) from exc
    else:
        try:
            coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
        except Exception as exc:
            raise argparse.ArgumentTypeError(
                f"Invalid decimal coordinates: {ra_text} {dec_text}"
            ) from exc

    return float(coord.ra.deg), float(coord.dec.deg)


def parse_coordinates(value: str) -> list[tuple[float, float]]:
    """Parse up to MAX_SOURCES RA/Dec pairs from one normalized string."""
    parts = str(value).strip().replace(",", " ").split()

    if len(parts) < 2 or len(parts) % 2:
        raise argparse.ArgumentTypeError(
            "Coordinates require one or more RA Dec pairs, e.g. "
            "-c 17:05:35.520 -23:27:21.60 "
            "17:05:36.100 -23:27:10.0"
        )

    n_sources = len(parts) // 2
    if n_sources > MAX_SOURCES:
        raise argparse.ArgumentTypeError(
            f"A maximum of {MAX_SOURCES} sources is supported; got {n_sources}."
        )

    return [
        parse_coordinate_pair(parts[i], parts[i + 1])
        for i in range(0, len(parts), 2)
    ]


def parse_fov(value: str) -> float:
    """Parse the square FOV side in arcminutes."""
    text = str(value).strip().lower()
    match = re.fullmatch(r"([0-9]+(?:\.[0-9]+)?)\s*(?:arcmin)?", text)
    if not match:
        raise argparse.ArgumentTypeError(
            "FOV must be a positive number in arcminutes, e.g. -f 5 or -f 3.7."
        )

    fov = float(match.group(1))
    if fov <= 0:
        raise argparse.ArgumentTypeError("FOV must be positive.")
    return fov


def parse_angle(value: str) -> float:
    """Parse a finite position angle in degrees."""
    try:
        angle = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Position angle must be a number in degrees."
        ) from exc

    if not np.isfinite(angle):
        raise argparse.ArgumentTypeError("Position angle must be finite.")
    return angle


def parse_slit_length(value: str) -> float:
    """Parse a positive slit length in arcsec."""
    try:
        length = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Slit length must be a number in arcseconds."
        ) from exc

    if not np.isfinite(length) or length <= 0.0:
        raise argparse.ArgumentTypeError(
            "Slit length must be a positive number in arcseconds."
        )
    return length


def parse_slit_width(value: str) -> float:
    """Parse a slit width in arcsec, constrained to a practical range."""
    try:
        width = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Slit width must be a number in arcseconds."
        ) from exc

    if not np.isfinite(width) or not MIN_SLIT_WIDTH_ARCSEC <= width <= MAX_SLIT_WIDTH_ARCSEC:
        raise argparse.ArgumentTypeError(
            f"Slit width must be between {MIN_SLIT_WIDTH_ARCSEC:g} and "
            f"{MAX_SLIT_WIDTH_ARCSEC:g} arcsec."
        )
    return width


def canonical_survey(name: str) -> str:
    """Resolve survey aliases to their canonical names."""
    return SURVEY_ALIASES.get(name.lower(), name.lower())


def resolve_primary_hips_band(survey: str, band: str) -> tuple[str, str]:
    """Return normalized band and CDS HiPS ID for a primary HiPS survey."""
    key = str(band).strip().lower()
    key = HIPS_BAND_ALIASES.get(survey, {}).get(key, key)
    mapping = HIPS_PRIMARY[survey]
    if key not in mapping:
        raise ValueError(
            f"{SURVEYS[survey].label} band must be one of: {HIPS_BAND_HELP[survey]}"
        )
    return key, mapping[key]


def resolve_markers(requested: list[str] | None, n_sources: int) -> list[str]:
    """Return exactly one marker per source, defaulting missing ones to x."""
    markers = list(requested or [])
    if len(markers) > n_sources:
        raise argparse.ArgumentTypeError(
            f"Received {len(markers)} markers but only {n_sources} sources."
        )
    markers.extend(["x"] * (n_sources - len(markers)))
    return markers


def resolve_colors(requested: list[str] | None, n_sources: int) -> list[str]:
    """Resolve source colors, using a color-blind-friendly default palette."""
    defaults = [
        COLORBLIND_PALETTE[index % len(COLORBLIND_PALETTE)]
        for index in range(n_sources)
    ]
    colors = list(requested or [])

    if not colors:
        return defaults
    if len(colors) > n_sources:
        raise argparse.ArgumentTypeError(
            f"Received {len(colors)} colors but only {n_sources} sources."
        )

    invalid = [color for color in colors if not is_color_like(color)]
    if invalid:
        raise argparse.ArgumentTypeError(
            "Invalid Matplotlib color(s): " + ", ".join(invalid)
        )

    if len(colors) == 1:
        return colors * n_sources

    return colors + defaults[len(colors):]


def resolve_marker_sizes(
    requested: list[float] | None,
    n_sources: int,
) -> list[float]:
    """Resolve marker-size factors, supporting one global or per-source values."""
    sizes = list(requested or [])

    if not sizes:
        return [DEFAULT_MARKER_SIZE] * n_sources
    if len(sizes) > n_sources:
        raise argparse.ArgumentTypeError(
            f"Received {len(sizes)} marker sizes but only {n_sources} sources."
        )
    if any((not np.isfinite(size)) or size <= 0 for size in sizes):
        raise argparse.ArgumentTypeError("--marker-size values must be positive.")

    if len(sizes) == 1:
        return sizes * n_sources

    return sizes + [DEFAULT_MARKER_SIZE] * (n_sources - len(sizes))


def parse_faint_selection(values: list[str] | None) -> FaintSelection | None:
    """Interpret -F/--faint optional values.

    Accepted forms are: --faint, --faint N, --faint MAG_MIN MAG_MAX, and
    --faint N MAG_MIN MAG_MAX.  Magnitude cuts always refer to Gaia G_RP.
    """
    if values is None:
        return None

    if len(values) == 0:
        return FaintSelection(1, None, DEFAULT_GAIA_RP_MAX)

    if len(values) == 1:
        try:
            rank_value = float(values[0])
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                "With one value, --faint expects a positive integer source rank."
            ) from exc
        rank = int(rank_value)
        if not np.isfinite(rank_value) or rank_value != rank or rank < 1:
            raise argparse.ArgumentTypeError(
                "With one value, --faint expects a positive integer source rank."
            )
        return FaintSelection(rank, None, DEFAULT_GAIA_RP_MAX)

    if len(values) not in {2, 3}:
        raise argparse.ArgumentTypeError(
            "Use --faint, --faint N, --faint MAG_MIN MAG_MAX, or "
            "--faint N MAG_MIN MAG_MAX."
        )

    if len(values) == 2:
        rank = 1
        mag_tokens = values
    else:
        try:
            rank_value = float(values[0])
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                "The first value in '--faint N MAG_MIN MAG_MAX' must be a "
                "positive integer source rank."
            ) from exc
        rank = int(rank_value)
        if not np.isfinite(rank_value) or rank_value != rank or rank < 1:
            raise argparse.ArgumentTypeError(
                "The first value in '--faint N MAG_MIN MAG_MAX' must be a "
                "positive integer source rank."
            )
        mag_tokens = values[1:]

    try:
        mag_min, mag_max = (float(value) for value in mag_tokens)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Faint-source magnitude limits must be numbers."
        ) from exc

    if not np.isfinite(mag_min) or not np.isfinite(mag_max) or mag_min >= mag_max:
        raise argparse.ArgumentTypeError(
            "Require finite faint-source limits with MAG_MIN < MAG_MAX."
        )

    return FaintSelection(rank, mag_min, mag_max)


def parse_epoch(value: float) -> float:
    """Validate a Julian-year epoch used for Gaia proper-motion propagation."""
    epoch = float(value)
    if not np.isfinite(epoch) or epoch <= 0.0:
        raise argparse.ArgumentTypeError("--epoch must be a positive Julian year.")
    return epoch


# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------


def fov_to_pixels(fov_arcmin: float, pixel_scale: float) -> int:
    """Convert angular side length in arcmin to pixels."""
    return max(1, int(round(fov_arcmin * 60.0 / pixel_scale)))


def format_angular_offset(value_arcsec: float, *, signed: bool = False) -> str:
    """Format an angular offset in arcsec, switching to arcmin above 2 arcmin."""
    use_arcmin = abs(value_arcsec) > 120.0
    value = value_arcsec / 60.0 if use_arcmin else value_arcsec
    unit = "arcmin" if use_arcmin else "arcsec"
    prefix = "+" if signed and value >= 0 else ""
    return f"{prefix}{value:.2f} {unit}"


def relative_geometry_values(
    origin: SkyCoord,
    destination: SkyCoord,
) -> tuple[float, float, float, float | None, float | None]:
    """Return geometry for a move from origin to destination."""
    separation_arcsec = origin.separation(destination).arcsec
    delta_ra, delta_dec = origin.spherical_offsets_to(destination)

    if separation_arcsec < 1.0e-9:
        return (
            separation_arcsec,
            delta_ra.arcsec,
            delta_dec.arcsec,
            None,
            None,
        )

    pa_deg = float(origin.position_angle(destination).wrap_at(360.0 * u.deg).deg)
    opposite_pa_deg = (pa_deg + 180.0) % 360.0
    return (
        separation_arcsec,
        delta_ra.arcsec,
        delta_dec.arcsec,
        pa_deg,
        opposite_pa_deg,
    )


def source_to_target_offset_text(
    source: SkyCoord, target: SkyCoord, label: str
) -> str:
    """Format the chart's on-sky move from a source to T."""
    _, delta_ra_arcsec, delta_dec_arcsec, _, _ = relative_geometry_values(
        source, target
    )
    return (
        f"{label} -> T offsets\n"
        f"DeltaRA={format_angular_offset(delta_ra_arcsec, signed=True)}\n"
        f"DeltaDec={format_angular_offset(delta_dec_arcsec, signed=True)}"
    )


def report_relative_geometry(coords_deg: list[tuple[float, float]]) -> None:
    """Report geometry for moves from secondary sources to the main target."""
    if len(coords_deg) <= 1:
        return

    target = SkyCoord(
        ra=coords_deg[0][0] * u.deg,
        dec=coords_deg[0][1] * u.deg,
        frame="icrs",
    )

    print()
    print("Relative geometry to T:")
    print("  (Offsets are from each source to T; +RA=east, +Dec=north.)")
    print()

    for idx, (ra, dec) in enumerate(coords_deg[1:], start=2):
        source = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
        (
            separation_arcsec,
            delta_ra_arcsec,
            delta_dec_arcsec,
            pa_deg,
            opposite_pa_deg,
        ) = relative_geometry_values(source, target)

        if pa_deg is None:
            pa_text = "undefined (coincident positions)"
        else:
            pa_text = (
                f"{pa_deg:.2f} deg; "
                f"opposite PA={opposite_pa_deg:.2f} deg"
            )

        print(f"s{idx}:")
        print(f"    separation={format_angular_offset(separation_arcsec)};")
        print(
            f"    DeltaRA={format_angular_offset(delta_ra_arcsec, signed=True)}; "
            f"DeltaDec={format_angular_offset(delta_dec_arcsec, signed=True)};"
        )
        print(f"    PA(s{idx}->T)={pa_text}")
        print()


def _optional_float(value) -> float | None:
    """Return a finite float from an Astropy/VOTable cell, or None."""
    if np.ma.is_masked(value):
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if np.isfinite(result) else None


def _format_gaia_mag(value: float | None) -> str:
    return "--" if value is None else f"{value:.2f}"


def format_epoch_value(epoch: float) -> str:
    """Format a Julian year compactly while retaining useful precision."""
    if abs(epoch - round(epoch)) < 1.0e-9:
        return f"{epoch:.1f}"
    return f"{epoch:.3f}".rstrip("0").rstrip(".")


def propagate_gaia_position(
    ra_deg: float,
    dec_deg: float,
    pmra_masyr: float | None,
    pmdec_masyr: float | None,
    ref_epoch: float | None,
    epoch: float,
) -> tuple[SkyCoord, bool]:
    """Propagate a Gaia ICRS position using the catalog proper motion.

    The propagation is performed as a two-dimensional angular offset on the
    celestial sphere.  This avoids requiring an assumed distance or radial
    velocity, which are unnecessary for acquisition-star positioning and can
    otherwise trigger ERFA ``pmsafe`` warnings for Gaia sources without full
    6-D astrometry.
    """
    base = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    if pmra_masyr is None or pmdec_masyr is None:
        return base, False

    reference_epoch = 2016.0 if ref_epoch is None else ref_epoch
    dt_years = epoch - reference_epoch
    east_mas = pmra_masyr * dt_years
    north_mas = pmdec_masyr * dt_years

    offset_mas = float(np.hypot(east_mas, north_mas))
    if not np.isfinite(offset_mas):
        return base, False
    if offset_mas == 0.0:
        return base, True

    # Position angle is measured North through East, matching Astropy's
    # directional_offset_by convention. Gaia pmra is mu_alpha* = dRA/dt*cos(dec).
    pa_rad = float(np.arctan2(east_mas, north_mas))
    try:
        propagated = base.directional_offset_by(pa_rad * u.rad, offset_mas * u.mas)
    except Exception:
        return base, False

    return propagated.icrs, True


def query_gaia_candidates(
    target: SkyCoord,
    selection: FaintSelection,
    epoch: float,
) -> list[GaiaCandidate]:
    """Query Gaia DR3 and return magnitude-selected stars sorted by separation."""
    radius_deg = GAIA_SEARCH_RADIUS_ARCSEC / 3600.0
    mag_clause = f"phot_rp_mean_mag < {selection.mag_max:.8f}"
    if selection.mag_min is not None:
        mag_clause += f" AND phot_rp_mean_mag > {selection.mag_min:.8f}"

    query = f"""
SELECT source_id, ra, dec, pmra, pmdec, ref_epoch,
       phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag
FROM {GAIA_DR3_TABLE}
WHERE 1 = CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', {target.ra.deg:.10f}, {target.dec.deg:.10f}, {radius_deg:.10f})
)
AND {mag_clause}
"""

    print(
        "Querying Gaia DR3 around T: "
        f"radius={GAIA_SEARCH_RADIUS_ARCSEC:g} arcsec, G_RP selection="
        + (
            f"{selection.mag_min:g}<G_RP<{selection.mag_max:g}"
            if selection.mag_min is not None
            else f"G_RP<{selection.mag_max:g}"
        )
    )
    try:
        response = requests.post(
            GAIA_TAP_SYNC_URL,
            data={
                "REQUEST": "doQuery",
                "LANG": "ADQL",
                # CSV keeps the TAP column names explicit and avoids relying on
                # VOTable FIELD IDs, which can vary between TAP responses.
                "FORMAT": "csv",
                "MAXREC": "2000",
                "QUERY": query,
            },
            timeout=DEFAULT_TIMEOUT,
        )
        response.raise_for_status()
    except requests.Timeout as exc:
        raise RuntimeError(
            f"Gaia DR3 TAP query timed out after {DEFAULT_TIMEOUT} seconds."
        ) from exc
    except requests.RequestException as exc:
        raise RuntimeError(f"Gaia DR3 TAP query failed: {exc}") from exc

    try:
        reader = csv.DictReader(StringIO(response.text))
        if not reader.fieldnames:
            raise ValueError("missing CSV header")
        normalized_headers = {
            str(name).strip().lower() for name in reader.fieldnames if name is not None
        }
        required_headers = {"source_id", "ra", "dec", "phot_rp_mean_mag"}
        if not required_headers.issubset(normalized_headers):
            raise ValueError(
                "unexpected Gaia CSV columns: "
                + ", ".join(sorted(normalized_headers))
            )
        rows = list(reader)
    except Exception as exc:
        preview = response.text[:300].replace("\n", " ")
        raise RuntimeError(
            "Gaia DR3 TAP returned an unreadable response. "
            f"Response begins: {preview!r}"
        ) from exc

    candidates: list[GaiaCandidate] = []
    for raw_row in rows:
        # Normalize keys defensively because TAP services may vary the case of
        # returned column names.
        row = {
            str(key).strip().lower(): value
            for key, value in raw_row.items()
            if key is not None
        }

        ra = _optional_float(row.get("ra"))
        dec = _optional_float(row.get("dec"))
        rp = _optional_float(row.get("phot_rp_mean_mag"))
        if ra is None or dec is None or rp is None:
            continue

        coord, propagated = propagate_gaia_position(
            ra,
            dec,
            _optional_float(row.get("pmra")),
            _optional_float(row.get("pmdec")),
            _optional_float(row.get("ref_epoch")),
            epoch,
        )
        separation_arcsec = float(target.separation(coord).arcsec)

        # The TAP cone is defined at the Gaia catalog epoch. Re-apply the
        # requested 90-arcsec limit after propagating positions to --epoch.
        if separation_arcsec > GAIA_SEARCH_RADIUS_ARCSEC:
            continue

        source_raw = str(row.get("source_id", "")).strip()
        if not source_raw:
            continue
        # Preserve Gaia's 64-bit source identifier exactly. Some tabular
        # serializers may append a harmless '.0'; remove only that exact suffix
        # rather than converting through float, which could lose integer precision.
        if re.fullmatch(r"[+-]?\d+\.0+", source_raw):
            source_id = source_raw.split(".", 1)[0]
        else:
            source_id = source_raw

        candidates.append(
            GaiaCandidate(
                source_id=source_id,
                g_mag=_optional_float(row.get("phot_g_mean_mag")),
                bp_mag=_optional_float(row.get("phot_bp_mean_mag")),
                rp_mag=rp,
                coord=coord,
                propagated=propagated,
                separation_arcsec=separation_arcsec,
            )
        )

    candidates.sort(key=lambda item: item.separation_arcsec)
    return candidates


def format_gaia_coord(coord: SkyCoord) -> str:
    """Return a compact sexagesimal Gaia coordinate string."""
    ra_text = coord.ra.to_string(unit=u.hourangle, sep=":", precision=3, pad=True)
    dec_text = coord.dec.to_string(
        unit=u.deg, sep=":", precision=2, alwayssign=True, pad=True
    )
    return f"{ra_text} {dec_text}"


def report_gaia_candidates(
    target: SkyCoord,
    candidates: list[GaiaCandidate],
    selection: FaintSelection,
    epoch: float,
) -> None:
    """Print up to five closest Gaia candidates plus a selected source beyond five."""
    if not candidates:
        return

    selected_index = selection.rank - 1
    selected_label = f"G{selection.rank}"
    shown_indices = list(range(min(MAX_GAIA_TERMINAL_SOURCES, len(candidates))))
    if selected_index not in shown_indices:
        shown_indices.append(selected_index)

    print()
    print("Gaia DR3 acquisition candidates around T:")
    if selection.mag_min is None:
        cut_text = f"G_RP < {selection.mag_max:g}"
    else:
        cut_text = f"{selection.mag_min:g} < G_RP < {selection.mag_max:g}"
    print(
        f"  Selection: {cut_text}; radius={GAIA_SEARCH_RADIUS_ARCSEC:g} arcsec; "
        f"epoch={format_epoch_value(epoch)}."
    )
    if len(candidates) > MAX_GAIA_TERMINAL_SOURCES:
        suffix = (
            f" Showing the 5 closest plus selected {selected_label}."
            if selected_index >= MAX_GAIA_TERMINAL_SOURCES
            else " Showing the 5 closest."
        )
        print(f"  {len(candidates)} candidates found.{suffix}")
    else:
        print(f"  {len(candidates)} candidate(s) found.")
    print("  (Offsets are from each Gaia star to T; +RA=east, +Dec=north.)")
    print()

    for index in shown_indices:
        candidate = candidates[index]
        label = f"G{index + 1}"
        selected_text = " [selected]" if index == selected_index else ""
        print(f"{label}{selected_text}:")
        print(
            f"    Gaia DR3 {candidate.source_id}; "
            f"G={_format_gaia_mag(candidate.g_mag)} mag; "
            f"G_BP={_format_gaia_mag(candidate.bp_mag)} mag; "
            f"G_RP={candidate.rp_mag:.2f} mag"
        )
        pm_note = "" if candidate.propagated else " [no PM; catalog position retained]"
        print(
            f"    RA Dec (ICRS, epoch {format_epoch_value(epoch)}): "
            f"{format_gaia_coord(candidate.coord)}{pm_note}"
        )

        (
            separation_arcsec,
            delta_ra_arcsec,
            delta_dec_arcsec,
            pa_deg,
            opposite_pa_deg,
        ) = relative_geometry_values(candidate.coord, target)

        if pa_deg is None:
            pa_text = "undefined (coincident positions)"
        else:
            pa_text = f"{pa_deg:.2f} deg; opposite PA={opposite_pa_deg:.2f} deg"

        print(f"    separation={format_angular_offset(separation_arcsec)};")
        print(
            f"    DeltaRA={format_angular_offset(delta_ra_arcsec, signed=True)}; "
            f"DeltaDec={format_angular_offset(delta_dec_arcsec, signed=True)};"
        )
        print(f"    PA({label}->T)={pa_text}")
        print()


def download_file(url: str, output: Path, timeout: int = DEFAULT_TIMEOUT) -> None:
    """Download a URL to disk."""
    print(f"Requesting: {url}")
    with requests.get(url, stream=True, timeout=timeout) as response:
        response.raise_for_status()
        with output.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=DOWNLOAD_CHUNK_SIZE):
                if chunk:
                    handle.write(chunk)


def save_requested_fits(
    source: Path,
    requested: str | None,
    retrieval: RetrievalInfo,
) -> Path | None:
    """Persist a downloaded FITS only when it is suitable for user retention."""
    if requested is None:
        return None

    destination = Path(requested).expanduser()
    if retrieval.from_hips2fits:
        print(
            "WARNING: --fits was requested, but this image was created via CDS "
            "HiPS2FITS as a finding-chart / visualization product, not as the "
            "preferred image for performing survey photometry. To reduce the risk "
            "of misuse, the requested FITS file will not be created or overwritten.",
            file=sys.stderr,
        )
        return None

    shutil.copy2(source, destination)
    print(f"Downloaded FITS: {destination}")
    return destination


def format_fov_tag(fov_arcmin: float) -> str:
    return f"{fov_arcmin:g}".replace(".", "p") + "arcmin"


def format_coordinate_tag(ra: float, dec: float) -> str:
    ra_tag = f"{ra:.3f}".replace(".", "p")
    dec_tag = f"{dec:+.3f}".replace(".", "p")
    return f"{ra_tag}_{dec_tag}"


def normalize_pa(angle: float) -> float:
    pa = round(angle % 360.0, 1)
    return 0.0 if pa >= 360.0 else pa


def format_pa_value(angle: float) -> str:
    return f"{angle:.1f}".rstrip("0").rstrip(".").replace(".", "p")


def build_output_name(
    name: str,
    spec: SurveySpec,
    band: str,
    ra: float,
    dec: float,
    fov_arcmin: float,
    extension: str,
    angle: float | None,
) -> Path:
    """Construct a descriptive output filename."""
    safe_ident = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "fc"
    safe_band = str(band).replace(" ", "").lower()

    parts = [safe_ident, spec.tag]
    if safe_band:
        parts.append(safe_band)
    parts.extend((format_coordinate_tag(ra, dec), format_fov_tag(fov_arcmin)))
    if angle is not None:
        parts.append(f"PA{format_pa_value(angle)}")

    return Path("_".join(parts) + f".{extension}")


def resolve_vhs_band(band: str) -> tuple[str, str]:
    """Return normalized VHS display band and ESO archive filter name."""
    key = band.strip().lower()
    aliases = {
        "y": ("Y", "Y"),
        "j": ("J", "J"),
        "h": ("H", "H"),
        "k": ("K", "Ks"),
        "ks": ("K", "Ks"),
        "k_s": ("K", "Ks"),
    }
    if key not in aliases:
        raise ValueError("VHS band must be one of: Y J H K (Ks/K_s accepted for K)")
    return aliases[key]


def query_vhs_products(
    ra: float,
    dec: float,
    radius_deg: float,
    eso_filter: str,
    require_full_field: bool = True,
) -> list[dict[str, object]]:
    """Query ESO ObsCore for VHS images covering the requested region."""
    if require_full_field:
        spatial_constraint = (
            "CONTAINS("
            f"CIRCLE('ICRS',{ra:.10f},{dec:.10f},{radius_deg:.10f}),"
            "s_region)=1"
        )
    else:
        spatial_constraint = (
            "CONTAINS("
            f"POINT('ICRS',{ra:.10f},{dec:.10f}),"
            "s_region)=1"
        )

    query = f"""
        SELECT TOP 50
            dp_id, filter, t_exptime, publication_date,
            s_resolution, s_fov, n_obs, access_estsize, dataproduct_subtype
        FROM ivoa.ObsCore
        WHERE obs_collection='VHS'
          AND dataproduct_type='image'
          AND filter='{eso_filter}'
          AND {spatial_constraint}
        ORDER BY n_obs DESC, t_exptime DESC, s_fov DESC, publication_date DESC
    """

    print(
        "Querying ESO archive for VHS "
        f"{eso_filter}-band image coverage at RA={ra:.8f}, Dec={dec:+.8f}"
    )
    response = requests.post(
        ESO_TAP_SYNC_URL,
        data={
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "votable",
            "MAXREC": "50",
            "QUERY": query,
        },
        timeout=DEFAULT_TIMEOUT,
    )
    response.raise_for_status()

    try:
        table = parse_single_table(BytesIO(response.content)).to_table()
    except Exception as exc:
        preview = response.text[:300].replace("\n", " ")
        raise RuntimeError(
            "ESO TAP returned an unreadable response for the VHS query. "
            f"Response begins: {preview!r}"
        ) from exc

    rows: list[dict[str, object]] = []
    for row in table:
        item: dict[str, object] = {}
        for name in table.colnames:
            value = row[name]
            if np.ma.is_masked(value):
                item[name] = None
            elif isinstance(value, bytes):
                item[name] = value.decode("utf-8", errors="replace").strip()
            elif isinstance(value, str):
                item[name] = value.strip()
            else:
                item[name] = value
        rows.append(item)
    return rows


def get_vhs_eso(
    ra: float,
    dec: float,
    fov_arcmin: float,
    band: str,
    output: Path,
) -> None:
    """Retrieve a VHS FITS cutout from the ESO archive using TAP + SODA.

    ``fov_arcmin`` is the side of the enlarged input square used by the finding
    chart.  A SODA circular cutout with radius half that side is sufficient to
    cover the final rotated square because the input side has already been
    enlarged by sqrt(2).
    """
    _, eso_filter = resolve_vhs_band(band)
    radius_deg = fov_arcmin / 120.0

    rows = query_vhs_products(
        ra, dec, radius_deg, eso_filter, require_full_field=True
    )
    if not rows:
        print(
            "WARNING: no single VHS image fully contains the requested field; "
            "retrying with products that contain the target position.",
            file=sys.stderr,
        )
        rows = query_vhs_products(
            ra, dec, radius_deg, eso_filter, require_full_field=False
        )

    if not rows:
        raise RuntimeError(
            f"ESO archive returned no VHS {eso_filter}-band image at this position."
        )

    last_error: Exception | None = None
    for row in rows:
        dp_id = row.get("dp_id")
        if dp_id is None:
            continue
        dp_id = str(dp_id).strip()
        if not dp_id:
            continue

        print(f"Requesting ESO SODA VHS cutout: ID={dp_id}")
        try:
            with requests.get(
                ESO_SODA_SYNC_URL,
                params={
                    "ID": dp_id,
                    "POS": f"CIRCLE {ra:.10f} {dec:.10f} {radius_deg:.10f}",
                },
                stream=True,
                timeout=DEFAULT_TIMEOUT,
            ) as response:
                response.raise_for_status()
                content_type = response.headers.get("Content-Type", "").lower()
                if "text/" in content_type or "html" in content_type:
                    preview = response.text[:300].replace("\n", " ")
                    raise RuntimeError(
                        "ESO SODA returned text instead of a FITS cutout: "
                        f"{preview!r}"
                    )
                with output.open("wb") as handle:
                    for chunk in response.iter_content(
                        chunk_size=DOWNLOAD_CHUNK_SIZE
                    ):
                        if chunk:
                            handle.write(chunk)

            # Verify that SODA actually returned a readable celestial FITS image.
            read_fits_image(output)
            return
        except Exception as exc:
            last_error = exc
            try:
                output.unlink(missing_ok=True)
            except OSError:
                pass
            print(
                f"WARNING: VHS product {dp_id} could not provide a usable cutout; "
                "trying another overlapping product. "
                f"Reason: {exc}",
                file=sys.stderr,
            )

    raise RuntimeError(
        "ESO archive found VHS products, but none yielded a usable FITS cutout. "
        f"Last error: {last_error}"
    )


# -----------------------------------------------------------------------------
# Survey retrieval
# -----------------------------------------------------------------------------


def get_hips2fits(
    ra: float,
    dec: float,
    side: int,
    fov_arcmin: float,
    hips_id: str,
    output: Path,
) -> None:
    """Generate a mosaicked FITS cutout from a CDS HiPS survey."""
    print(
        f"Requesting CDS HiPS2FITS: hips={hips_id}, "
        f"size={fov_arcmin:g} arcmin, pixels={side}x{side}"
    )

    if hips2fits is None:
        raise RuntimeError(
            "astroquery.hips2fits is unavailable; upgrade astroquery to use the "
            "CDS HiPS2FITS backend."
        )

    # Use the same timeout policy as the other network retrieval backends.
    hips2fits.timeout = DEFAULT_TIMEOUT
    result = hips2fits.query(
        hips=hips_id,
        width=side,
        height=side,
        projection="TAN",
        ra=Longitude(ra * u.deg),
        dec=Latitude(dec * u.deg),
        fov=Angle(fov_arcmin * u.arcmin),
        coordsys="icrs",
        format="fits",
    )

    try:
        result.writeto(output, overwrite=True)
    finally:
        result.close()


def get_ps1_stsci(ra: float, dec: float, side: int, band: str, output: Path) -> None:
    """Retrieve a native Pan-STARRS1 FITS cutout from STScI."""
    band = band.lower()
    if band not in {"g", "r", "i", "z", "y"}:
        raise ValueError("PS1 band must be one of: g r i z y")

    pos = f"{ra:.8f},{dec:.8f}"
    query_url = (
        "https://ps1images.stsci.edu/cgi-bin/ps1cutouts"
        f"?filter={band}&pos={pos}&size={side}"
    )

    print(f"Requesting PS1 cutout page: {query_url}")
    response = requests.get(query_url, timeout=DEFAULT_TIMEOUT)
    response.raise_for_status()

    match = re.search(
        r'href=["\']([^"\']+)["\'][^>]*>\s*FITS(?:\s+|-)?cutout',
        response.text,
        flags=re.IGNORECASE,
    )

    if match:
        fits_url = match.group(1)
    else:
        hrefs = re.findall(r'href=["\']([^"\']+)["\']', response.text, re.I)
        candidates = [
            href
            for href in hrefs
            if "fits" in href.lower() and "cutout" in href.lower()
        ]
        if not candidates:
            raise RuntimeError(
                "Could not find a FITS cutout link in the Pan-STARRS response."
            )
        fits_url = candidates[0]

    fits_url = urljoin(
        "https://ps1images.stsci.edu/",
        fits_url.replace("&amp;", "&"),
    )
    download_file(fits_url, output)


def get_legacy(
    ra: float,
    dec: float,
    side: int,
    band: str,
    output: Path,
    layer: str,
) -> None:
    """Retrieve a DESI Legacy Survey FITS cutout."""
    band = band.lower()
    if band not in {"g", "r", "i", "z"}:
        raise ValueError("Legacy Survey band must be one of: g r i z")

    url = (
        "https://www.legacysurvey.org/viewer/fits-cutout"
        f"?ra={ra:.8f}&dec={dec:.8f}&layer={layer}"
        f"&pixscale={SURVEYS['legacy'].pixel_scale}&size={side}&bands={band}"
    )
    download_file(url, output)


def get_skyview(
    ra: float,
    dec: float,
    side: int,
    fov_arcmin: float,
    survey_name: str,
    output: Path,
) -> None:
    """Retrieve a FITS image using NASA SkyView."""
    position = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
    angular_size = fov_arcmin * u.arcmin

    print(
        f"Requesting SkyView: survey={survey_name}, "
        f"size={fov_arcmin:g} arcmin, pixels={side}x{side}"
    )

    images = SkyView.get_images(
        position=position,
        survey=[survey_name],
        width=angular_size,
        height=angular_size,
        pixels=f"{side},{side}",
        coordinates="J2000",
        projection="Tan",
    )

    if not images:
        raise RuntimeError(f"SkyView returned no image for survey '{survey_name}'.")

    images[0].writeto(output, overwrite=True)


def get_dss2(
    ra: float,
    dec: float,
    side: int,
    fov_arcmin: float,
    band: str,
    output: Path,
) -> None:
    """Retrieve a DSS2 FITS cutout via NASA SkyView."""
    aliases = {
        "r": "DSS2 Red",
        "red": "DSS2 Red",
        "b": "DSS2 Blue",
        "blue": "DSS2 Blue",
        "i": "DSS2 IR",
        "ir": "DSS2 IR",
    }
    key = band.lower()
    if key not in aliases:
        raise ValueError("DSS2 band must be one of: red/r, blue/b, ir/i")

    get_skyview(ra, dec, side, fov_arcmin, aliases[key], output)


def get_2mass_skyview(
    ra: float,
    dec: float,
    side: int,
    fov_arcmin: float,
    band: str,
    output: Path,
) -> None:
    """Retrieve a 2MASS FITS cutout via NASA SkyView."""
    band = band.upper()
    if band not in {"J", "H", "K"}:
        raise ValueError("2MASS band must be one of: J H K")

    get_skyview(ra, dec, side, fov_arcmin, f"2MASS-{band}", output)


def get_skymapper(
    ra: float,
    dec: float,
    fov_arcmin: float,
    band: str,
    output: Path,
) -> None:
    """Retrieve a SkyMapper DR4 FITS cutout through its SIAP service."""
    band = band.lower()
    if band not in {"u", "v", "g", "r", "i", "z"}:
        raise ValueError("SkyMapper band must be one of: u v g r i z")

    if fov_arcmin >= SKYMAPPER_MAX_CUTOUT_ARCMIN:
        raise ValueError("SkyMapper DR4 cutouts must be smaller than 10 arcmin.")

    size_deg = fov_arcmin / 60.0

    def query(intersect: str) -> list[dict[str, str]]:
        params = {
            "POS": f"{ra:.8f},{dec:.8f}",
            "SIZE": f"{size_deg:.8f}",
            "BAND": band,
            "FORMAT": "image/fits",
            "INTERSECT": intersect,
            "VERB": "3",
            "RESPONSEFORMAT": "CSV",
        }
        print(
            f"Requesting SkyMapper DR4: band={band}, "
            f"size={fov_arcmin:g} arcmin, intersect={intersect}"
        )
        response = requests.get(
            SKYMAPPER_SIAP_URL,
            params=params,
            timeout=DEFAULT_TIMEOUT,
        )
        response.raise_for_status()
        return list(csv.DictReader(StringIO(response.text)))

    rows = query("COVERS")
    if not rows:
        print(
            "SkyMapper: no single image fully covers the requested field; "
            "retrying with INTERSECT=CENTER."
        )
        rows = query("CENTER")

    if not rows:
        raise RuntimeError(
            f"SkyMapper DR4 returned no {band}-band image for this position."
        )

    row = {
        str(key).strip().lower(): str(value).strip()
        for key, value in rows[0].items()
        if key is not None and value is not None
    }

    fits_url = next(
        (
            row[key]
            for key in ("get_fits", "get_image", "access_url", "access_reference")
            if row.get(key)
        ),
        None,
    )
    if fits_url is None:
        raise RuntimeError(
            "SkyMapper returned metadata but no FITS download URL."
        )

    fits_url = urljoin(SKYMAPPER_BASE_URL, fits_url.replace("&amp;", "&"))
    download_file(fits_url, output)


def retrieve_image(
    survey: str,
    ra: float,
    dec: float,
    side: int,
    fov_arcmin: float,
    final_fov_arcmin: float,
    pixel_scale_arcsec: float,
    display_pa: float,
    band: str,
    output: Path,
    legacy_layer: str,
) -> RetrievalInfo:
    """Retrieve the survey image and report which backend produced it.

    Pan-STARRS1 and 2MASS prefer their previous backends and use HiPS2FITS
    only as a fallback. SDSS, GALEX, DES, VIKING, and AllWISE are retrieved
    directly through CDS HiPS2FITS. VHS is retrieved from the ESO archive via
    TAP discovery plus a SODA FITS cutout.
    """
    if survey == "ps1":
        band = band.lower()
        if band not in HIPS_PS1:
            raise ValueError("PS1 band must be one of: g r i z y")

        native_error: Exception | None = None
        try:
            get_ps1_stsci(ra, dec, side, band, output)
            coverage = assess_fits_coverage(
                output, ra, dec, final_fov_arcmin, pixel_scale_arcsec, display_pa
            )
            report_native_coverage("STScI Pan-STARRS final-field", coverage)
            if coverage.acceptable:
                return RetrievalInfo("STScI Pan-STARRS cutout")

            print(
                "WARNING: the STScI Pan-STARRS cutout has substantial missing "
                "edge/central coverage; retrying with CDS HiPS2FITS.",
                file=sys.stderr,
            )
        except Exception as exc:
            native_error = exc
            print(
                "WARNING: STScI Pan-STARRS retrieval failed; retrying with "
                f"CDS HiPS2FITS. Reason: {exc}",
                file=sys.stderr,
            )

        hips_candidate = output.with_name(output.stem + "_hips.fits")
        try:
            get_hips2fits(
                ra, dec, side, fov_arcmin, HIPS_PS1[band], hips_candidate
            )
            hips_candidate.replace(output)
            return RetrievalInfo("CDS HiPS2FITS (Pan-STARRS1)", True)
        except Exception as hips_exc:
            try:
                hips_candidate.unlink(missing_ok=True)
            except OSError:
                pass

            if output.exists() and native_error is None:
                print(
                    "WARNING: CDS HiPS2FITS fallback also failed. The finding "
                    "chart will use the incomplete STScI Pan-STARRS image. "
                    f"Reason: {hips_exc}",
                    file=sys.stderr,
                )
                return RetrievalInfo(
                    "STScI Pan-STARRS cutout (incomplete; HiPS fallback failed)"
                )

            raise RuntimeError(
                "Pan-STARRS retrieval failed with both STScI and CDS "
                f"HiPS2FITS. STScI error: {native_error}; "
                f"HiPS2FITS error: {hips_exc}"
            ) from hips_exc

    if survey == "2mass":
        band = band.upper()
        if band not in HIPS_2MASS:
            raise ValueError("2MASS band must be one of: J H K")

        native_error: Exception | None = None
        try:
            get_2mass_skyview(ra, dec, side, fov_arcmin, band, output)
            coverage = assess_fits_coverage(
                output, ra, dec, final_fov_arcmin, pixel_scale_arcsec, display_pa
            )
            report_native_coverage("NASA SkyView 2MASS final-field", coverage)
            if coverage.acceptable:
                return RetrievalInfo("NASA SkyView (2MASS)")

            print(
                "WARNING: the SkyView 2MASS cutout has substantial missing "
                "edge/central coverage; retrying with CDS HiPS2FITS.",
                file=sys.stderr,
            )
        except Exception as exc:
            native_error = exc
            print(
                "WARNING: NASA SkyView 2MASS retrieval failed; retrying with "
                f"CDS HiPS2FITS. Reason: {exc}",
                file=sys.stderr,
            )

        hips_candidate = output.with_name(output.stem + "_hips.fits")
        try:
            get_hips2fits(
                ra, dec, side, fov_arcmin, HIPS_2MASS[band], hips_candidate
            )
            hips_candidate.replace(output)
            return RetrievalInfo("CDS HiPS2FITS (2MASS)", True)
        except Exception as hips_exc:
            try:
                hips_candidate.unlink(missing_ok=True)
            except OSError:
                pass

            if output.exists() and native_error is None:
                print(
                    "WARNING: CDS HiPS2FITS fallback also failed. The finding "
                    "chart will use the incomplete SkyView 2MASS image. "
                    f"Reason: {hips_exc}",
                    file=sys.stderr,
                )
                return RetrievalInfo(
                    "NASA SkyView (2MASS; incomplete; HiPS fallback failed)"
                )

            raise RuntimeError(
                "2MASS retrieval failed with both NASA SkyView and CDS "
                f"HiPS2FITS. SkyView error: {native_error}; "
                f"HiPS2FITS error: {hips_exc}"
            ) from hips_exc

    if survey == "vhs":
        get_vhs_eso(ra, dec, fov_arcmin, band, output)
        return RetrievalInfo("ESO Archive TAP + SODA (VHS)")

    if survey in HIPS_PRIMARY:
        _, hips_id = resolve_primary_hips_band(survey, band)
        get_hips2fits(ra, dec, side, fov_arcmin, hips_id, output)
        return RetrievalInfo(f"CDS HiPS2FITS ({SURVEYS[survey].label})", True)

    if survey == "legacy":
        get_legacy(ra, dec, side, band, output, legacy_layer)
        return RetrievalInfo("DESI Legacy Survey cutout service")
    if survey == "skymapper":
        get_skymapper(ra, dec, fov_arcmin, band, output)
        return RetrievalInfo("SkyMapper DR4 SIAP service")
    if survey == "dss2":
        get_dss2(ra, dec, side, fov_arcmin, band, output)
        return RetrievalInfo("NASA SkyView (DSS2)")

    raise ValueError(f"Unsupported survey: {survey}")


# -----------------------------------------------------------------------------
# FITS/WCS handling and image scaling
# -----------------------------------------------------------------------------


def normalize_deprecated_pc_keywords(header: fits.Header) -> fits.Header:
    """Convert legacy PC001001-style matrix keys before Astropy parses WCS."""
    normalized = header.copy()
    for keyword in list(normalized):
        match = re.fullmatch(r"PC(\d{3})(\d{3})", keyword)
        if match is None:
            continue
        standard_keyword = f"PC{int(match.group(1))}_{int(match.group(2))}"
        if standard_keyword in normalized:
            # Keep an explicitly supplied standard-form value if both exist.
            del normalized[keyword]
        else:
            normalized.rename_keyword(keyword, standard_keyword)
    return normalized


def read_fits_image(path: Path) -> tuple[np.ndarray, WCS]:
    """Read the first numeric 2-D image with a celestial WCS."""
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is None:
                continue

            try:
                data = np.asarray(hdu.data, dtype=float)
            except (TypeError, ValueError):
                continue

            while data.ndim > 2:
                data = data[0]
            data = np.squeeze(data)
            if data.ndim != 2:
                continue

            full_wcs = WCS(normalize_deprecated_pc_keywords(hdu.header))
            if not full_wcs.has_celestial:
                continue

            return data, full_wcs.celestial

    raise RuntimeError("The FITS file contains no usable 2-D image/WCS.")


def image_pixel_scale_arcsec(wcs: WCS) -> float:
    """Return the geometric-mean celestial pixel scale of a FITS image."""
    scales = proj_plane_pixel_scales(wcs.celestial)
    axis_units = wcs.celestial.wcs.cunit
    arcsec_scales = [
        (float(scale) * u.Unit(unit)).to_value(u.arcsec)
        for scale, unit in zip(scales, axis_units)
    ]
    pixel_scale = float(np.sqrt(np.prod(arcsec_scales)))
    if not np.isfinite(pixel_scale) or pixel_scale <= 0:
        raise ValueError("FITS image has no valid celestial pixel scale.")
    return pixel_scale


def assess_array_coverage(data: np.ndarray) -> CoverageInfo:
    """Assess whether non-finite pixels form a substantial missing edge region.

    The fallback criterion is intentionally conservative: isolated NaNs are
    ignored, while non-finite regions connected to an image border are counted.
    A badly missing central region also triggers the fallback independently.
    """
    if data.ndim != 2 or data.size == 0:
        return CoverageInfo(0.0, 1.0, 1.0, False)

    finite = np.isfinite(data)
    finite_fraction = float(np.mean(finite))
    invalid = ~finite

    if not np.any(invalid):
        return CoverageInfo(1.0, 0.0, 0.0, True)

    # Grow only through invalid pixels starting from invalid border pixels.
    # This distinguishes the broad edge stripes/gaps that motivate HiPS from
    # isolated bad pixels inside an otherwise useful cutout.
    seed = np.zeros_like(invalid, dtype=bool)
    seed[0, :] = invalid[0, :]
    seed[-1, :] = invalid[-1, :]
    seed[:, 0] |= invalid[:, 0]
    seed[:, -1] |= invalid[:, -1]
    edge_connected = binary_propagation(seed, mask=invalid)
    edge_missing_fraction = float(np.mean(edge_connected))

    ny, nx = data.shape
    half_y = max(1, int(round(0.5 * CENTER_CHECK_FRACTION * ny)))
    half_x = max(1, int(round(0.5 * CENTER_CHECK_FRACTION * nx)))
    cy, cx = ny // 2, nx // 2
    y0, y1 = max(0, cy - half_y), min(ny, cy + half_y + 1)
    x0, x1 = max(0, cx - half_x), min(nx, cx + half_x + 1)
    center_missing_fraction = float(np.mean(invalid[y0:y1, x0:x1]))

    acceptable = (
        finite_fraction > 0.0
        and edge_missing_fraction <= MAX_EDGE_MISSING_FRACTION
        and center_missing_fraction <= MAX_CENTER_MISSING_FRACTION
    )
    return CoverageInfo(
        finite_fraction,
        edge_missing_fraction,
        center_missing_fraction,
        acceptable,
    )


def assess_fits_coverage(
    path: Path,
    ra: float,
    dec: float,
    final_fov_arcmin: float,
    pixel_scale_arcsec: float,
    display_pa: float,
) -> CoverageInfo:
    """Assess coverage specifically over the field that will appear in the chart.

    This avoids switching to HiPS merely because the deliberately enlarged input
    cutout contains NaNs outside the final displayed field.
    """
    data, input_wcs = read_fits_image(path)
    finite_input = np.isfinite(data).astype(np.uint8)

    side = fov_to_pixels(final_fov_arcmin, pixel_scale_arcsec)
    output_wcs = make_pa_wcs(
        ra, dec, side, pixel_scale_arcsec, display_pa
    )
    covered = np.zeros((side, side), dtype=bool)

    for y0 in range(0, side, REPROJECT_BLOCK_SIZE):
        y1 = min(side, y0 + REPROJECT_BLOCK_SIZE)
        yy, xx = np.mgrid[y0:y1, 0:side]
        world_ra, world_dec = output_wcs.pixel_to_world_values(xx, yy)
        xin, yin = input_wcs.world_to_pixel_values(world_ra, world_dec)
        sampled = map_coordinates(
            finite_input,
            [yin, xin],
            order=0,
            mode="constant",
            cval=0,
            prefilter=False,
        )
        covered[y0:y1, :] = sampled > 0

    # Reuse the same edge/center logic on the actual field that will be shown.
    coverage_image = np.where(covered, 1.0, np.nan)
    return assess_array_coverage(coverage_image)


def report_native_coverage(label: str, coverage: CoverageInfo) -> None:
    """Print a compact diagnostic for a preferred native/fallback backend."""
    print(
        f"{label} coverage: {100.0 * coverage.finite_fraction:.2f}% finite; "
        f"edge-connected missing={100.0 * coverage.edge_missing_fraction:.2f}%; "
        f"center missing={100.0 * coverage.center_missing_fraction:.2f}%"
    )


def make_pa_wcs(
    ra: float,
    dec: float,
    side: int,
    pixel_scale_arcsec: float,
    angle: float,
) -> WCS:
    """Create a TAN WCS with the requested sky PA pointing upward."""
    theta = np.deg2rad(angle)
    scale_deg = pixel_scale_arcsec / 3600.0

    out_wcs = WCS(naxis=2)
    out_wcs.wcs.crpix = [(side + 1.0) / 2.0, (side + 1.0) / 2.0]
    out_wcs.wcs.crval = [ra, dec]
    out_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    out_wcs.wcs.cunit = ["deg", "deg"]
    out_wcs.wcs.cd = scale_deg * np.array(
        [
            [-np.cos(theta), np.sin(theta)],
            [np.sin(theta), np.cos(theta)],
        ]
    )
    return out_wcs


def reproject_to_pa(
    data: np.ndarray,
    input_wcs: WCS,
    ra: float,
    dec: float,
    side: int,
    pixel_scale_arcsec: float,
    angle: float,
) -> tuple[np.ndarray, WCS]:
    """Resample an image onto a square WCS rotated to the requested PA."""
    output_wcs = make_pa_wcs(
        ra,
        dec,
        side,
        pixel_scale_arcsec,
        angle,
    )
    output = np.full((side, side), np.nan, dtype=float)

    for y0 in range(0, side, REPROJECT_BLOCK_SIZE):
        y1 = min(side, y0 + REPROJECT_BLOCK_SIZE)
        yy, xx = np.mgrid[y0:y1, 0:side]

        world_ra, world_dec = output_wcs.pixel_to_world_values(xx, yy)
        xin, yin = input_wcs.world_to_pixel_values(world_ra, world_dec)

        output[y0:y1, :] = map_coordinates(
            data,
            [yin, xin],
            order=1,
            mode="constant",
            cval=np.nan,
            prefilter=False,
        )

    return output, output_wcs


def get_stretch(name: str):
    stretches = {
        "asinh": AsinhStretch(a=0.08),
        "linear": LinearStretch(),
        "sqrt": SqrtStretch(),
        "log": LogStretch(a=1000.0),
    }
    return stretches[name]


def compute_intensity_limits(
    finite: np.ndarray,
    scale: str,
    low_percentile: float,
    high_percentile: float,
    user_vmin: float | None,
    user_vmax: float | None,
) -> tuple[float, float]:
    """Determine robust display limits."""
    if finite.size == 0:
        raise RuntimeError("The image has no finite pixels in the requested field.")

    stride = max(1, finite.size // MAX_SCALING_SAMPLE)
    sample = finite[::stride]
    p_low, p_high = np.nanpercentile(sample, [low_percentile, high_percentile])

    if scale == "auto":
        _, median, std = sigma_clipped_stats(sample, sigma=3.0, maxiters=5)
        if not np.isfinite(std) or std <= 0:
            vmin, vmax = p_low, p_high
        else:
            vmin = max(median - 2.0 * std, p_low)
            vmax = min(median + 20.0 * std, p_high)
            if vmin >= vmax:
                vmin, vmax = p_low, p_high

    elif scale == "zscale":
        try:
            vmin, vmax = ZScaleInterval(contrast=0.25).get_limits(sample)
        except Exception:
            vmin, vmax = p_low, p_high

    elif scale == "percentile":
        vmin, vmax = p_low, p_high

    else:
        raise ValueError(f"Unsupported scale mode: {scale}")

    if user_vmin is not None:
        vmin = user_vmin
    if user_vmax is not None:
        vmax = user_vmax

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        raise RuntimeError(
            f"Invalid display limits: vmin={vmin!r}, vmax={vmax!r}. "
            "Try different --low/--high or --vmin/--vmax values."
        )

    return float(vmin), float(vmax)


# -----------------------------------------------------------------------------
# Markers and annotations
# -----------------------------------------------------------------------------


def source_pixel_position(
    wcs: WCS,
    coord: SkyCoord,
    nx: int,
    ny: int,
) -> tuple[float, float] | None:
    x, y = wcs.world_to_pixel(coord)
    if not (np.isfinite(x) and np.isfinite(y)):
        return None
    if not (0.0 <= x < nx and 0.0 <= y < ny):
        return None
    return float(x), float(y)


def marker_radius(nx: int, ny: int, size_factor: float) -> float:
    """Characteristic radius of compact non-slit markers in pixels."""
    return MARKER_SCALE * size_factor * max(min(nx, ny) * 0.04, 12.0)


def add_cross(
    ax,
    x: float,
    y: float,
    nx: int,
    ny: int,
    color: str,
    size_factor: float,
) -> None:
    size = min(nx, ny)
    gap = MARKER_SCALE * size_factor * max(size * 0.018, 4.0)
    arm = MARKER_SCALE * size_factor * max(size * 0.075, 18.0)
    style = dict(color=color, linewidth=1.7, solid_capstyle="butt", zorder=10)

    ax.plot([x - arm, x - gap], [y, y], **style)
    ax.plot([x + gap, x + arm], [y, y], **style)
    ax.plot([x, x], [y - arm, y - gap], **style)
    ax.plot([x, x], [y + gap, y + arm], **style)


def add_circle(
    ax,
    x: float,
    y: float,
    nx: int,
    ny: int,
    color: str,
    size_factor: float,
) -> None:
    ax.add_patch(
        Circle(
            (x, y),
            radius=marker_radius(nx, ny, size_factor),
            fill=False,
            edgecolor=color,
            linewidth=1.8,
            zorder=10,
        )
    )


def add_diamond(
    ax,
    x: float,
    y: float,
    nx: int,
    ny: int,
    color: str,
    size_factor: float,
) -> None:
    radius = marker_radius(nx, ny, size_factor)
    ax.add_patch(
        Polygon(
            [
                (x, y + radius),
                (x + radius, y),
                (x, y - radius),
                (x - radius, y),
            ],
            closed=True,
            fill=False,
            edgecolor=color,
            linewidth=1.7,
            zorder=10,
        )
    )


def add_x_marker(
    ax,
    x: float,
    y: float,
    nx: int,
    ny: int,
    color: str,
    size_factor: float,
) -> None:
    radius = marker_radius(nx, ny, size_factor)
    style = dict(color=color, linewidth=1.7, solid_capstyle="butt", zorder=10)
    ax.plot([x - radius, x + radius], [y - radius, y + radius], **style)
    ax.plot([x - radius, x + radius], [y + radius, y - radius], **style)


def add_slit(
    ax,
    x: float,
    y: float,
    pixel_scale_arcsec: float,
    slit_width_arcsec: float,
    slit_length_arcsec: float,
    color: str,
    size_factor: float,
) -> None:
    """Draw a vertical slit and circle its exact center."""
    slit_width_pix = slit_width_arcsec / pixel_scale_arcsec
    slit_length_pix = slit_length_arcsec / pixel_scale_arcsec

    ax.add_patch(
        Rectangle(
            (x - slit_width_pix / 2.0, y - slit_length_pix / 2.0),
            slit_width_pix,
            slit_length_pix,
            fill=False,
            edgecolor=color,
            linewidth=0.9,
            clip_on=True,
            zorder=10,
        )
    )

    # --marker-size changes only the center marker, not the physical slit
    # length or width.
    center_radius_pix = size_factor * max(2.5 / pixel_scale_arcsec, 4.0)
    ax.add_patch(
        Circle(
            (x, y),
            radius=center_radius_pix,
            fill=False,
            edgecolor=color,
            linewidth=1.2,
            zorder=11,
        )
    )


def add_sources(
    ax,
    wcs: WCS,
    coords: list[SkyCoord],
    markers: list[str],
    colors: list[str],
    marker_sizes: list[float],
    pixel_scale_arcsec: float,
    slit_width_arcsec: float,
    slit_length_arcsec: float,
    nx: int,
    ny: int,
) -> None:
    """Draw and label all sources."""
    base_label_offset = max(min(nx, ny) * 0.02, 8.0)

    for idx, (coord, marker, color, size_factor) in enumerate(
        zip(coords, markers, colors, marker_sizes),
        start=1,
    ):
        position = source_pixel_position(wcs, coord, nx, ny)
        if position is None:
            continue

        x, y = position
        label = "T" if idx == 1 else f"s{idx}"

        if marker == "cross":
            add_cross(ax, x, y, nx, ny, color, size_factor)
        elif marker == "circle":
            add_circle(ax, x, y, nx, ny, color, size_factor)
        elif marker == "diamond":
            add_diamond(ax, x, y, nx, ny, color, size_factor)
        elif marker == "x":
            add_x_marker(ax, x, y, nx, ny, color, size_factor)
        elif marker == "slit":
            add_slit(
                ax,
                x,
                y,
                pixel_scale_arcsec,
                slit_width_arcsec,
                slit_length_arcsec,
                color,
                size_factor,
            )

        label_offset = base_label_offset * max(0.75, min(size_factor, 2.0))
        ax.text(
            x + label_offset,
            y + label_offset,
            label,
            fontsize=10,
            fontweight="bold",
            color=color,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=1.0),
            clip_on=True,
            zorder=12,
        )


def add_acquisition_star(
    ax,
    wcs: WCS,
    coord: SkyCoord,
    label: str,
    nx: int,
    ny: int,
) -> None:
    """Draw the selected Gaia acquisition star as a labelled x marker."""
    position = source_pixel_position(wcs, coord, nx, ny)
    if position is None:
        return

    x, y = position
    add_x_marker(ax, x, y, nx, ny, ACQUISITION_COLOR, DEFAULT_MARKER_SIZE)
    label_offset = max(min(nx, ny) * 0.02, 8.0)
    ax.text(
        x + label_offset,
        y + label_offset,
        label,
        fontsize=10,
        fontweight="bold",
        color=ACQUISITION_COLOR,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=1.0),
        clip_on=True,
        zorder=12,
    )


def add_compass(
    ax,
    wcs: WCS,
    coord: SkyCoord,
    nx: int,
    ny: int,
    fov_arcmin: float,
) -> None:
    """Add North and East arrows using the WCS."""
    separation = max(fov_arcmin * 0.10, 0.05) * u.arcmin

    x0, y0 = wcs.world_to_pixel(coord)
    north = coord.directional_offset_by(0.0 * u.deg, separation)
    east = coord.directional_offset_by(90.0 * u.deg, separation)
    xn, yn = wcs.world_to_pixel(north)
    xe, ye = wcs.world_to_pixel(east)

    if not all(np.isfinite(v) for v in (x0, y0, xn, yn, xe, ye)):
        return

    anchor = np.array([0.82 * nx, 0.80 * ny], dtype=float)
    vectors = (
        (np.array([xn - x0, yn - y0], dtype=float), "N"),
        (np.array([xe - x0, ye - y0], dtype=float), "E"),
    )

    for vector, label in vectors:
        endpoint = anchor + vector
        ax.annotate(
            "",
            xy=endpoint,
            xytext=anchor,
            arrowprops=dict(arrowstyle="-|>", color=COMPASS_COLOR, linewidth=1.4),
            zorder=12,
        )
        label_pos = endpoint + 0.08 * vector
        ax.text(
            label_pos[0],
            label_pos[1],
            label,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color=COMPASS_COLOR,
            zorder=13,
        )


# -----------------------------------------------------------------------------
# Title layout
# -----------------------------------------------------------------------------


def format_title_coord(coord: SkyCoord) -> tuple[str, str]:
    """Return compact sexagesimal RA and Dec strings for title annotations."""
    ra_text = coord.ra.to_string(
        unit=u.hourangle,
        sep=":",
        precision=3,
        pad=True,
    )
    dec_text = coord.dec.to_string(
        unit=u.deg,
        sep=":",
        precision=2,
        alwayssign=True,
        pad=True,
    )
    return ra_text, dec_text


def text_fits_axis(fig, ax, text: str, fontsize: float = 10.0) -> bool:
    """Return True when a single text line fits comfortably inside the axes."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_width = ax.get_window_extent(renderer=renderer).width

    probe = ax.text(
        0.0,
        0.0,
        text,
        fontsize=fontsize,
        alpha=0.0,
        transform=ax.transAxes,
    )
    text_width = probe.get_window_extent(renderer=renderer).width
    probe.remove()
    return text_width <= 0.98 * axes_width


def choose_heading_fontsize(fig, ax, text: str) -> float:
    """Choose the largest standard heading size that fits on one line."""
    for fontsize in (14.0, 13.0, 12.0, 11.0):
        if text_fits_axis(fig, ax, text, fontsize=fontsize):
            return fontsize
    return 10.0


def build_info_lines(
    fig,
    ax,
    skycoords: list[SkyCoord],
    markers: list[str],
    slit_width_arcsec: float,
    slit_length_arcsec: float,
) -> list[str]:
    """Build at most two compact source-information lines.

    Each source is treated as one indivisible text block.  Slit information is
    attached directly to the source that carries the slit marker.  Source
    blocks are packed from left to right onto the first line and then the
    second line.  If additional sources remain after the second line is full,
    the second line ends with ``[...]``.
    """

    def source_token(idx: int, coord: SkyCoord, marker: str) -> str:
        ra_text, dec_text = format_title_coord(coord)

        if idx == 1:
            token = f"T: RA {ra_text}   Dec {dec_text}"
        else:
            # RA/Dec are obvious from the paired sexagesimal values, so avoid
            # repeating the labels for every secondary source and save space.
            token = f"s{idx}: {ra_text} {dec_text}"

        if marker == "slit":
            token += (
                f" [slit {slit_length_arcsec:g}\" x "
                f"{slit_width_arcsec:g}\"]"
            )

        return token

    tokens = [
        source_token(idx, coord, marker)
        for idx, (coord, marker) in enumerate(zip(skycoords, markers), start=1)
    ]

    separator = "   |   "

    def compose(parts: list[str], line_number: int, truncated: bool = False) -> str:
        pieces = list(parts)
        if truncated:
            pieces.append("[...]")
        text = separator.join(pieces)
        return ("    " + text) if line_number == 2 and text else text

    # Pack as many complete source blocks as possible onto line 1.
    line1_parts: list[str] = []
    next_index = 0
    while next_index < len(tokens):
        candidate = compose(line1_parts + [tokens[next_index]], line_number=1)
        if line1_parts and not text_fits_axis(fig, ax, candidate):
            break

        # Always keep T on line 1 even in the unlikely case that the complete
        # T block itself is wider than the available axis.
        line1_parts.append(tokens[next_index])
        next_index += 1
        if not text_fits_axis(fig, ax, compose(line1_parts, line_number=1)):
            break

    line1 = compose(line1_parts, line_number=1)
    if next_index >= len(tokens):
        return [line1]

    # Continue packing source blocks onto the indented second line.
    line2_parts: list[str] = []
    while next_index < len(tokens):
        candidate = compose(line2_parts + [tokens[next_index]], line_number=2)
        if not text_fits_axis(fig, ax, candidate):
            break
        line2_parts.append(tokens[next_index])
        next_index += 1

    if next_index < len(tokens):
        # Reserve room for the truncation indicator. Remove the last displayed
        # source(s), if necessary, until ``[...]`` also fits cleanly.
        while line2_parts and not text_fits_axis(
            fig,
            ax,
            compose(line2_parts, line_number=2, truncated=True),
        ):
            line2_parts.pop()

        line2 = compose(line2_parts, line_number=2, truncated=True)
    else:
        line2 = compose(line2_parts, line_number=2)

    return [line1, line2]


# -----------------------------------------------------------------------------
# Chart rendering
# -----------------------------------------------------------------------------


def make_finding_chart(
    fits_path: Path,
    output: Path,
    coords_deg: list[tuple[float, float]],
    center_deg: tuple[float, float],
    acquisition_coord_deg: tuple[float, float] | None,
    acquisition_label: str | None,
    markers: list[str],
    colors: list[str],
    marker_sizes: list[float],
    fov_arcmin: float,
    pixel_scale_arcsec: float,
    survey_label: str,
    band: str,
    angle: float | None,
    slit_width_arcsec: float,
    slit_length_arcsec: float,
    scale: str,
    stretch_name: str,
    low_percentile: float,
    high_percentile: float,
    user_vmin: float | None,
    user_vmax: float | None,
    invert: bool,
) -> None:
    data, wcs = read_fits_image(fits_path)
    ra_center, dec_center = center_deg

    # Use a controlled final orientation for every survey/backend.  When the
    # user does not specify --angle, PA=0 means North up and East left.
    display_pa = 0.0 if angle is None else angle
    output_side = fov_to_pixels(fov_arcmin, pixel_scale_arcsec)
    data, wcs = reproject_to_pa(
        data,
        wcs,
        ra_center,
        dec_center,
        output_side,
        pixel_scale_arcsec,
        display_pa,
    )

    ny, nx = data.shape
    finite = data[np.isfinite(data)]
    vmin, vmax = compute_intensity_limits(
        finite,
        scale,
        low_percentile,
        high_percentile,
        user_vmin,
        user_vmax,
    )
    print(f"Display limits: vmin={vmin:.6g}  vmax={vmax:.6g}")

    norm = ImageNormalize(
        vmin=vmin,
        vmax=vmax,
        stretch=get_stretch(stretch_name),
        clip=True,
    )

    skycoords = [
        SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
        for ra, dec in coords_deg
    ]
    target = skycoords[0]
    center_coord = SkyCoord(
        ra=ra_center * u.deg, dec=dec_center * u.deg, frame="icrs"
    )
    acquisition_coord = (
        SkyCoord(
            ra=acquisition_coord_deg[0] * u.deg,
            dec=acquisition_coord_deg[1] * u.deg,
            frame="icrs",
        )
        if acquisition_coord_deg is not None
        else None
    )

    fig = plt.figure(figsize=(8.2, 8.2))
    ax = fig.add_subplot(111, projection=wcs)
    ax.imshow(
        data,
        origin="lower",
        cmap="gray" if invert else "gray_r",
        norm=norm,
        interpolation="nearest",
    )

    # Keep the requested sky field as the fixed plotting area.  Long slits or
    # other annotations may extend beyond the image, but must not expand the
    # WCS axes or change the square field geometry.
    ax.set_xlim(-0.5, nx - 0.5)
    ax.set_ylim(-0.5, ny - 0.5)
    ax.set_autoscale_on(False)

    add_sources(
        ax,
        wcs,
        skycoords,
        markers,
        colors,
        marker_sizes,
        pixel_scale_arcsec,
        slit_width_arcsec,
        slit_length_arcsec,
        nx,
        ny,
    )
    if acquisition_coord is not None and acquisition_label is not None:
        add_acquisition_star(
            ax, wcs, acquisition_coord, acquisition_label, nx, ny
        )
    add_compass(ax, wcs, center_coord, nx, ny, fov_arcmin)

    # Show the move from the selected Gaia star, or s2, to T.
    offset_source = None
    offset_label = None
    if acquisition_coord is not None and acquisition_label is not None:
        offset_source = acquisition_coord
        offset_label = acquisition_label
    elif len(skycoords) >= 2:
        offset_source = skycoords[1]
        offset_label = "s2"

    if offset_source is not None and offset_label is not None:
        offset_text = source_to_target_offset_text(
            offset_source, target, offset_label
        )
        ax.text(
            0.02,
            0.02,
            offset_text,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=9,
            color="black",
            bbox=dict(facecolor="white", edgecolor="0.5", alpha=0.78, pad=3.0),
            zorder=14,
        )

    ax.coords[0].set_axislabel("RA (J2000)")
    ax.coords[1].set_axislabel("Dec (J2000)")
    ax.coords.grid(color="0.65", linestyle=":", linewidth=0.6, alpha=0.65)

    field_heading = f"{survey_label}  |  {fov_arcmin:g}x{fov_arcmin:g} arcmin"
    if band:
        field_heading = (
            f"{survey_label}  |  band {band}  |  "
            f"{fov_arcmin:g}x{fov_arcmin:g} arcmin"
        )
    if angle is not None:
        pa_text = f"{angle:.1f}".rstrip("0").rstrip(".")
        field_heading += f"  |  PA {pa_text} deg"

    heading_fontsize = choose_heading_fontsize(fig, ax, field_heading)
    fig.suptitle(
        field_heading,
        fontsize=heading_fontsize,
        y=0.972,
    )

    info_lines = build_info_lines(
        fig,
        ax,
        skycoords,
        markers,
        slit_width_arcsec,
        slit_length_arcsec,
    )
    title_artist = ax.set_title("\n".join(info_lines), fontsize=10, pad=10)
    if len(info_lines) > 1:
        title_artist.set_multialignment("left")

    fig.savefig(
        output,
        format=output.suffix.lstrip("."),
        dpi=220,
        bbox_inches="tight",
        pad_inches=0.08,
        facecolor="white",
    )
    plt.close(fig)


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create annotated astronomical finding charts from public sky surveys."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -c 17:05:35.520 -23:27:21.60
  %(prog)s -c 256.398 -23.456 -f 5 -n target1 -s skymapper
  %(prog)s -c 256.398 -23.456 -i local_image.fits -n target1
  %(prog)s -c 256.398 -23.456 -f 10 -s allwise -b w1
  %(prog)s -c 256.398 -23.456 -f 6 -s vhs -b K
  %(prog)s -c 17:05:35.520 -23:27:21.60 -f 5 -a 45 -m slit --slit-length 180 --slit-width 0.8
  %(prog)s -c 17:05:35.520 -23:27:21.60 17:05:36.100 -23:27:10.0 \
           17:05:34.800 -23:27:40.0 -f 5 -a 45 -m x slit diamond
  %(prog)s -c 17:05:35.520 -23:27:21.60 17:05:36.100 -23:27:10.0 \
           -f 5 --color red orange --marker-size 1.2 0.8
  %(prog)s -c 17:05:35.520 -23:27:21.60 -f 5 -s legacy --fits field.fits
  %(prog)s -c 17:05:35.520 -23:27:21.60 --faint
  %(prog)s -c 17:05:35.520 -23:27:21.60 -F 2 14 18 --epoch 2026.5
        """,
    )

    basic = parser.add_argument_group("Basic field options")
    basic.add_argument(
        "-c",
        "--coordinates",
        metavar="RA DEC [RA2 DEC2 ...]",
        required=True,
        help=(
            "One or more RA/Dec pairs. The first source is T; additional "
            "sources are s2, s3, ... (maximum 10)."
        ),
    )
    basic.add_argument(
        "-f",
        "--fov",
        default=DEFAULT_FOV_ARCMIN,
        metavar="ARCMIN",
        help=f"Square field side in arcminutes (default: {DEFAULT_FOV_ARCMIN:g}).",
    )
    basic.add_argument(
        "-n",
        "--name",
        default="fc",
        help="Name/prefix for the output filename (default: fc).",
    )
    basic.add_argument(
        "-i",
        "--image",
        metavar="FILENAME",
        default=None,
        help="Use a local FITS image with celestial WCS instead of downloading a survey image.",
    )

    survey_group = parser.add_argument_group("Survey options")
    survey_group.add_argument(
        "-s",
        "--survey",
        choices=SURVEY_CHOICES,
        default=None,
        help="Survey to use (default: ps1).",
    )
    survey_group.add_argument(
        "--legacy-layer",
        default=None,
        help="Legacy Survey viewer layer; used only with -s legacy (default: ls-dr10).",
    )
    survey_group.add_argument(
        "-b",
        "--band",
        default=None,
        help="Filter/band; survey-dependent default.",
    )

    geometry = parser.add_argument_group("Orientation and markers")
    geometry.add_argument(
        "-a",
        "--angle",
        type=parse_angle,
        default=None,
        help=(
            "Position angle in degrees, North through East. The chart is "
            "rotated so this PA points upward. If omitted, PA=0 is used for "
            "display (North up, East left)."
        ),
    )
    geometry.add_argument(
        "-F",
        "--faint",
        "--faint-source",
        nargs="*",
        default=None,
        metavar="VALUE",
        help=(
            "Faint-target acquisition mode using Gaia DR3 stars within 90 arcsec. "
            "Forms: --faint; --faint N; --faint MAG_MIN MAG_MAX; or "
            "--faint N MAG_MIN MAG_MAX. Magnitude cuts use G_RP. The selected "
            "Gaia star becomes the chart center and sets the automatic PA unless "
            "--angle is supplied."
        ),
    )
    geometry.add_argument(
        "--epoch",
        type=float,
        default=DEFAULT_GAIA_EPOCH,
        metavar="YEAR",
        help=(
            "Julian year for Gaia proper-motion propagation in faint-source mode "
            f"(default: {DEFAULT_GAIA_EPOCH:.1f})."
        ),
    )
    geometry.add_argument(
        "-m",
        "--marker",
        nargs="+",
        type=str.lower,
        choices=ALLOWED_MARKERS,
        default=None,
        help=(
            "One marker per source if desired: x, cross, circle, diamond, slit. "
            "Missing markers default to x."
        ),
    )
    geometry.add_argument(
        "--color",
        nargs="+",
        default=None,
        help=(
            "Marker color(s), using Matplotlib color names/specifications. One "
            "value applies to every source; multiple values apply in source order. "
            "Defaults use a color-blind-friendly palette; explicit --color values override it."
        ),
    )
    geometry.add_argument(
        "--marker-size",
        nargs="+",
        type=float,
        default=None,
        metavar="FACTOR",
        help=(
            "Marker-size factor(s) relative to the default. One value applies to "
            "every source; multiple values apply in source order. For slit markers, "
            "this changes only the central circle, not the physical slit dimensions."
        ),
    )
    geometry.add_argument(
        "--slit-length",
        type=parse_slit_length,
        default=DEFAULT_SLIT_LENGTH_ARCSEC,
        metavar="ARCSEC",
        help=(
            "Slit length in arcsec; must be positive "
            f"(default: {DEFAULT_SLIT_LENGTH_ARCSEC:g}, i.e. 4 arcmin)."
        ),
    )
    geometry.add_argument(
        "--slit-width",
        type=parse_slit_width,
        default=DEFAULT_SLIT_WIDTH_ARCSEC,
        metavar="ARCSEC",
        help=(
            f"Slit width in arcsec, from {MIN_SLIT_WIDTH_ARCSEC:g} to "
            f"{MAX_SLIT_WIDTH_ARCSEC:g} "
            f"(default: {DEFAULT_SLIT_WIDTH_ARCSEC:g})."
        ),
    )

    display = parser.add_argument_group("Image display")
    display.add_argument(
        "--scale",
        choices=("auto", "zscale", "percentile"),
        default="auto",
        help="Intensity-limit method (default: auto).",
    )
    display.add_argument(
        "--stretch",
        choices=("asinh", "linear", "sqrt", "log"),
        default="asinh",
        help="Display stretch (default: asinh).",
    )
    display.add_argument(
        "--invert",
        action="store_true",
        help="Invert the grayscale image to use a dark sky background.",
    )
    display.add_argument(
        "--low",
        type=float,
        default=1.0,
        help="Lower percentile; also controls contrast in auto mode (default: 1).",
    )
    display.add_argument(
        "--high",
        type=float,
        default=99.5,
        help=(
            "Upper percentile; also controls contrast in auto mode "
            "(default: 99.5)."
        ),
    )
    display.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Absolute lower display limit; overrides the automatic value.",
    )
    display.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Absolute upper display limit; overrides the automatic value.",
    )

    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "-e",
        "--extension",
        choices=("png", "pdf", "jpg", "eps"),
        default="png",
        help="Finding-chart output format (default: png).",
    )
    output_group.add_argument(
        "--fits",
        metavar="FILENAME",
        default=None,
        help=(
            "Save the downloaded survey FITS image to FILENAME when a native "
            "survey/fallback backend is used. HiPS2FITS visualization products "
            "are intentionally not saved. With --angle, saved FITS files are "
            "the original enlarged downloads before PA reprojection."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    return parser


def main() -> int:
    parser = make_parser()
    args = parser.parse_args(normalize_coordinate_argv(sys.argv[1:]))

    try:
        coords = parse_coordinates(args.coordinates)
        fov_arcmin = parse_fov(args.fov)
        faint_selection = parse_faint_selection(args.faint)
        epoch = parse_epoch(args.epoch)
        markers = resolve_markers(args.marker, len(coords))
        colors = resolve_colors(args.color, len(coords))
        marker_sizes = resolve_marker_sizes(args.marker_size, len(coords))
    except argparse.ArgumentTypeError as exc:
        parser.error(str(exc))

    # In faint-source mode a slit is the default marker for T, but any marker
    # explicitly supplied by the user takes precedence.
    if faint_selection is not None and args.marker is None:
        markers[0] = "slit"

    image_path = Path(args.image).expanduser() if args.image is not None else None
    if image_path is not None:
        if args.fits is not None:
            parser.error("--fits saves a downloaded image and cannot be used with --image.")
        if args.survey is not None or args.band is not None or args.legacy_layer is not None:
            parser.error("--survey, --band, and --legacy-layer cannot be used with --image.")
        if not image_path.is_file():
            parser.error(f"FITS image does not exist: {image_path}")
        try:
            _, image_wcs = read_fits_image(image_path)
            pixel_scale = image_pixel_scale_arcsec(image_wcs)
        except Exception as exc:
            parser.error(f"Cannot use FITS image {image_path}: {exc}")
        survey = "image"
        spec = SurveySpec("User FITS image", "image", "", pixel_scale)
        band = ""
    else:
        survey = canonical_survey(args.survey or "ps1")
        spec = SURVEYS[survey]
        band = args.band if args.band is not None else spec.default_band
    user_angle = normalize_pa(args.angle) if args.angle is not None else None

    if not 0.0 <= args.low < args.high <= 100.0:
        parser.error("Require 0 <= --low < --high <= 100.")
    if args.vmin is not None and args.vmax is not None and args.vmin >= args.vmax:
        parser.error("--vmin must be smaller than --vmax.")

    ra_main, dec_main = coords[0]
    target_coord = SkyCoord(
        ra=ra_main * u.deg, dec=dec_main * u.deg, frame="icrs"
    )

    gaia_candidates: list[GaiaCandidate] = []
    selected_gaia: GaiaCandidate | None = None
    acquisition_label: str | None = None
    center_coord = target_coord
    automatic_pa: float | None = None

    if faint_selection is not None:
        try:
            gaia_candidates = query_gaia_candidates(
                target_coord, faint_selection, epoch
            )
        except Exception as exc:
            print(f"ERROR: Gaia DR3 query failed: {exc}", file=sys.stderr)
            return 1

        if not gaia_candidates:
            if faint_selection.mag_min is None:
                cut_text = f"G_RP < {faint_selection.mag_max:g}"
            else:
                cut_text = (
                    f"{faint_selection.mag_min:g} < G_RP < "
                    f"{faint_selection.mag_max:g}"
                )
            print(
                "ERROR: no Gaia DR3 stars satisfying "
                f"{cut_text} were found within "
                f"{GAIA_SEARCH_RADIUS_ARCSEC:g} arcsec of T.",
                file=sys.stderr,
            )
            return 1

        if faint_selection.rank > len(gaia_candidates):
            print(
                f"ERROR: --faint {faint_selection.rank} requested, but only "
                f"{len(gaia_candidates)} Gaia DR3 star(s) satisfying the "
                "magnitude selection were found within 90 arcsec of T.",
                file=sys.stderr,
            )
            return 1

        selected_gaia = gaia_candidates[faint_selection.rank - 1]
        acquisition_label = f"G{faint_selection.rank}"
        center_coord = selected_gaia.coord

        # Use the selected acquisition star -> target direction so T lies along
        # the vertical slit axis. An explicit --angle overrides this value.
        if selected_gaia.separation_arcsec < 1.0e-9:
            automatic_pa = 0.0
            print(
                "WARNING: T and the selected Gaia acquisition star are "
                "coincident; using PA=0 deg.",
                file=sys.stderr,
            )
        else:
            automatic_pa = float(
                center_coord.position_angle(target_coord)
                .wrap_at(360.0 * u.deg)
                .deg
            )

    angle = user_angle if user_angle is not None else automatic_pa
    display_pa = 0.0 if angle is None else angle

    # Downloaded square images need extra coverage for PA reprojection.
    download_fov = fov_arcmin * np.sqrt(2.0)
    download_side = fov_to_pixels(download_fov, spec.pixel_scale)
    output_side = fov_to_pixels(fov_arcmin, spec.pixel_scale)

    if survey == "skymapper" and download_fov >= SKYMAPPER_MAX_CUTOUT_ARCMIN:
        suffix = " after orientation enlargement"
        parser.error(
            "SkyMapper DR4 restricts cutouts to <10 arcmin on a side; "
            f"this request needs {download_fov:.3g} arcmin{suffix}."
        )

    if selected_gaia is not None:
        if selected_gaia.separation_arcsec > 0.5 * fov_arcmin * 60.0:
            print(
                "WARNING: Requested FOV is too small to include T and the "
                "acquisition star.",
                file=sys.stderr,
            )
        if (
            markers[0] == "slit"
            and selected_gaia.separation_arcsec > 0.5 * args.slit_length
        ):
            print(
                "WARNING: Requested slit is too short to include both T and "
                "the acquisition star.",
                file=sys.stderr,
            )

    output = build_output_name(
        args.name,
        spec,
        band,
        ra_main,
        dec_main,
        fov_arcmin,
        args.extension,
        angle,
    )

    print(f"Main target: RA={ra_main:.8f} deg  Dec={dec_main:+.8f} deg")
    for idx, (ra, dec) in enumerate(coords[1:], start=2):
        print(f"  s{idx}: RA={ra:.8f} deg  Dec={dec:+.8f} deg")

    if selected_gaia is not None and acquisition_label is not None:
        print(
            f"Faint-source acquisition: {acquisition_label} = Gaia DR3 "
            f"{selected_gaia.source_id}; chart centered on the propagated "
            f"epoch-{format_epoch_value(epoch)} position."
        )
        print(
            f"Automatic faint-source PA ({acquisition_label}->T): "
            f"{automatic_pa:.2f} deg"
        )
        if user_angle is not None:
            print(
                f"User --angle overrides automatic PA: {user_angle:.2f} deg"
            )

    print(f"Markers: {' '.join(markers)}")
    print(f"Colors: {' '.join(colors)}")
    print("Marker sizes: " + " ".join(f"{size:g}" for size in marker_sizes))
    if image_path is not None:
        print(f"Input FITS image: {image_path}")
        print(f"Input pixel scale: {spec.pixel_scale:.4g} arcsec/pixel")
    else:
        print(f"Survey: {spec.label}; band: {band}")
    print(f"FOV: {fov_arcmin:g} x {fov_arcmin:g} arcmin")
    print(f"Output sampling: {output_side} x {output_side} pixels")
    if image_path is None:
        print(
            f"Survey download: {download_fov:.4g} arcmin, "
            f"{download_side} x {download_side} pixels"
        )
    if angle is not None:
        print(f"Position angle: {angle:g} deg (North through East; PA points up)")
    else:
        print("Orientation: PA 0 deg (North up, East left)")
    if "slit" in markers:
        print(
            f"Slit: length={args.slit_length:g} arcsec, "
            f"width={args.slit_width:g} arcsec"
        )
    print(
        f"Intensity: scale={args.scale}, stretch={args.stretch}, "
        f"low={args.low:g}, high={args.high:g}"
    )
    if args.invert:
        print("Grayscale: inverted (dark background)")

    try:
        with tempfile.TemporaryDirectory(prefix="finding_chart_") as tmpdir:
            fits_path = image_path
            retrieval = None
            if fits_path is None:
                fits_path = Path(tmpdir) / "survey_cutout.fits"
                retrieval = retrieve_image(
                    survey,
                    center_coord.ra.deg,
                    center_coord.dec.deg,
                    download_side,
                    download_fov,
                    fov_arcmin,
                    spec.pixel_scale,
                    display_pa,
                    band,
                    fits_path,
                    args.legacy_layer or "ls-dr10",
                )
                print(f"Retrieval backend: {retrieval.backend}")

            acquisition_coord_deg = (
                (selected_gaia.coord.ra.deg, selected_gaia.coord.dec.deg)
                if selected_gaia is not None
                else None
            )
            make_finding_chart(
                fits_path,
                output,
                coords,
                (center_coord.ra.deg, center_coord.dec.deg),
                acquisition_coord_deg,
                acquisition_label,
                markers,
                colors,
                marker_sizes,
                fov_arcmin,
                spec.pixel_scale,
                spec.label,
                band,
                angle,
                args.slit_width,
                args.slit_length,
                args.scale,
                args.stretch,
                args.low,
                args.high,
                args.vmin,
                args.vmax,
                args.invert,
            )

            if retrieval is not None:
                save_requested_fits(fits_path, args.fits, retrieval)

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    # Print source-to-target geometry only after retrieval/output warnings have
    # already been emitted, so it appears as a compact final observing aid.
    sys.stderr.flush()
    if faint_selection is not None:
        report_gaia_candidates(
            target_coord, gaia_candidates, faint_selection, epoch
        )
    if len(coords) > 1:
        report_relative_geometry(coords)

    print(f"Finding chart: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
