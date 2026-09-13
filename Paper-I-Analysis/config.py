"""
Single source of truth for paths and constants used by every notebook in
Paper-I-Analysis (Singh et al., "Linear Polarimetry of the Solar Corona during
the Total Solar Eclipse of 8 April 2024").

Everything here was extracted from the legacy notebooks that produced the
paper figures; see README.md for the provenance of each value.
"""

import os
import tomllib
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent            # .../solar-corona-analysis/Paper-I-Analysis


def _resolve_data_root() -> Path:
    """Where the large, non-versioned eclipse data lives (Calibrated_Lights, legacy products).

    1. $SOLAR_ECLIPSE_DATA_ROOT
    2. local_paths.toml next to this file (git-ignored; see local_paths.example.toml)
    3. the SSD path used by the authors
    """
    env = os.environ.get("SOLAR_ECLIPSE_DATA_ROOT")
    if env:
        root = Path(env).expanduser()
    elif (ROOT / "local_paths.toml").exists():
        with open(ROOT / "local_paths.toml", "rb") as f:
            root = Path(tomllib.load(f)["data_root"]).expanduser()
    else:
        root = Path("/Volumes/WD_SSD/Solar Eclipse 2024 Data")
    if not (root / "Calibrated_Lights").is_dir():
        raise FileNotFoundError(
            f"Data root {root} has no Calibrated_Lights/ — set SOLAR_ECLIPSE_DATA_ROOT or "
            f"edit {ROOT / 'local_paths.toml'} (see local_paths.example.toml)")
    return root


DATA_ROOT = _resolve_data_root()                  # .../Solar Eclipse 2024 Data
SSD_ROOT = DATA_ROOT                              # backwards-compatible alias

# Input: APP-calibrated (dark, flat, debayer) light frames, 3-plane RGB float32
CALIBRATED_LIGHTS_DIR = DATA_ROOT / "Calibrated_Lights"

DATA_DIR = ROOT / "data"                          # committed CSVs + ephemeris
PRODUCTS_DIR = ROOT / "products"                  # everything regenerated
STACKED_DIR = PRODUCTS_DIR / "stacked"
HDR_DIR = PRODUCTS_DIR / "hdr"
FIGURES_DIR = ROOT / "figures"

EPHEMERIS_FILE = DATA_DIR / "de421.bsp"

# Regenerated tables
FRAME_TABLE_CSV = PRODUCTS_DIR / "frame_table.csv"
MOON_CENTERS_CSV = PRODUCTS_DIR / "moon_centers.csv"
SUN_MOON_CENTERS_CSV = PRODUCTS_DIR / "sun_moon_centers.csv"
PAPER_METADATA_CSV = PRODUCTS_DIR / "paper_metadata_table.csv"
CELESTIAL_NORTH_FIT_CSV = PRODUCTS_DIR / "celestial_north_fit.csv"
ERROR_BUDGET_CSV = PRODUCTS_DIR / "error_budget.csv"
REFERENCE_STARS_CSV = DATA_DIR / "reference_stars.csv"    # cached catalogue positions (notebook 03)
DOCS_DIR = ROOT / "docs"

# The tables that were actually used for the paper (copied from the SSD, for
# validation and as a fallback if steps 02/03 are skipped).
LEGACY_MOON_CENTERS_CSV = DATA_DIR / "sheet_for_centers_and_radius_and_SB.csv"
LEGACY_SUN_MOON_CENTERS_CSV = DATA_DIR / "sun_moon_centers_output.csv"
LEGACY_PAPER_METADATA_CSV = DATA_DIR / "paper_metadata_documentation.csv"

# Legacy products on the SSD, used only by 99_validate_against_legacy.ipynb
LEGACY_DIR = SSD_ROOT / "FINAL_PAPER_ANALYSIS"
LEGACY_STACKED_DIR = LEGACY_DIR / "Linear_Composites_Instrument_Frame_Chaitanya"
LEGACY_HDR_EXPNORM_DIR = LEGACY_DIR / "ExpNorm_HDR_Images_Paras"
LEGACY_HDR_LDIC_DIR = LEGACY_DIR / "LDIC_HDR_Images_Paras"
LEGACY_FIGURES = {
    "expnorm": LEGACY_DIR / "expnorm_hdr_solar_frame.png",
    "ldic": LEGACY_DIR / "ldic_hdr_solar_frame.png",
    "histogram": LEGACY_DIR / "histogram_normalized_exposures.png",
}

for _d in (PRODUCTS_DIR, STACKED_DIR, HDR_DIR, FIGURES_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------
# Observation
# --------------------------------------------------------------------------
# Dallas, Texas, USA (Section 2 of the paper)
SITE_LAT_DMS = "32°47'08.34'' N"
SITE_LON_DMS = "96°46'57.36'' W"
SITE_ALT_M = 130

# The camera clock was set to Indian Standard Time; DATE-OBS in the calibrated
# FITS headers is in that zone and is converted to UTC before any ephemeris.
CAMERA_CLOCK_TZ = "Asia/Kolkata"

# Mid-totality, used for the P-angle and for the Sun's angular radius.
TOTALITY_MID_UTC = "2024-04-08T18:42:00"

# Filter-wheel positions.  Position4 was the unpolarised / diamond-ring slot
# and is excluded from the polarimetric pipeline.
POLARIZER_POSITIONS = ("1", "2", "3")
EXCLUDED_POSITIONS = ("4",)
POLARIZER_ANGLE_DEG = {"1": 0.0, "2": +60.0, "3": -60.0}   # w.r.t. sensor x axis

# Exposure bracket, as 1/t denominators, ordered short -> long as in the legacy
# HDR notebooks.
INV_EXPOSURES = (800, 400, 200, 100, 50, 25, 13, 6, 3)
EXPOSURE_TIMES_S = tuple(1.0 / e for e in INV_EXPOSURES)

# Physical radii used for the angular-size computation in step 03
SUN_RADIUS_KM = 696340.0
MOON_RADIUS_KM = 1737.4

# --------------------------------------------------------------------------
# Image geometry
# --------------------------------------------------------------------------
IMAGE_SHAPE = (5504, 8256)                        # (rows, cols) of a calibrated frame

# Celestial North in the instrument frame: 168 deg CCW from the +x axis as the
# image is displayed in camera orientation (row 0 at the top).  Found manually
# in Stellarium / Affinity Photo (Fig. 5 of the paper) and confirmed to 0.1 deg
# by the star astrometry in 03_celestial_north_from_stars.ipynb.
#
# Convention used everywhere in this pipeline, in ARRAY coordinates (x = column
# to the right, y = row downwards):
#     North unit vector  n = (cos a, -sin a)
#     East  unit vector  e = (-sin a, -cos a)
#     a sky offset of (separation, position angle PA east of north) maps to
#     pixel offset  sep * (sin PA * e + cos PA * n)
CELESTIAL_NORTH_DEG = 168.0
IMAGE_NORTH_ANGLE_FROM_Y_DEG = 78.0               # the legacy code's parameterisation (see below)

# How the Sun-Moon displacement is rotated into the image in step 04:
#   "astrometric" : the convention above (validated on the background stars) — default
#   "legacy"      : dx = sep*sin(PA-78), dy = sep*cos(PA-78) as in the legacy
#                   find_sun_center notebook.  Not consistent with the star field;
#                   the resulting 3-11 px Sun-centre error differs between polariser
#                   positions and was compensated downstream by the by-eye
#                   LEGACY_CHANNEL_SHIFTS.  Use it only to reproduce the submitted figures.
SUN_CENTER_ROTATION = "astrometric"

# Background stars used to verify the orientation (notebook 03).  Positions are
# fetched from Sesame once and cached in data/reference_stars.csv.
REFERENCE_STARS = ("zet Psc", "88 Psc")
ZET_PSC_SEPARATION_ARCSEC = 23.0                  # zeta Psc A-B, blended in our PSF
STAR_SEARCH_MIN_RADIUS_PX = 1800                  # ignore the corona when looking for stars

# Solar P-angle on the day (SunPy); 00_constants.ipynb recomputes it.
SOLAR_P_ANGLE_DEG = -26.24
SOLAR_NORTH_ANGLE_DEG = CELESTIAL_NORTH_DEG + SOLAR_P_ANGLE_DEG   # 141.76
SOLAR_FRAME_PHI_DEG = -SOLAR_NORTH_ANGLE_DEG                      # -141.76

PLATE_SCALE_ARCSEC_PER_PIX = 1.49

# Sun centre in every stacked / HDR product: the array centre (0-based, as the
# legacy LDIC notebook defined it with `len(img[0]) // 2`, `len(img) // 2`).
SUN_CENTER_XY = (IMAGE_SHAPE[1] // 2, IMAGE_SHAPE[0] // 2)         # (4128, 2752)

# Radius of the black disk drawn over the Moon in the "moon-masked" figures
MOON_MASK_RADIUS_PX = 695

# --------------------------------------------------------------------------
# Pipeline parameters
# --------------------------------------------------------------------------
# Step 03 — which Moon-centre table to feed the Sun-centre computation:
#   "legacy"     : data/sheet_for_centers_and_radius_and_SB.csv, the CHT run
#                  that produced the paper (default; reproduces the paper exactly)
#   "recomputed" : products/moon_centers.csv written by 02_moon_center_cht.ipynb
#                  (current scikit-image gives +/-1 px differences)
MOON_CENTERS_SOURCE = "legacy"

# Step 04 — frames whose CHT lunar radius lies outside mean +/- N*sigma are dropped
RADIUS_FILTER_NSIGMA = 1.0

# Step 04 — frames additionally left out of the stacks that produced the paper.
# Recovered by non-negative least-squares fits of the legacy stacks onto their
# candidate frames (see README, "The stacking step"); C. Gandhi's off-line run
# excluded these on top of the 1-sigma filter.  Set to {} to stack all 130
# frames that pass the radius filter.
EXTRA_EXCLUDED_FRAMES = {
    "DSC_3810-Position1-cal.fits": "1/50 s; first frame of totality (18:40:58 UT), 3.4x brighter than its siblings (diamond ring)",
    "DSC_3868-Position2-cal.fits": "1/3 s; first Position2 frame, lunar radius 666 px vs 669-672 px in the group",
    "DSC_3935-Position3-cal.fits": "1/100 s; lunar radius 663 px vs 675-677 px in the group (plate scale / Sun centre off)",
    "DSC_3952-Position3-cal.fits": "1/200 s; lunar radius 667 px vs 676 px in the group",
    "DSC_3869-Position2-cal.fits": "1/800 s; excluded in the legacy stack (lunar radius 680 px)",
    "DSC_3887-Position2-cal.fits": "1/800 s; excluded in the legacy stack (lunar radius 680 px)",
}

# Step 02 — circular Hough transform (two-pass: coarse on a downsampled image,
# then +/-100 px refinement at full resolution)
CHT_RESCALE_FACTOR = 10
CHT_RADII_RANGE_PX = np.arange(100, 1000, 10)
CHT_REFINE_HALFWIDTH_PX = 100

# Step 04 — synthetic luminance (ITU-R BT.709 / sRGB) and stacking
LUMINANCE_WEIGHTS = (0.2126, 0.7152, 0.0722)      # R, G, B
STACK_INTERPOLATION_ORDER = 1                     # bilinear sub-pixel shift

# Steps 06/07 — inter-polariser pixel offsets (dy, dx) applied to the stacked
# images before HDR construction.  The legacy notebooks used by-eye values "to
# align the plumes and make them white"; notebook 10 shows they compensate the
# legacy Sun-centre rotation error, so with the astrometric convention no shift
# is applied.
LEGACY_CHANNEL_SHIFTS = {"1": (3, 2), "2": (0, 0), "3": (-4, -7)}
CHANNEL_SHIFTS = (LEGACY_CHANNEL_SHIFTS if SUN_CENTER_ROTATION == "legacy"
                  else {"1": (0, 0), "2": (0, 0), "3": (0, 0)})

# Notebook 10 — annulus (in solar radii, plate scale below) used for the
# frame-to-stack residual statistics
ERROR_ANNULUS_RSUN = (1.2, 2.0)

# Step 06 — Linear Digital Image Composer (Druckmullerova 2014), Section 3.3.2
LDIC_PARAMS = dict(
    num_angular_segments=60,
    trig_poly_order_k=4,
    trig_poly_order_q=4,
    wf_low_reject_percent=2.0,
    wf_low_full_percent=15.0,
    wf_high_start_reject_percent=65.0,
    wf_high_reject_percent=80.0,
    regression_weight_threshold=0.8,
    min_valid_segments=15,
    fit_intercept=False,
)
LDIC_TARGET_MAX = 65535.0

# Step 07 — display composite + MGN (Morgan & Druckmuller 2014)
MGN_SIGMAS = (10, 20, 40, 80)
COMPOSITE_GAMMA = 0.4
COMPOSITE_LOWER_PCT = 0.5
COMPOSITE_UPPER_PCT = 99.9
COMPOSITE_BOOST_PCT = 99.8
COLORWHEEL_SATURATION = 0.45
FIGURE_DPI = 300
