"""
Shared helpers for the Paper-I-Analysis notebooks.

Every function here is a de-duplicated copy of code from the legacy notebooks
on the SSD that produced the paper figures (see README.md, "Provenance").  The
algorithms are unchanged; only names, docstrings and the config plumbing are new.

Sections
    1. Frame inventory / FITS I/O
    2. Circular Hough transform lunar-limb fit            (step 02)
    3. Skyfield Sun-centre from the Moon centre           (step 03)
    4. Alignment and stacking                             (step 04)
    5. HDR: exposure normalisation and LDIC               (steps 05, 06)
    6. Polarimetry: Malus-basis <-> Stokes, frame rotation
    7. Display composites, MGN, figure annotations        (steps 07, 08)
    8. Validation helpers                                 (step 99)
"""

from __future__ import annotations

import datetime as dt
import math
import os
import re
from pathlib import Path
from zoneinfo import ZoneInfo

import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from dateutil import parser as date_parser
from matplotlib.projections import get_projection_class
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import ndimage
from skimage.feature import canny
from skimage.transform import hough_circle, hough_circle_peaks, resize

import config

# ==========================================================================
# 1. Frame inventory / FITS I/O
# ==========================================================================

_FRAME_RE = re.compile(r"^(DSC_\d+)-Position(\d)-cal\.fits$")


def list_calibrated_frames(lights_dir: Path | None = None) -> list[str]:
    """Basenames of all ``DSC_*-Position?-cal.fits`` frames, sorted."""
    lights_dir = Path(lights_dir or config.CALIBRATED_LIGHTS_DIR)
    names = [f for f in os.listdir(lights_dir) if _FRAME_RE.match(f)]
    return sorted(names)


def position_from_filename(filename: str) -> str:
    """'DSC_3775-Position4-cal.fits' -> '4'.  (Lenient: the legacy tables also
    list a stray 'DSC_3770-Position4-cal-luminance.fits'.)"""
    m = re.search(r"Position(\d)", os.path.basename(filename))
    if not m:
        raise ValueError(f"Unexpected frame name: {filename}")
    return m.group(1)


def read_header(filename: str, lights_dir: Path | None = None) -> fits.Header:
    return fits.getheader(Path(lights_dir or config.CALIBRATED_LIGHTS_DIR) / filename)


def inverse_exposure(header: fits.Header) -> int:
    """EXPTIME -> nearest 1/t denominator (0.01 -> 100, 0.33333 -> 3)."""
    return int(np.round(1.0 / header["EXPTIME"]))


def luminance(rgb: np.ndarray, weights=config.LUMINANCE_WEIGHTS) -> np.ndarray:
    """Synthetic luminance L = wR*R + wG*G + wB*B from a (3, H, W) array."""
    wr, wg, wb = weights
    return (wr * rgb[0] + wg * rgb[1] + wb * rgb[2]).astype(np.float32)


def build_frame_table(
    sun_moon_centers_csv: Path | None = None,
    lights_dir: Path | None = None,
    nsigma: float = config.RADIUS_FILTER_NSIGMA,
) -> pd.DataFrame:
    """
    Frame inventory used by every downstream step (legacy:
    statistical_analysis_figures.ipynb cells 3-4).

    Reads the Sun/Moon-centre table, drops the excluded filter-wheel positions,
    attaches exposure times from the FITS headers and flags frames whose CHT
    lunar radius lies within ``nsigma`` of the mean.
    """
    csv = Path(sun_moon_centers_csv or config.SUN_MOON_CENTERS_CSV)
    t = pd.read_csv(csv)
    t["position"] = [position_from_filename(f) for f in t["filename"]]
    t = t[~t["position"].isin(config.EXCLUDED_POSITIONS)].copy()

    exp = [read_header(f, lights_dir)["EXPTIME"] for f in t["filename"]]
    t["exposure_time"] = exp
    t["inverse_exposure_time"] = np.round(1.0 / t["exposure_time"]).astype(int)

    mean_r = t["moon_radius"].mean()
    std_r = t["moon_radius"].std(ddof=0)          # np.std, as in the legacy code
    t["radius_ok"] = (t["moon_radius"] - mean_r).abs() <= nsigma * std_r
    # frames left out of the paper's stacks on top of the radius filter
    t["manual_ok"] = ~t["filename"].isin(config.EXTRA_EXCLUDED_FRAMES)
    t["use"] = t["radius_ok"] & t["manual_ok"]
    t.attrs.update(moon_radius_mean=mean_r, moon_radius_std=std_r)
    return t.reset_index(drop=True)


def paper_metadata_table(frame_table: pd.DataFrame) -> pd.DataFrame:
    """Appendix table of usable frames (legacy cell 7)."""
    df = frame_table[frame_table["radius_ok"]].copy()
    # legacy: sort by position then exposure (longest first); filename added as a
    # deterministic tie-breaker (the legacy order of ties followed os.listdir)
    df = df.sort_values(by=["position", "exposure_time", "filename"], ascending=[True, False, True])
    pos_map = {p: f"{int(a):+d}$^\\circ$".replace("+0", "0") for p, a in config.POLARIZER_ANGLE_DEG.items()}
    df["Polarization"] = df["position"].map(pos_map)
    df["Exposure"] = "1/" + df["inverse_exposure_time"].astype(int).astype(str)
    out = pd.DataFrame({
        "Filename": df["filename"],
        "Polarization": df["Polarization"],
        "Exposure": df["Exposure"],
        "Moon Radius (px)": df["moon_radius"].round(2),
        "Sun Radius (px)": df["sun_radius"].round(2),
    })
    return out.reset_index(drop=True)


# ==========================================================================
# 2. Circular Hough transform lunar-limb fit   (legacy: automated_sheet_for_centers_radius_and_SB.py)
# ==========================================================================

def find_moon_center_radius(
    img: np.ndarray,
    img_rescale_factor: int = config.CHT_RESCALE_FACTOR,
    radii_range: np.ndarray = config.CHT_RADII_RANGE_PX,
    refine_halfwidth: int = config.CHT_REFINE_HALFWIDTH_PX,
    verbose: bool = False,
):
    """
    Two-pass circular Hough transform on a 2-D image.

    Pass 1: Canny edges on a 10x downsampled image, CHT over ``radii_range``.
    Pass 2: crop 2R around the coarse centre at full resolution and refine the
            radius within +/- ``refine_halfwidth`` px in 1 px steps.
    Returns (xc, yc, radius) in full-resolution pixels.
    """
    img = img / np.max(img)
    lowres = resize(img, (img.shape[0] / img_rescale_factor, img.shape[1] / img_rescale_factor))
    edges = canny(lowres)
    hough_radii = radii_range / img_rescale_factor
    hough_res = hough_circle(edges, hough_radii)
    _, cx, cy, radii = hough_circle_peaks(hough_res, hough_radii, total_num_peaks=1)
    if verbose:
        print(f"coarse: centre=({cx[0]*img_rescale_factor}, {cy[0]*img_rescale_factor}) "
              f"radius={radii[0]*img_rescale_factor}")

    crop_dist = int(radii[0] * 2 * img_rescale_factor)
    y0 = cy[0] * img_rescale_factor - crop_dist
    x0 = cx[0] * img_rescale_factor - crop_dist
    crop = img[y0:cy[0] * img_rescale_factor + crop_dist, x0:cx[0] * img_rescale_factor + crop_dist]
    crop_edges = canny(crop)
    crop_radii = np.arange(radii[0] * img_rescale_factor - refine_halfwidth,
                           radii[0] * img_rescale_factor + refine_halfwidth, 1)
    crop_res = hough_circle(crop_edges, crop_radii)
    _, ccx, ccy, cradii = hough_circle_peaks(crop_res, crop_radii, total_num_peaks=1)
    xc = ccx[0] + x0
    yc = ccy[0] + y0
    if verbose:
        print(f"refined: centre=({xc}, {yc}) radius={cradii[0]}")
    return float(xc), float(yc), float(cradii[0])


# ==========================================================================
# 3. Skyfield Sun-centre from the Moon centre   (legacy: find_sun_center.ipynb cell 9)
# ==========================================================================

def position_angle_deg(ra1_deg, dec1_deg, ra2_deg, dec2_deg) -> float:
    """Position angle of object 2 relative to object 1, degrees E of N, [0, 360)."""
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1_deg, dec1_deg, ra2_deg, dec2_deg))
    dra = ra2 - ra1
    y = math.sin(dra) * math.cos(dec2)
    x = math.cos(dec1) * math.sin(dec2) - math.sin(dec1) * math.cos(dec2) * math.cos(dra)
    return (math.degrees(math.atan2(y, x)) + 360) % 360


def dms_to_decimal(dms_str: str) -> float:
    """"32°47'08.34'' N" -> +32.7856...;  W and S are negative."""
    s = dms_str.strip().upper()
    sign = -1.0 if ("S" in s or "W" in s) else 1.0
    parts = re.findall(r"\d+\.?\d*", re.sub(r"[NSEW]", "", s))
    if len(parts) != 3:
        raise ValueError(f"Cannot parse DMS string {dms_str!r}")
    d, m, sec = map(float, parts)
    return sign * (d + m / 60.0 + sec / 3600.0)


def fits_timestamp_utc(header: fits.Header, camera_tz: str = config.CAMERA_CLOCK_TZ) -> str:
    """
    DATE-AVG / DATE-OBS(+TIME-OBS) interpreted in the camera's clock zone and
    converted to an ISO UTC string (millisecond precision, no offset suffix).
    """
    date_obs, time_obs, date_avg = header.get("DATE-OBS"), header.get("TIME-OBS"), header.get("DATE-AVG")
    if date_avg:
        stamp = date_avg
    elif date_obs and "T" in date_obs:
        stamp = date_obs
    elif date_obs and time_obs:
        stamp = f"{date_obs.split('T')[0]}T{time_obs}"
    else:
        raise ValueError("No usable timestamp keyword in header")
    naive = date_parser.parse(stamp)
    aware = naive.replace(tzinfo=ZoneInfo(camera_tz))
    utc = aware.astimezone(dt.timezone.utc)
    return utc.isoformat(sep="T", timespec="milliseconds").replace("+00:00", "")


class Ephemeris:
    """Thin wrapper around Skyfield for the Sun/Moon geometry at the site."""

    def __init__(self, kernel: Path | None = None,
                 lat_dms: str = config.SITE_LAT_DMS, lon_dms: str = config.SITE_LON_DMS,
                 alt_m: float = config.SITE_ALT_M):
        from skyfield.api import Loader, Topos, load
        kernel = Path(kernel or config.EPHEMERIS_FILE)
        self.ts = load.timescale()
        if kernel.exists():
            self.eph = load(str(kernel))
        else:                                   # first run on a fresh clone: fetch DE421 (16 MB) into data/
            self.eph = Loader(str(kernel.parent))(kernel.name)
        self.observer = Topos(latitude_degrees=dms_to_decimal(lat_dms),
                              longitude_degrees=dms_to_decimal(lon_dms), elevation_m=alt_m)
        self.earth, self.sun, self.moon = self.eph["earth"], self.eph["sun"], self.eph["moon"]

    def _here(self, utc_iso: str):
        t = self.ts.from_astropy(Time(utc_iso, format="isot", scale="utc"))
        return (self.earth + self.observer).at(t)

    def apparent_sun_radec(self, utc_iso: str):
        """Apparent (of-date) RA, Dec of the Sun in degrees."""
        ra, dec, _ = self._here(utc_iso).observe(self.sun).apparent().radec(epoch="date")
        return ra._degrees, dec.degrees

    def star_offsets(self, utc_iso: str, stars: pd.DataFrame, body: str = "sun"):
        """
        For a table with columns name, ra_deg, dec_deg (ICRS), return a DataFrame of
        separation (arcsec) and position angle (deg E of N) of each star from the
        topocentric apparent Sun or Moon (``body``) at ``utc_iso``.
        """
        from skyfield.api import Star
        here = self._here(utc_iso)
        sun = here.observe(self.sun if body == "sun" else self.moon).apparent()
        ra_s, dec_s, _ = sun.radec(epoch="date")
        rows = []
        for _, r in stars.iterrows():
            st = here.observe(Star(ra_hours=r.ra_deg / 15.0, dec_degrees=r.dec_deg)).apparent()
            ra, dec, _ = st.radec(epoch="date")
            rows.append(dict(name=r["name"], sep_arcsec=sun.separation_from(st).degrees * 3600.0,
                             pa_deg=position_angle_deg(ra_s._degrees, dec_s.degrees, ra._degrees, dec.degrees)))
        return pd.DataFrame(rows)

    def sun_moon_geometry(self, utc_iso: str):
        """
        Returns (moon_radius_arcsec, sun_radius_arcsec, separation_arcsec,
        PA_of_sun_relative_to_moon_deg) for the given UTC time.
        """
        t = self.ts.from_astropy(Time(utc_iso, format="isot", scale="utc"))
        here = (self.earth + self.observer).at(t)
        astro_moon, astro_sun = here.observe(self.moon), here.observe(self.sun)
        app_moon, app_sun = astro_moon.apparent(), astro_sun.apparent()
        ra_m, dec_m, _ = app_moon.radec()
        ra_s, dec_s, _ = app_sun.radec()
        pa = position_angle_deg(ra_m.degrees, dec_m.degrees, ra_s.degrees, dec_s.degrees)
        sep_arcsec = app_moon.separation_from(app_sun).degrees * 3600.0
        rad2as = (180.0 / math.pi) * 3600.0
        moon_r = math.asin(config.MOON_RADIUS_KM / astro_moon.distance().km) * rad2as
        sun_r = math.asin(config.SUN_RADIUS_KM / astro_sun.distance().km) * rad2as
        return moon_r, sun_r, sep_arcsec, pa


def north_east_vectors(alpha_deg: float = config.CELESTIAL_NORTH_DEG, hand: int = 1):
    """
    Unit vectors of celestial North and East in ARRAY coordinates (x right, y down)
    for a North direction ``alpha_deg`` CCW from +x as displayed in camera orientation.
    ``hand`` = -1 mirrors East (only used when fitting the handedness in notebook 03).
    """
    a = math.radians(alpha_deg)
    n = np.array([math.cos(a), -math.sin(a)])
    e = np.array([-math.sin(a), -math.cos(a)]) * hand
    return n, e


def sky_offset_to_pixels(sep_arcsec, pa_deg, px_per_arcsec, alpha_deg=config.CELESTIAL_NORTH_DEG, hand=1):
    """(separation, PA east of north) -> (dx, dy) pixel offset in array coordinates."""
    n, e = north_east_vectors(alpha_deg, hand)
    pa = math.radians(pa_deg)
    v = sep_arcsec * px_per_arcsec * (math.sin(pa) * e + math.cos(pa) * n)
    return float(v[0]), float(v[1])


def sun_center_from_moon(moon_r_arcsec, sep_arcsec, pa_deg,
                         moon_xc_px, moon_yc_px, moon_r_px,
                         convention: str = config.SUN_CENTER_ROTATION,
                         alpha_deg: float = config.CELESTIAL_NORTH_DEG,
                         image_north_from_y_deg: float = config.IMAGE_NORTH_ANGLE_FROM_Y_DEG):
    """
    Sun-disk centre in pixels from the Moon centre, the plate scale implied by
    the Moon's pixel radius, and the Sun-Moon displacement vector rotated into
    the image frame.  Returns (sun_xc, sun_yc, px_per_arcsec).

    convention="astrometric": :func:`sky_offset_to_pixels` (validated on stars).
    convention="legacy":      dx = sep*sin(PA-78), dy = sep*cos(PA-78) as in the
                              legacy find_sun_center notebook (see config).
    """
    if moon_r_px <= 0 or moon_r_arcsec <= 0:
        raise ValueError("Moon radius must be positive")
    px_per_arcsec = moon_r_px / moon_r_arcsec
    if convention == "legacy":
        pa_img = math.radians(pa_deg - image_north_from_y_deg)
        dx = sep_arcsec * px_per_arcsec * math.sin(pa_img)
        dy = sep_arcsec * px_per_arcsec * math.cos(pa_img)
    elif convention == "astrometric":
        dx, dy = sky_offset_to_pixels(sep_arcsec, pa_deg, px_per_arcsec, alpha_deg)
    else:
        raise ValueError(f"unknown convention {convention!r}")
    return moon_xc_px + dx, moon_yc_px + dy, px_per_arcsec


# ==========================================================================
# 3b. Background stars: detection, PSF fits, orientation fit   (notebook 03)
# ==========================================================================

def reference_star_positions(cache: Path | None = None, radius_deg: float = 1.8, vmax: float = 8.5,
                             utc_iso: str = config.TOTALITY_MID_UTC) -> pd.DataFrame:
    """
    Catalogue stars (Simbad, V < ``vmax``) within ``radius_deg`` of the Sun at
    ``utc_iso``: columns name, ra_deg, dec_deg (ICRS), vmag.  Cached in
    data/reference_stars.csv so the notebooks work offline afterwards.
    """
    cache = Path(cache or config.REFERENCE_STARS_CSV)
    if cache.exists():
        return pd.read_csv(cache)
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.simbad import Simbad
    ra, dec = Ephemeris().apparent_sun_radec(utc_iso)
    sim = Simbad()
    sim.add_votable_fields("V", "otype")
    t = sim.query_region(SkyCoord(ra=ra * u.deg, dec=dec * u.deg), radius=radius_deg * u.deg).to_pandas()
    vcol = [c for c in t.columns if c.lower() in ("v", "flux_v")][0]
    t = t[t[vcol] < vmax].sort_values(vcol)
    df = pd.DataFrame(dict(name=t["main_id"].str.replace(r"\s+", " ", regex=True).str.strip("* ").str.strip(),
                           ra_deg=t["ra"].astype(float), dec_deg=t["dec"].astype(float), vmag=t[vcol].astype(float)))
    df.to_csv(cache, index=False)
    return df.reset_index(drop=True)


def detect_point_sources(img, center_xy=config.SUN_CENTER_XY, min_radius=config.STAR_SEARCH_MIN_RADIUS_PX,
                         hp_size=21, smooth=1.2, nsigma=10.0, box=31, max_sources=10):
    """
    Compact sources in ``img`` outside ``min_radius`` from ``center_xy``:
    median high-pass, light Gaussian smoothing, local maxima above ``nsigma``
    robust noise.  Returns a DataFrame (x, y, snr) sorted by SNR.
    """
    img = np.asarray(img, dtype=np.float32)
    yy, xx = np.mgrid[:img.shape[0], :img.shape[1]]
    r = np.hypot(xx - center_xy[0], yy - center_xy[1])
    hp = img - ndimage.median_filter(img, size=hp_size)
    sm = ndimage.gaussian_filter(hp, smooth)
    valid = (r > min_radius) & (img > 0)
    noise = 1.4826 * np.median(np.abs(sm[valid & (r > min_radius + 1000)]))
    peaks = (sm == ndimage.maximum_filter(sm, size=box)) & valid & (sm > nsigma * noise)
    ys, xs = np.nonzero(peaks)
    df = pd.DataFrame(dict(x=xs, y=ys, snr=sm[ys, xs] / noise)).sort_values("snr", ascending=False)
    df.attrs["noise"] = float(noise)
    return df.head(max_sources).reset_index(drop=True), hp


def _gauss2d(coords, amp, x0, y0, sx, sy, bg):
    x, y = coords
    return (amp * np.exp(-0.5 * (((x - x0) / sx) ** 2 + ((y - y0) / sy) ** 2)) + bg).ravel()


def _gauss2d_double(coords, amp1, x1, y1, amp2, x2, y2, s, bg):
    x, y = coords
    g = amp1 * np.exp(-0.5 * (((x - x1) / s) ** 2 + ((y - y1) / s) ** 2))
    g += amp2 * np.exp(-0.5 * (((x - x2) / s) ** 2 + ((y - y2) / s) ** 2))
    return (g + bg).ravel()


def fit_gaussian_2d(cutout, x0=None, y0=None, double=False, sep_px=None):
    """
    Fit a 2-D Gaussian (or two Gaussians of common width, initial separation
    ``sep_px``) to ``cutout``.  Returns a dict with centroid(s), sigma, FWHM,
    amplitude, background, model image and fit uncertainties (1-sigma).
    """
    from scipy.optimize import curve_fit
    cut = np.asarray(cutout, dtype=np.float64)
    h, w = cut.shape
    yy, xx = np.mgrid[:h, :w]
    bg0 = np.median(cut)
    amp0 = cut.max() - bg0
    x0 = w / 2 if x0 is None else x0
    y0 = h / 2 if y0 is None else y0
    if not double:
        p0 = [amp0, x0, y0, 3.0, 3.0, bg0]
        popt, pcov = curve_fit(_gauss2d, (xx, yy), cut.ravel(), p0=p0, maxfev=20000)
        err = np.sqrt(np.diag(pcov))
        sigma = float(np.hypot(popt[3], popt[4]) / np.sqrt(2))
        return dict(x=float(popt[1]), y=float(popt[2]), x_err=float(err[1]), y_err=float(err[2]),
                    sigma=sigma, fwhm=2.3548 * sigma, amp=float(popt[0]), bg=float(popt[5]),
                    model=_gauss2d((xx, yy), *popt).reshape(h, w), components=1)
    sep = sep_px or 15.0
    p0 = [amp0, x0 - sep / 4, y0, amp0 * 0.7, x0 + sep / 4, y0, 3.0, bg0]
    popt, pcov = curve_fit(_gauss2d_double, (xx, yy), cut.ravel(), p0=p0, maxfev=40000)
    err = np.sqrt(np.diag(pcov))
    a1, x1, y1, a2, x2, y2, sg, bg = popt
    xc = (a1 * x1 + a2 * x2) / (a1 + a2)          # flux-weighted centroid of the blend
    yc = (a1 * y1 + a2 * y2) / (a1 + a2)
    return dict(x=float(xc), y=float(yc), x_err=float(np.hypot(err[1], err[4]) / 2), y_err=float(np.hypot(err[2], err[5]) / 2),
                x1=float(x1), y1=float(y1), x2=float(x2), y2=float(y2), sep_px=float(np.hypot(x2 - x1, y2 - y1)),
                sigma=float(sg), fwhm=2.3548 * float(sg), amp=float(max(a1, a2)), bg=float(bg),
                model=_gauss2d_double((xx, yy), *popt).reshape(h, w), components=2)


def fit_image_rotation_to_stars(measured_xy, sky_offsets, center_xy,
                                alpha0=config.CELESTIAL_NORTH_DEG, scale0=config.PLATE_SCALE_ARCSEC_PER_PIX):
    """
    Least-squares fit of the North angle (deg) and plate scale (arcsec/px) — for
    both handedness values — that map the (sep, PA) sky offsets from a reference
    body (Sun or Moon, whose pixel position is ``center_xy``) onto the measured
    star pixel positions.
    Returns the best solution as a dict (alpha_deg, alpha_err, scale, scale_err,
    hand, residuals (N,2) px, rms_px).
    """
    from scipy.optimize import least_squares
    meas = np.asarray(measured_xy, dtype=float)
    seps = np.asarray(sky_offsets["sep_arcsec"], dtype=float)
    pas = np.asarray(sky_offsets["pa_deg"], dtype=float)
    c = np.asarray(center_xy, dtype=float)

    def model(p, hand):
        alpha, scale = p
        return np.array([c + np.array(sky_offset_to_pixels(s, pa, 1.0 / scale, alpha, hand)) for s, pa in zip(seps, pas)])

    best = None
    for hand in (1, -1):
        res = least_squares(lambda p: (model(p, hand) - meas).ravel(), x0=[alpha0, scale0])
        resid = model(res.x, hand) - meas
        rms = float(np.sqrt(np.mean(resid ** 2)))
        dof = max(resid.size - 2, 1)
        cov = np.linalg.inv(res.jac.T @ res.jac) * (np.sum(resid ** 2) / dof)
        sol = dict(alpha_deg=float(res.x[0]) % 360, alpha_err=float(np.sqrt(cov[0, 0])),
                   scale=float(res.x[1]), scale_err=float(np.sqrt(cov[1, 1])), hand=hand,
                   residuals=resid, rms_px=rms)
        if best is None or rms < best["rms_px"]:
            best = sol
    return best


def measure_interchannel_shift(ref, img, center_xy=config.SUN_CENTER_XY, radius_rsun=1.4, rsun_px=None,
                               half=400, n_patches=12, dog_sigma=4.0, upsample=20):
    """
    Residual registration (dy, dx) of ``img`` relative to ``ref``: plain cross-
    correlation (upsampled) of difference-of-Gaussians filtered, Hann-windowed
    patches placed around the Sun at ``radius_rsun``; the median over patches is
    returned together with the 75th percentile of the patch scatter.
    (A masked annulus does not work: the mask edges dominate the correlation.)
    """
    from skimage.registration import phase_cross_correlation
    rsun_px = rsun_px or (959.0 / config.PLATE_SCALE_ARCSEC_PER_PIX)

    def dog(x):
        x = np.asarray(x, dtype=np.float32)
        hp = ndimage.gaussian_filter(x, 1.5) - ndimage.gaussian_filter(x, dog_sigma)
        return (hp / (ndimage.gaussian_filter(np.abs(hp), 3 * dog_sigma) + 1e-9)).astype(np.float32)

    A, B = dog(ref), dog(img)
    w = np.outer(np.hanning(2 * half), np.hanning(2 * half))
    out = []
    for ang in np.linspace(0, 360, n_patches, endpoint=False):
        x0 = int(center_xy[0] + radius_rsun * rsun_px * np.cos(np.radians(ang)))
        y0 = int(center_xy[1] + radius_rsun * rsun_px * np.sin(np.radians(ang)))
        if y0 - half < 0 or x0 - half < 0 or y0 + half > A.shape[0] or x0 + half > A.shape[1]:
            continue
        pa = A[y0 - half:y0 + half, x0 - half:x0 + half] * w
        pb = B[y0 - half:y0 + half, x0 - half:x0 + half] * w
        sh, _, _ = phase_cross_correlation(pa, pb, upsample_factor=upsample, normalization=None)
        out.append(sh)
    out = np.array(out)
    med = np.median(out, axis=0)
    spread = np.percentile(np.abs(out - med), 75, axis=0)
    return float(med[0]), float(med[1]), float(np.hypot(*spread))


def exposure_ratio_table(planes_by_exp: dict) -> pd.DataFrame:
    """max(short)/max(next longer) per step of the bracket; ~0.5 per doubling = linear, ~1 = saturated."""
    exps = list(planes_by_exp)
    rows = []
    for a, b in zip(exps[:-1], exps[1:]):
        rows.append(dict(step=f"1/{a}->1/{b}", expected=b / a,
                         ratio=float(np.max(planes_by_exp[a]) / np.max(planes_by_exp[b]))))
    return pd.DataFrame(rows)


# ==========================================================================
# 4. Alignment and stacking   (step 04; recipe recovered from the legacy products,
#    confirmed by C. Gandhi: translate to the Sun centre, take the mean)
# ==========================================================================

def shift_to_center(img2d: np.ndarray, xc: float, yc: float,
                    target_xy=config.SUN_CENTER_XY,
                    order: int = config.STACK_INTERPOLATION_ORDER) -> np.ndarray:
    """Translate ``img2d`` so that pixel (xc, yc) lands on ``target_xy``; zero fill."""
    tx, ty = target_xy
    return ndimage.shift(img2d, (ty - yc, tx - xc), order=order, mode="constant", cval=0.0)


def stack_frames(filenames, centers: pd.DataFrame, lights_dir: Path | None = None,
                 dtype=np.float32) -> np.ndarray:
    """
    Mean of the Sun-centred luminance images of ``filenames``.
    ``centers`` must be indexed by filename with columns sun_xc, sun_yc.
    """
    lights_dir = Path(lights_dir or config.CALIBRATED_LIGHTS_DIR)
    acc = None
    for f in filenames:
        rgb = fits.getdata(lights_dir / f).astype(np.float32)
        lum = luminance(rgb)
        del rgb
        row = centers.loc[f]
        shifted = shift_to_center(lum, row["sun_xc"], row["sun_yc"]).astype(dtype)
        acc = shifted if acc is None else acc + shifted
    return acc / len(filenames)


def stacked_filename(inv_exp: int) -> Path:
    return config.STACKED_DIR / f"stacked_exp{inv_exp}.fits"


def load_stacked(inv_exp: int, apply_channel_shifts: bool = True) -> np.ndarray:
    """
    (3, H, W) float array [pol1, pol2, pol3] for exposure 1/inv_exp s, with the
    empirical inter-polariser shifts applied (as the legacy HDR notebooks did at
    load time).
    """
    d = fits.getdata(stacked_filename(inv_exp)).astype(np.float64)
    if apply_channel_shifts:
        for i, pos in enumerate(config.POLARIZER_POSITIONS):
            d[i] = ndimage.shift(d[i], config.CHANNEL_SHIFTS[pos], order=1)
    return d


def hdr_filename(method: str, pos: str) -> Path:
    return config.HDR_DIR / f"hdr_{method}_pol{pos}.fits"


def hdr_header(template: fits.Header, pos: str, method: str) -> fits.Header:
    """2-D header for one polariser plane derived from a stacked-image header."""
    h = template.copy()
    h["NAXIS"] = 2
    if "NAXIS3" in h:
        h.remove("NAXIS3")
    h["FILT-1"] = f"Position{pos}"
    h["HDRMETH"] = (method, "HDR construction method")
    return h


# ==========================================================================
# 5. HDR construction
# ==========================================================================

def exposure_normalization_stacking(images, exposure_times) -> np.ndarray:
    """Mean over exposures of clip(img, 0) / t   (paper eqs. 2-3)."""
    normalized = [np.clip(img.astype(np.float64), 0.0, None) / t
                  for img, t in zip(images, exposure_times)]
    return np.mean(normalized, axis=0)


# ---- LDIC (Druckmullerova 2014), legacy make_ldic_hdr_from_stacked_exposures ----

def get_polar_coordinates(height, width, center_x=None, center_y=None):
    if center_x is None:
        center_x = width / 2.0 - 0.5
    if center_y is None:
        center_y = height / 2.0 - 0.5
    x, y = np.meshgrid(np.arange(width), np.arange(height))
    xr, yr = x - center_x, y - center_y
    return np.sqrt(xr ** 2 + yr ** 2), np.arctan2(yr, xr)


def weight_function_ldic(pixel_values, low_reject_abs, low_full_abs,
                         high_start_reject_abs, high_reject_abs):
    """LDIC weight w(f): 0 below low_reject / above high_reject, 1 on the plateau, linear ramps between."""
    w = np.zeros_like(pixel_values, dtype=np.float64)
    m = (pixel_values >= low_reject_abs) & (pixel_values < low_full_abs)
    if low_full_abs > low_reject_abs:
        w[m] = (pixel_values[m] - low_reject_abs) / (low_full_abs - low_reject_abs)
    else:
        w[(pixel_values >= low_reject_abs) & (pixel_values >= low_full_abs)] = 1.0
    w[(pixel_values >= low_full_abs) & (pixel_values < high_start_reject_abs)] = 1.0
    m = (pixel_values >= high_start_reject_abs) & (pixel_values < high_reject_abs)
    if high_reject_abs > high_start_reject_abs:
        w[m] = 1.0 - (pixel_values[m] - high_start_reject_abs) / (high_reject_abs - high_start_reject_abs)
    w[pixel_values < low_reject_abs] = 0.0
    w[pixel_values >= high_reject_abs] = 0.0
    return np.clip(w, 0.0, 1.0)


def fit_trigonometric_polynomial(angles_rad, values, order):
    """Least-squares y = c0 + sum_p (c_{2p-1} cos p*phi + c_{2p} sin p*phi), order reduced adaptively."""
    n = len(angles_rad)
    if n == 0:
        return np.zeros(1 + 2 * order)
    eff = min(order, max(0, (n - 1) // 2))
    A = np.ones((n, 1 + 2 * eff))
    for p in range(1, eff + 1):
        A[:, 2 * p - 1] = np.cos(p * angles_rad)
        A[:, 2 * p] = np.sin(p * angles_rad)
    try:
        c, *_ = np.linalg.lstsq(A, values, rcond=None)
    except np.linalg.LinAlgError:
        c, eff = np.array([np.mean(values)]), 0
    full = np.zeros(1 + 2 * order)
    full[:1 + 2 * eff] = c
    return full


def evaluate_trigonometric_polynomial(coeffs, order, phi_map):
    out = np.full_like(phi_map, coeffs[0], dtype=np.float64)
    if order > 0 and len(coeffs) == 1 + 2 * order:
        for p in range(1, order + 1):
            out += coeffs[2 * p - 1] * np.cos(p * phi_map) + coeffs[2 * p] * np.sin(p * phi_map)
    return out


def prepare_images_for_ldic(images, target_max=config.LDIC_TARGET_MAX):
    """Clip negatives and apply ONE global scale so the brightest pixel of the set == target_max."""
    clipped = [np.clip(img.astype(np.float64), 0.0, None) for img in images]
    gmax = max(img.max() for img in clipped)
    if gmax <= 0:
        return clipped, target_max
    s = target_max / gmax
    return [np.clip(img * s, 0.0, target_max) for img in clipped], target_max


def ldic_hdr_stacking(images_in, exposure_times, images_ref=None,
                      max_pixel_value_input=config.LDIC_TARGET_MAX,
                      num_angular_segments=60, trig_poly_order_k=4, trig_poly_order_q=4,
                      wf_low_reject_percent=1.0, wf_low_full_percent=5.0,
                      wf_high_start_reject_percent=97.0, wf_high_reject_percent=99.9,
                      regression_weight_threshold=0.8, min_valid_segments=15,
                      sun_center_x=None, sun_center_y=None, fit_intercept=False,
                      verbose=False, return_diagnostics=False):
    """
    Linear Digital Image Composer (paper Section 3.3.2, eqs. 4-5).

    ``images_ref`` (same length as ``images_in``) is the master luminance used
    to evaluate the weight function so that all polariser channels hand over
    between exposures at the same pixels.  Longest exposure first is index 0.
    """
    from scipy.stats import linregress

    images = [img.astype(np.float64) for img in images_in]
    H, W = images[0].shape
    refs = images if images_ref is None else [r.astype(np.float64) for r in images_ref]

    order = sorted(zip(exposure_times, images, refs), key=lambda x: x[0], reverse=True)
    images = [p[1] for p in order]
    refs = [p[2] for p in order]

    g_cum = np.zeros((H, W)); w_cum = np.zeros((H, W))
    _, phi = get_polar_coordinates(H, W, sun_center_x, sun_center_y)
    phi_pos = np.where(phi < 0, phi + 2 * np.pi, phi)

    lo_rej = wf_low_reject_percent / 100 * max_pixel_value_input
    lo_full = wf_low_full_percent / 100 * max_pixel_value_input
    hi_start = wf_high_start_reject_percent / 100 * max_pixel_value_input
    hi_rej = wf_high_reject_percent / 100 * max_pixel_value_input

    step = 2 * np.pi / num_angular_segments
    seg_centers = np.arange(num_angular_segments) * step + step / 2
    diag = []

    for i, (f, fref) in enumerate(zip(images, refs)):
        cur_lo_rej, cur_lo_full = (0.0, 0.0) if i == 0 else (lo_rej, lo_full)
        shortest = i == len(images) - 1
        cur_hi_start = max_pixel_value_input if shortest else hi_start
        cur_hi_rej = max_pixel_value_input if shortest else hi_rej
        w = weight_function_ldic(fref, cur_lo_rej, cur_lo_full, cur_hi_start, cur_hi_rej)

        with np.errstate(divide="ignore", invalid="ignore"):
            g_avg = np.where(w_cum > 0, g_cum / w_cum, 0.0)

        if i == 0:
            k_c = np.zeros(1 + 2 * trig_poly_order_k); k_c[0] = 1.0
            q_c = np.zeros(1 + 2 * trig_poly_order_q)
        else:
            ks, qs, angs = [], [], []
            for s in range(num_angular_segments):
                a0, a1 = s * step, (s + 1) * step
                if s == num_angular_segments - 1:
                    m = (phi_pos >= a0) & (phi_pos <= 2 * np.pi)
                else:
                    m = (phi_pos >= a0) & (phi_pos < a1)
                if not np.any(m):
                    continue
                X, Y, Wm = f[m], g_avg[m], w[m]
                valid = (Wm >= regression_weight_threshold) & (Y > lo_full)
                X, Y = X[valid], Y[valid]
                if len(X) < 2 or np.all(X == X[0]):
                    continue
                if fit_intercept:
                    k, q, *_ = linregress(X, Y)
                else:
                    k, q = np.sum(X * Y) / np.sum(X ** 2), 0.0   # sky-subtracted: fit through origin
                if k <= 0 or not (np.isfinite(k) and np.isfinite(q)):
                    continue
                ks.append(k); qs.append(q); angs.append(seg_centers[s])
            if verbose:
                print(f"  exposure {i}: {len(angs)}/{num_angular_segments} valid segments")
            if len(angs) < min_valid_segments:
                continue
            angs = np.array(angs)
            k_c = fit_trigonometric_polynomial(angs, np.array(ks), trig_poly_order_k)
            q_c = fit_trigonometric_polynomial(angs, np.array(qs), trig_poly_order_q)

        k_map = evaluate_trigonometric_polynomial(k_c, trig_poly_order_k, phi)
        q_map = evaluate_trigonometric_polynomial(q_c, trig_poly_order_q, phi)
        g_cum += w * (k_map * f + q_map)
        w_cum += w
        if return_diagnostics:
            diag.append(dict(exposure_time=order[i][0], k_coeffs=k_c.copy(),
                             k_at=evaluate_trigonometric_polynomial(k_c, trig_poly_order_k, seg_centers),
                             segment_angles=seg_centers.copy(),
                             segment_k=(np.array(ks) if i > 0 else np.ones(0)),
                             segment_k_angles=(np.array(angs) if i > 0 and len(angs) else np.ones(0)),
                             weight_fraction=float(np.mean(w > 0)), n_valid_segments=(len(angs) if i > 0 else num_angular_segments)))

    with np.errstate(divide="ignore", invalid="ignore"):
        g = np.where(w_cum > 0, g_cum / w_cum, 0.0)
    return (g, diag) if return_diagnostics else g


# ==========================================================================
# 6. Polarimetry
# ==========================================================================

def stokes_to_mzp(I, Q, U, phi_deg):
    """Stokes I,Q,U -> intensities behind polarisers at phi, phi+60, phi-60 (degrees)."""
    phi = np.radians(phi_deg)
    angles = [phi, phi + np.radians(60), phi - np.radians(60)]
    return tuple(0.5 * (I + Q * np.cos(2 * th) + U * np.sin(2 * th)) for th in angles)


def mzp_to_stokes(I_phi, I_phi_p60, I_phi_m60, phi_deg):
    """Inverse of :func:`stokes_to_mzp`."""
    phi = np.radians(phi_deg)
    c, s = np.cos(2 * phi), np.sin(2 * phi)
    I = (2 / 3) * (I_phi + I_phi_p60 + I_phi_m60)
    delta = I_phi - 0.5 * (I_phi_p60 + I_phi_m60)
    diff = I_phi_p60 - I_phi_m60
    Q = (4 / 3) * (delta * c - (np.sqrt(3) / 2) * diff * s)
    U = (4 / 3) * (delta * s + (np.sqrt(3) / 2) * diff * c)
    return I, Q, U


def stokes_to_polarization(I, Q, U):
    """Degree and angle (deg) of linear polarisation."""
    return np.sqrt(Q ** 2 + U ** 2) / I, np.degrees(0.5 * np.arctan2(U, Q))


def rotate_to_solar_frame(pol1, pol2, pol3, solar_frame_phi_deg=config.SOLAR_FRAME_PHI_DEG):
    """
    Re-express the instrument-frame triplet (0, +60, -60 deg) as the triplet an
    observer would measure with the reference polariser aligned to solar North.
    Returns a (3, H, W) stack.
    """
    I, Q, U = mzp_to_stokes(pol1, pol2, pol3, 0.0)
    return np.stack(stokes_to_mzp(I, Q, U, solar_frame_phi_deg), axis=0)


# ==========================================================================
# 7. Display composites, MGN, annotations   (legacy: *_hdr_paper_figure_Paras.ipynb)
# ==========================================================================

def make_hdr_composite(stack, lower_pct=config.COMPOSITE_LOWER_PCT, upper_pct=config.COMPOSITE_UPPER_PCT,
                       gamma=config.COMPOSITE_GAMMA, boost_pct=config.COMPOSITE_BOOST_PCT, eps=1e-8):
    """
    Colour = per-channel polarisation fraction (ch / total);
    luminance = percentile-clipped total intensity, gamma-stretched.
    """
    total = stack[0] + stack[1] + stack[2] + eps
    R, G, B = stack[0] / total, stack[1] / total, stack[2] / total
    p_low, p_high = np.percentile(total, [lower_pct, upper_pct])
    L = np.clip(total, p_low, p_high)
    L = ((L - p_low) / (p_high - p_low)) ** gamma
    RGB = np.stack([R * L, G * L, B * L], axis=-1)
    RGB = RGB / np.percentile(RGB, boost_pct)
    return np.clip(RGB, 0, 1)


def make_mgn_polarimetric_composite(stack, sigma=config.MGN_SIGMAS, gamma=config.COMPOSITE_GAMMA,
                                    boost_pct=config.COMPOSITE_BOOST_PCT):
    """MGN applied to the total linear intensity, recombined with the polarisation-fraction colours."""
    from sunkit_image.enhance import mgn
    total = stack[0] + stack[1] + stack[2] + 1e-8
    R, G, B = stack[0] / total, stack[1] / total, stack[2] / total
    L = np.clip(mgn(total / np.max(total), sigma=list(sigma)), 0, 1) ** gamma
    RGB = np.stack([R * L, G * L, B * L], axis=-1)
    RGB = RGB / np.percentile(RGB, boost_pct)
    return np.clip(RGB, 0, 1)


def asinh_stretch(img, black_point_percentile=1.0, scale_percentile=99.5):
    """Asinh stretch for quick-look display of a single HDR channel."""
    img = np.clip(img, 0, None)
    nz = img[img > 0]
    if len(nz) == 0:
        return img
    black, scale = np.percentile(nz, black_point_percentile), np.percentile(nz, scale_percentile)
    norm = np.clip(img - black, 0, None) / (scale - black + 1e-30)
    return np.clip(np.arcsinh(10.0 * norm) / np.arcsinh(10.0), 0, 1)


def normalise_channel(arr):
    """1-99 percentile stretch of the positive pixels, for quick-look RGB previews."""
    lo, hi = np.percentile(arr[arr > 0], [1, 99]) if (arr > 0).any() else (0, 1)
    return np.clip((arr - lo) / (hi - lo + 1e-30), 0, 1)


def draw_compass_arrow(ax, center, angle_deg, length, label, color="white"):
    a = -np.deg2rad(angle_deg)
    ex, ey = center[0] + length * np.cos(a), center[1] + length * np.sin(a)
    ax.annotate("", xy=(ex, ey), xytext=center, arrowprops=dict(arrowstyle="->", color=color, lw=1.5))
    ax.text(center[0] + length * 1.2 * np.cos(a), center[1] + length * 1.2 * np.sin(a), label,
            color=color, fontsize=14, fontweight="bold", ha="center", va="center")


def add_discrete_color_wheel(parent_ax, align_angle_deg, image_gamma=config.COMPOSITE_GAMMA,
                             corona_saturation=config.COLORWHEEL_SATURATION):
    """Inset polar axes showing which colour a given polarisation angle maps to."""
    polar_class = get_projection_class("polar")
    ax = inset_axes(parent_ax, width="22%", height="22%", loc="upper left",
                    bbox_to_anchor=(0.0, -0.05, 1, 1), bbox_transform=parent_ax.transAxes,
                    axes_class=polar_class, borderpad=2)
    ax.set_theta_offset(np.radians(align_angle_deg))
    ax.set_theta_direction(1)
    ax.axis("off")
    for angle in np.arange(0, 360, 30):
        t = np.radians(angle)
        rgb = np.array([0.5 * (1 + np.cos(2 * (-t))),
                        0.5 * (1 + np.cos(2 * (-t - np.radians(60)))),
                        0.5 * (1 + np.cos(2 * (-t + np.radians(60))))])
        rgb = (rgb / rgb.max()) ** image_gamma
        hsv = mcolors.rgb_to_hsv(rgb); hsv[1] *= corona_saturation
        rgb = mcolors.hsv_to_rgb(hsv)
        ax.annotate("", xy=(t, 1.05), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="->", color=rgb, lw=2.5, shrinkA=0, shrinkB=0, mutation_scale=10))
        ax.text(t, 1.35, f"{angle}°", color=rgb, ha="center", va="center", fontsize=8, fontweight="bold")
    ax.plot(0, 0, "o", color="white", markersize=8, zorder=5)
    ax.set_ylim(0, 1.25)


def plot_solar_frame_panels(hdr_plot, mgn_plot, titles, out_png, moon_mask=False,
                            center_xy=config.SUN_CENTER_XY,
                            solar_north_angle=config.SOLAR_NORTH_ANGLE_DEG, dpi=config.FIGURE_DPI):
    """The two-panel paper figure (HDR composite over its MGN-enhanced version)."""
    import matplotlib.pyplot as plt
    cx, cy = center_xy
    compass = (cx + 2000, cy + 1500)
    fig, ax = plt.subplots(2, 1, figsize=(8, 11), constrained_layout=True)
    for i, (a, img, title) in enumerate(zip(ax, [hdr_plot, mgn_plot], titles)):
        a.imshow(img)
        a.set_title(title, fontsize=14, fontweight="bold", pad=10)
        a.axis("off")
        a.set_xlim(cx - 3000, cx + 3000)
        a.set_ylim(cy + 2000, cy - 1500)
        if moon_mask:
            a.add_artist(plt.Circle((cx, cy), config.MOON_MASK_RADIUS_PX, color="black", fill=True, zorder=5))
        draw_compass_arrow(a, compass, solar_north_angle, 500, "N")
        draw_compass_arrow(a, compass, solar_north_angle - 90, 500, "W")
        a.plot(compass[0], compass[1], "o", color="white", markersize=8, zorder=5)
        if i == 0:
            add_discrete_color_wheel(a, solar_north_angle)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    return fig


def plot_exposure_histograms(hists, counts, orientation="vertical", out_stem=None, dpi=config.FIGURE_DPI):
    """
    Fig. 8 of the paper.  ``hists[pos]`` = list of normalised-intensity samples, one
    per column (nine exposures, then LDIC HDR, then Exp-Norm HDR); ``counts[pos][e]`` =
    frames stacked.  orientation="vertical" is the paper layout (intensity on y,
    counts on x, 11 panels in a row); "horizontal" is the conventional layout
    (intensity on x, log counts on y) in a 4 x 3 grid.
    """
    import matplotlib.pyplot as plt
    colors = {"1": "red", "2": "green", "3": "blue"}
    labels = {p: f"Pol {p} ({config.POLARIZER_ANGLE_DEG[p]:+.0f}°)".replace("+0°", "0°") for p in config.POLARIZER_POSITIONS}
    exps = list(config.INV_EXPOSURES)
    titles = [f"1/{e} s" for e in exps] + ["LDIC HDR", "Exp. Norm HDR"]
    n = len(titles)
    if orientation == "vertical":
        fig, axes = plt.subplots(1, n, figsize=(24, 10), sharey=True)
        axes = list(axes)
    else:
        fig, grid = plt.subplots(3, 4, figsize=(20, 12), sharex=True, sharey=True)
        axes = list(grid.flat)
    for i, ax in enumerate(axes[:n]):
        for pos in config.POLARIZER_POSITIONS:
            ax.hist(hists[pos][i], bins=150, color=colors[pos], alpha=0.75, histtype="step", linewidth=2.0,
                    label=labels[pos] if i == 0 else "",
                    orientation="horizontal" if orientation == "vertical" else "vertical")
        if orientation == "vertical":
            ax.set_title(titles[i], fontsize=17, fontweight="bold", pad=35)
            if i < len(exps):
                for x, pos in zip((0.15, 0.50, 0.85), config.POLARIZER_POSITIONS):
                    ax.text(x, 1.02, f"{counts[pos][exps[i]]}", color=colors[pos], transform=ax.transAxes,
                            ha="center", va="bottom", fontsize=16, fontweight="bold")
            ax.set_xscale("log"); ax.set_ylim(0, 1.05)
            ax.set_xlabel("Counts", fontsize=18, fontweight="bold")
            if i == 0:
                ax.set_ylabel("Normalized Pixel Intensity", fontsize=18, fontweight="bold", labelpad=12)
            ax.grid(True, axis="y", alpha=0.4, linestyle="--", linewidth=1.0)
        else:
            ax.set_title(titles[i], fontsize=16, fontweight="bold", pad=8)
            if i < len(exps):
                for k, pos in enumerate(config.POLARIZER_POSITIONS):
                    ax.text(0.97, 0.95 - 0.09 * k, f"{counts[pos][exps[i]]} frames", color=colors[pos],
                            transform=ax.transAxes, ha="right", va="top", fontsize=13, fontweight="bold")
            ax.set_yscale("log"); ax.set_xlim(0, 1.05)
            ax.grid(True, axis="x", alpha=0.4, linestyle="--", linewidth=1.0)
        for sp in ax.spines.values():
            sp.set_linewidth(1.5)
        ax.tick_params(axis="both", which="major", labelsize=14, width=1.5, length=6)
        ax.tick_params(axis="both", which="minor", width=1.0, length=3)
    if orientation == "vertical":
        axes[0].legend(loc="upper right", fontsize=14, framealpha=0.95, edgecolor="black", fancybox=False)
        plt.tight_layout(); plt.subplots_adjust(wspace=0.0)
    else:
        handles, lbls = axes[0].get_legend_handles_labels()
        leg_ax = axes[n]
        leg_ax.axis("off")
        leg_ax.legend(handles, lbls, loc="center", fontsize=16, framealpha=0.95, edgecolor="black", fancybox=False)
        for ax in grid[-1, :]:
            ax.set_xlabel("Normalized Pixel Intensity", fontsize=16, fontweight="bold")
        for ax in grid[:, 0]:
            ax.set_ylabel("Counts", fontsize=16, fontweight="bold")
        plt.tight_layout()
    if out_stem is not None:
        fig.savefig(f"{out_stem}.png", bbox_inches="tight", dpi=dpi)
        fig.savefig(f"{out_stem}.pdf", bbox_inches="tight")
    return fig


# ==========================================================================
# 8. Validation helpers
# ==========================================================================

def compare_arrays(a: np.ndarray, b: np.ndarray, region=None) -> dict:
    """Correlation and residual statistics between two same-shaped arrays."""
    if region is not None:
        a, b = a[region], b[region]
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    diff = np.abs(a - b)
    scale = max(np.abs(b).max(), 1e-30)
    return dict(corr=float(np.corrcoef(a, b)[0, 1]),
                median_abs_diff=float(np.median(diff)),
                p99_abs_diff=float(np.percentile(diff, 99)),
                max_abs_diff=float(diff.max()),
                max_rel_to_peak=float(diff.max() / scale))


def compare_images(png_a: Path, png_b: Path, size=(300, 400)) -> float:
    """Pearson correlation of two PNGs after resizing both to ``size``."""
    from PIL import Image
    a = np.asarray(Image.open(png_a).convert("RGB").resize(size), dtype=float)
    b = np.asarray(Image.open(png_b).convert("RGB").resize(size), dtype=float)
    return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])
