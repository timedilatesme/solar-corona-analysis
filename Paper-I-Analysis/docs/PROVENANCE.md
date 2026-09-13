# Provenance of the Paper I figures and the legacy code

Written when the scattered notebooks on the SSD were consolidated into this
folder (September 2026). It records **which legacy files produced which paper
figure**, how the missing stacking step was recovered, and the validation
numbers of the first reproduction. For the current layout and how to run, see
[../README.md](../README.md); for the things that turned out to be wrong or
inconsistent, see [DISCREPANCIES.md](DISCREPANCIES.md).

Notebook numbers below use the current numbering (00 constants · 01 inventory ·
02 CHT · 03 celestial North from stars · 04 Sun centre · 05 select/align/stack ·
06 exp-norm HDR · 07 LDIC HDR · 08 paper figures · 09 histogram · 10 error budget · 99 validation).

## Pipeline and provenance

| step | paper section | what | legacy source (on the SSD) |
|---|---|---|---|
| calibration | 2 | dark, flat, debayer → RGB FITS | **Astro Pixel Processor 2.0.0-beta29** (GUI, not python); output = `Calibrated_Lights/` |
| 01 | 2 | frame census | `categorize_data.ipynb` |
| 02 | 3.2 (1) | two-pass CHT on the RGB sum | `automated_sheet_for_centers_radius_and_SB.py` |
| 04 | 3.2 (2–3) | Skyfield Sun–Moon geometry, plate scale from the lunar radius, image-North rotation (78° from +y = 168° from +x) | `aligning_with_gemini_code/find_sun_center.ipynb` cells 9, 11 |
| 05 | 2, 3.2 | 1-σ lunar-radius filter (130/136 frames); luminance 0.2126R+0.7152G+0.0722B; bilinear shift of the Sun centre to (4128, 2752); **mean** over frames | filter/table: `FINAL_PAPER_ANALYSIS/statistical_analysis_figures.ipynb`; stacking: run off-line by C. Gandhi (recipe confirmed by him and by reconstruction, see below) |
| 06 | 3.3.1 | exposure-normalised HDR | `FINAL_PAPER_ANALYSIS/make_exposure_norm_HDR_from_stacked_exposures_Chaitanya.ipynb` |
| 07 | 3.3.2 | LDIC HDR with master-luminance weights | `FINAL_PAPER_ANALYSIS/make_ldic_hdr_from_stacked_exposures_Chaitanya.ipynb` |
| 08 | 3.4 | rotation to solar frame (−141.76°), colour composite, MGN σ = 10/20/40/80, compass + colour wheel | `FINAL_PAPER_ANALYSIS/{exposure_norm,ldic}_hdr_paper_figure_Paras.ipynb` |
| 09 | 3.2 | histogram figure | `FINAL_PAPER_ANALYSIS/statistical_analysis_figures.ipynb` cell 11 |

Constants that came from outside python:

* **Celestial North = 168°** CCW from the image +x axis — measured by hand in
  Stellarium + Affinity Photo (`../finding_celestial_north_with_stellarium.afphoto`,
  paper Fig. 5).
* **P-angle = −26.24°** — SunPy; recomputed in `00_constants` (−26.240°). Celestial North is re-derived from the stars in notebook 03 (168.0 ± 0.1°).
* **Inter-polariser shifts** `(dy, dx)` = pol1 (3, 2), pol2 (0, 0), pol3 (−4, −7) —
  set by eye in the legacy HDR notebooks "to align the plumes and make them
  white"; applied to the stacked images before HDR construction. **Explained** in
  DISCREPANCIES.md §4: they compensate the legacy Sun-centre rotation convention
  and are not applied in the default (astrometric) mode.

### Which legacy data went into the paper

`FINAL_PAPER_ANALYSIS/` holds parallel `_Shivam` and `_Chaitanya` chains. The
paper figures come from the **Chaitanya** chain: the LaTeX references
`expnorm_hdr_solar_frame.png` / `ldic_hdr_solar_frame.png`, which are written by
the `*_Paras.ipynb` figure notebooks from `ExpNorm_HDR_Images_Paras/` and
`LDIC_HDR_Images_Paras/`, which are built from
`Linear_Composites_Instrument_Frame_Chaitanya/` (byte-identical to
`Sun_Linear_Composite/Linear_Composite_*.fits`, April 2025). Pixel correlation of
the PNGs on disk with the inline outputs embedded in each notebook: 0.998 with
the Paras notebooks, 0.72 with the Shivam one. The Shivam chain (rect-cropped
5356×4121 stacks) is not used here.

### The stacking step (notebook 05)

Linearly translated the images as per the Sun's centres, then took the mean:

| candidate | corr | median \|Δ\| |
|---|---|---|
| luminance, bilinear sub-pixel shift to (W/2, H/2), **mean** | **0.999983** | 2×10⁻⁵ |
| same, median | 0.99986 | 3×10⁻⁵ |
| RGB mean / sum / green only | ≤ 0.998 | — |

The remaining 10⁻⁵-level residual is float32 summation / interpolation detail.
**The manuscript says "median stacking"; the data and the author say mean** —
the text should be corrected.

**Frame membership.** Fitting each legacy stack onto its candidate frames
(non-negative least squares; the notebook-05 validation cell shows the per-plane result)
reveals that six frames which pass the 1-σ radius filter were nevertheless
*not* in the paper's stacks: `DSC_3810` (1/50 s pol1, the first frame of
totality and 3.4× brighter than its siblings — diamond ring), `DSC_3868`
(1/3 s pol2), `DSC_3935` (1/100 s pol3), `DSC_3952` (1/200 s pol3), `DSC_3869`
and `DSC_3887` (1/800 s pol2); all but the first have a lunar radius that
disagrees with the rest of their group by 3–10 px. They are listed with reasons
in `config.EXTRA_EXCLUDED_FRAMES`; with that list the regenerated stacks match
the legacy ones at corr ≥ 0.9997 on every plane. **So the paper's "130 useful
frames" and the per-column frame counts in the histogram figure describe the
radius filter, whereas the stacks actually contain 124 frames** — the
regenerated histogram figure annotates the true counts. Set
`EXTRA_EXCLUDED_FRAMES = {}` to stack all 130.

### Validation of the first reproduction (legacy geometry, notebook 99)

Every regenerated product is compared with its legacy counterpart; the report is
written to `products/validation_report.csv`. See the bottom of this README for
the numbers from the last full run.

## Notes for the co-authors

* Stacking is a **mean**, not a median (Section 3.2 of the manuscript).
* The stacks behind Figs. 6–8 contain **124** frames, not 130: six frames that
  pass the 1-σ lunar-radius filter were additionally left out (list and reasons
  in `config.EXTRA_EXCLUDED_FRAMES`). Either the text/figure counts should say
  124, or the stacks should be rebuilt with all 130 (`EXTRA_EXCLUDED_FRAMES = {}`).
* The histogram figure previously mixed `LDIC_HDR_Images_Paras` with
  `HDR_Exposure_Normalization_Shivam` for its two HDR panels; here both panels use
  this pipeline's own products, so the whole figure has one provenance.
* Plate scale from the lunar radius is 1.506″/px (mean CHT radius 672.3 px vs
  1012.2″); the paper quotes the optical value 1.49″/px.

## Validation numbers of the first reproduction (2026-09-13, legacy geometry)

Regenerated product vs the file that went into the paper (Pearson correlation;
stacks compared on the central 2500×3500 px window, HDR on the full frame):

| product | corr |
|---|---|
| stacked exposures, 27 planes | 0.99967 – 0.999996 (min: 1/800 s pol2) |
| `hdr_expnorm_pol{1,2,3}` | 0.999994 / 0.999988 / 0.999986 |
| `hdr_ldic_pol{1,2,3}` | 0.999991 / 0.999994 / 0.999954 |
| `expnorm_hdr_solar_frame.png` (Fig. 6) | 0.998 |
| `ldic_hdr_solar_frame.png` (Fig. 7) | 0.999 |
| `histogram_normalized_exposures_vertical.png` (Fig. 8) | 0.79 — same curves; the frame-count annotations now show the true stack membership and the Exp-Norm panel uses this pipeline's product instead of `HDR_Exposure_Normalization_Shivam` |
| Sun-centre table | max Δ 0.001 px (CSV rounding) |
| paper metadata table | identical up to the order of equal-exposure rows |

Notebook run times on an M-series Mac with the SSD: 01 10 s · 02 ~75 min · 03 ~17 min ·
04 40 s · 05 5 min · 06 1 min · 07 5 min · 08 2 min · 09 1 min · 10 ~10 min · 99 30 s.
