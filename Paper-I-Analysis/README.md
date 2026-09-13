# Paper-I-Analysis

End-to-end, self-contained analysis for

> Singh, Sharma, Arora, Burman & Gandhi — *Linear Polarimetry of the Solar Corona
> during the Total Solar Eclipse of 8 April 2024: Multi-Band Observations and HDR
> Imaging* (Paper I).

Starting from the APP-calibrated light frames in `../Calibrated_Lights/`, the
numbered notebooks regenerate every table and figure in the paper. Nothing in
this folder imports from any other folder on the SSD; the only external inputs
are the calibrated frames (read-only) and, for the validation notebook, the
legacy products in `../FINAL_PAPER_ANALYSIS/`.

```
Paper-I-Analysis/
├── config.py            paths + every constant (168°, −26.24°, exposure bracket, LDIC thresholds, …)
├── utils.py             all helpers (CHT, Skyfield, stacking, HDR, Stokes, MGN, plotting), one copy each
├── run_all.sh           executes the notebooks in order with the project environment
├── requirements.txt
├── notebooks/
│   ├── 00_constants.ipynb               prints config, recomputes the P-angle (SunPy) and plate scale
│   ├── 01_frame_inventory.ipynb         header census of Calibrated_Lights → products/frame_inventory.csv
│   ├── 02_moon_center_cht.ipynb         circular Hough transform lunar limb → products/moon_centers.csv   (~75 min)
│   ├── 03_sun_center_skyfield.ipynb     Sun centre from Moon centre + DE421 → products/sun_moon_centers.csv
│   ├── 04_select_align_stack.ipynb      1-σ radius filter, Sun-centric alignment, mean stack → products/stacked/
│   ├── 05_hdr_exposure_normalization.ipynb  → products/hdr/hdr_expnorm_pol{1,2,3}.fits
│   ├── 06_hdr_ldic.ipynb                    → products/hdr/hdr_ldic_pol{1,2,3}.fits
│   ├── 07_paper_figures.ipynb           Figs. 6 & 7 → figures/{expnorm,ldic}_hdr_solar_frame.png
│   ├── 08_histogram_figure.ipynb        Fig. 8 → figures/histogram_normalized_exposures.png
│   └── 99_validate_against_legacy.ipynb compares every product with the files that went into the paper
├── data/                the CSV tables that produced the paper (copied from the SSD) + de421.bsp
├── products/            regenerated tables, stacks, HDR FITS, executed notebooks   (git-ignored)
├── figures/             regenerated figures                                        (git-ignored)
└── ARCHIVE.md           where every legacy notebook on the SSD lives and what it was for
```

## Running

```bash
# one-time environment (kept outside the exFAT SSD; uv cannot install into it)
uv venv --python 3.12 ~/.venvs/paper-i-analysis
uv pip install --python ~/.venvs/paper-i-analysis/bin/python -r requirements.txt

./run_all.sh            # everything, in order  (02 is slow; see below)
./run_all.sh 05 06 07   # a subset
```

Executed copies with outputs land in `products/executed/`; the notebooks in
`notebooks/` are committed without outputs. To work interactively, launch
`~/.venvs/paper-i-analysis/bin/jupyter lab` from this folder.

`02_moon_center_cht` takes ~30 s per frame. `config.MOON_CENTERS_SOURCE`
(default `"legacy"`) makes step 03 use the CHT table that produced the paper
(`data/sheet_for_centers_and_radius_and_SB.csv`), which is what reproduces the
paper exactly; set it to `"recomputed"` to use the table from notebook 02
(current scikit-image differs from the 2024 run by ±1 px on a few frames).

## Pipeline and provenance

| step | paper section | what | legacy source (on the SSD) |
|---|---|---|---|
| calibration | 2 | dark, flat, debayer → RGB FITS | **Astro Pixel Processor 2.0.0-beta29** (GUI, not python); output = `Calibrated_Lights/` |
| 01 | 2 | frame census | `categorize_data.ipynb` |
| 02 | 3.2 (1) | two-pass CHT on the RGB sum | `automated_sheet_for_centers_radius_and_SB.py` |
| 03 | 3.2 (2–3) | Skyfield Sun–Moon geometry, plate scale from the lunar radius, image-North rotation (78° from +y = 168° from +x) | `aligning_with_gemini_code/find_sun_center.ipynb` cells 9, 11 |
| 04 | 2, 3.2 | 1-σ lunar-radius filter (130/136 frames); luminance 0.2126R+0.7152G+0.0722B; bilinear shift of the Sun centre to (4128, 2752); **mean** over frames | filter/table: `FINAL_PAPER_ANALYSIS/statistical_analysis_figures.ipynb`; stacking: run off-line by C. Gandhi (recipe confirmed by him and by reconstruction, see below) |
| 05 | 3.3.1 | exposure-normalised HDR | `FINAL_PAPER_ANALYSIS/make_exposure_norm_HDR_from_stacked_exposures_Chaitanya.ipynb` |
| 06 | 3.3.2 | LDIC HDR with master-luminance weights | `FINAL_PAPER_ANALYSIS/make_ldic_hdr_from_stacked_exposures_Chaitanya.ipynb` |
| 07 | 3.4 | rotation to solar frame (−141.76°), colour composite, MGN σ = 10/20/40/80, compass + colour wheel | `FINAL_PAPER_ANALYSIS/{exposure_norm,ldic}_hdr_paper_figure_Paras.ipynb` |
| 08 | 3.2 | histogram figure | `FINAL_PAPER_ANALYSIS/statistical_analysis_figures.ipynb` cell 11 |

Constants that came from outside python:

* **Celestial North = 168°** CCW from the image +x axis — measured by hand in
  Stellarium + Affinity Photo (`../finding_celestial_north_with_stellarium.afphoto`,
  paper Fig. 5).
* **P-angle = −26.24°** — SunPy; recomputed in `00_constants` (−26.240°).
* **Inter-polariser shifts** `(dy, dx)` = pol1 (3, 2), pol2 (0, 0), pol3 (−4, −7) —
  set by eye in the legacy HDR notebooks "to align the plumes and make them
  white"; applied to the stacked images before HDR construction.

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

### The stacking step (04)

No script for the alignment + stacking existed on the SSD. The recipe was
recovered by reconstruction against `Linear_Composite_100.fits` and then
confirmed by C. Gandhi ("linearly translated the images as per the Sun's centres,
then took the mean"):

| candidate | corr | median \|Δ\| |
|---|---|---|
| luminance, bilinear sub-pixel shift to (W/2, H/2), **mean** | **0.999983** | 2×10⁻⁵ |
| same, median | 0.99986 | 3×10⁻⁵ |
| RGB mean / sum / green only | ≤ 0.998 | — |

The remaining 10⁻⁵-level residual is float32 summation / interpolation detail.
**The manuscript says "median stacking"; the data and the author say mean** —
the text should be corrected.

**Frame membership.** Fitting each legacy stack onto its candidate frames
(non-negative least squares, notebook-04 validation cell shows the result)
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

### Validation (notebook 99)

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

## Validation numbers (run of 2026-09-13, `products/validation_report.csv`)

Regenerated product vs the file that went into the paper (Pearson correlation;
stacks compared on the central 2500×3500 px window, HDR on the full frame):

| product | corr |
|---|---|
| stacked exposures, 27 planes | 0.99967 – 0.999996 (min: 1/800 s pol2) |
| `hdr_expnorm_pol{1,2,3}` | 0.999994 / 0.999988 / 0.999986 |
| `hdr_ldic_pol{1,2,3}` | 0.999991 / 0.999994 / 0.999954 |
| `expnorm_hdr_solar_frame.png` (Fig. 6) | 0.998 |
| `ldic_hdr_solar_frame.png` (Fig. 7) | 0.999 |
| `histogram_normalized_exposures.png` (Fig. 8) | 0.79 — same curves; the frame-count annotations now show the true stack membership and the Exp-Norm panel uses this pipeline's product instead of `HDR_Exposure_Normalization_Shivam` |
| Sun-centre table | max Δ 0.001 px (CSV rounding) |
| paper metadata table | identical up to the order of equal-exposure rows |

Notebook run times on an M-series Mac with the SSD: 01 10 s · 02 ~75 min ·
03 40 s · 04 5 min · 05 1 min · 06 5 min · 07 2 min · 08 40 s · 99 30 s.
