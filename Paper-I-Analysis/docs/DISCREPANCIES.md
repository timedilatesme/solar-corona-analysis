# Discrepancies found while consolidating the Paper I analysis

Everything below was found by reproducing the submitted figures from the
calibrated frames and comparing with the legacy notebooks/products on the SSD
(September 2026). Each item gives the evidence, the impact, and what this
repository does about it. Items marked **manuscript** need a change in the
text; items marked **decision** need a co-author decision.

## 1. Stacking is a mean, not a median — manuscript

*Evidence.* Reconstructing the paper's stacked composites
(`Linear_Composites_Instrument_Frame_Chaitanya/`) from the calibrated frames:
mean of the Sun-centred luminance frames gives corr 0.99998 with the paper's
files, median gives 0.99986; C. Gandhi confirmed "linearly translated … then
took the mean". Section 3.2 of the manuscript says "median stacking".
*Impact.* Wording only.
*Here.* `utils.stack_frames` averages; the notebook says so.

## 2. The stacks contain 124 frames, not 130 — manuscript / decision

*Evidence.* Six frames that pass the 1-σ lunar-radius filter are nevertheless
absent from the paper's stacks (non-negative least-squares fit of each legacy
stack onto its candidate frames gives them zero weight): `DSC_3810` (1/50 s
pol1 — first frame of totality, 3.4× brighter than its siblings, diamond
ring), `DSC_3868` (1/3 s pol2), `DSC_3935` (1/100 s pol3), `DSC_3952`
(1/200 s pol3), `DSC_3869` and `DSC_3887` (1/800 s pol2); all but the first have
a lunar radius 3–10 px off the rest of their group.
*Impact.* "130 useful frames" (Section 2) and the per-column counts in the
histogram figure describe the filter, not the stacks. Small effect on the
images (5 → 4 frames in a few stacks).
*Here.* `config.EXTRA_EXCLUDED_FRAMES` (with reasons) reproduces the paper's
stacks; the regenerated histogram annotates the true counts. Set it to `{}` to
stack all 130 and quote 130 — either way text and figure must agree.

## 3. Histogram figure mixed two data chains — manuscript (silently fixed)

*Evidence.* `statistical_analysis_figures.ipynb` drew the "Exp. Norm HDR" panel
from `HDR_Exposure_Normalization_Shivam/` (a different, rect-cropped alignment)
while the nine exposure panels and the LDIC panel came from the Chaitanya chain.
*Impact.* Cosmetic; the curves are similar.
*Here.* Both HDR panels come from this pipeline's own products.

## 4. Sun-centre rotation convention, and the by-eye channel shifts — decision

*Evidence.* The legacy `find_sun_center` notebook rotated the Sun–Moon
displacement into the image as `dx = sep·sin(PA−78°), dy = sep·cos(PA−78°)`.
Tested against the background stars (notebook 03: ζ Psc, 88 Psc, HD 7346 with
Simbad positions), this convention has the wrong handedness — it is the mirror
image of the one that fits the star field
(`n = (cos α, −sin α)`, `e = (−sin α, −cos α)` in array coordinates). Because
the Sun–Moon separation was only 27–40″ (18–27 px), the resulting Sun-centre
error is small, 3–11 px, **but it rotates with the Sun–Moon position angle and
therefore differs between the three polariser positions taken a minute apart:
pol1 (−4.0, −3.7) px and pol3 (+4.6, +4.3) px relative to pol2.** These are, to
within a couple of pixels, the "by-eye" shifts `(3,2)/(0,0)/(−4,−7)` that the
legacy HDR notebooks applied "to align the plumes and make them white".
*Impact.* The submitted figures are fine (the shifts compensated the error);
the *description* of the alignment in Section 3.2 does not mention that a
manual inter-channel shift was needed, and no such shift should be needed.
*Here.* Default `SUN_CENTER_ROTATION = "astrometric"` (star-validated) with
zero channel shifts; `"legacy"` reproduces the submitted products exactly.
Notebook 10 measures the residual inter-channel registration in both modes.

## 5. The CHT lunar radius is ~0.5 % too small — manuscript (plate scale)

*Evidence.* With the North angle and plate scale fitted to the stars and the
Moon held at its CHT centre, the stars sit only 1.7–2.2 px rms from their
predicted positions (notebook 03, `fit_rms_px`) — so the CHT Moon *centre* is
good to ~2 px. But the plate scale that fits the stars is 1.5005″/px whereas
the CHT lunar radius (668–676 px for 1012.2″) implies 1.506–1.52″/px: the CHT
circle is systematically 3–4 px (0.5 %) smaller than the lunar disk, as one
expects from an edge detector locking onto the inner side of the limb
gradient. At the per-frame lunar-radius scale the stars therefore appear
displaced by 10–20 px eastward (they lie 3000–4200″ east of the Moon), which is
what a naive check shows.
*Impact.* Negligible for the Sun centre (27–40″ × 0.5 % = 0.15″); 0.5 % on any
arcsec or R☉ scale derived from the lunar radius; the per-frame `sun_radius`
in the tables is 0.5 % too large in pixels.
*Here.* Notebook 10 uses the star-derived scale (1.500″/px) for its R☉
annuli and quotes the bias; `config.PLATE_SCALE_ARCSEC_PER_PIX` keeps the
paper's optical value for reference.

## 6. Plate scale: 1.49″/px (optical) vs 1.500″/px (measured) — manuscript

*Evidence.* The ζ Psc–88 Psc separation gives 1.5005 ± 0.001″/px (translation-
free; notebook 03). The paper quotes the nominal optical value 1.49″/px
(600 mm focal length, 4.35 µm pixels); the measured value corresponds to an
effective focal length of 598 mm.
*Impact.* 0.7 %; R☉ = 639 px rather than 643 px; any arcsec scale bar.
*Here.* `config.PLATE_SCALE_ARCSEC_PER_PIX` keeps the paper value for
reference; the alignment uses the per-frame value from the lunar radius (see
item 5) and the error budget uses the star value.

## 7. "Three background reference stars" — manuscript

*Evidence.* Two stars are identifiable in the frames at the required SNR:
ζ Psc (an unresolved 23″ double, A+B) and 88 Psc; HD 7346 (V = 8.5) is
marginal in the 1/3 s stack. The Stellarium overlay (Fig. 5) may have counted
ζ Psc A and B separately.
*Impact.* Wording.
*Here.* Notebook 03 uses the catalogue and reports what is detected.

## 8. CHT reproducibility across library versions — none

*Evidence.* Re-running the circular Hough transform with scikit-image 0.26
gives 1–3 px differences on ~20 frames and gross failures (> 10 px, the coarse
pass locking onto a wrong circle) on 6 of 145 frames relative to the 2024 run
(notebook 02 comparison cell; notebook 10: typical rms 0.9 px). The gross
failures are the frames with an inconsistent lunar radius, i.e. the ones the
1-σ filter (and `EXTRA_EXCLUDED_FRAMES`) already removes — the CHT is fragile
exactly where the frames are bad.
*Here.* `config.MOON_CENTERS_SOURCE = "legacy"` keeps the paper's table;
`"recomputed"` uses the new one.

## 9. Small legacy-data quirks — none

* `FRAME EXPOSURES.csv` on the SSD is an RTF document, not a CSV.
* `Calibrated_Lights/DSC_3770-Position4-cal-luminance.fits` is a stray test
  file that appears in the legacy tables; ignored here.
* Camera timestamps in the FITS headers are IST (`Asia/Kolkata`), converted to
  UTC before any ephemeris call (documented in `config.CAMERA_CLOCK_TZ`).
* The legacy CHT sheet has 311 rows (all frames of the day, including partial
  phases); only the 145 totality frames are on the SSD.
* 1/800 s pol2: the legacy stack is reproduced at corr 0.9997 rather than
  ≥ 0.9999 like every other plane — its membership (3 of 5 frames) was
  recovered, but two of the excluded frames could not be assigned a reason
  beyond "lunar radius 680 px".

## 10. P-angle and celestial North — confirmed, no discrepancy

P = −26.240° (SunPy, notebook 00) and α = 168.0 ± 0.2° from the star fits
(notebook 03) agree with the paper.
