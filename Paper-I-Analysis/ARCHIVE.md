# Legacy code on the SSD — where it is and what it was for

Everything below lives under `/Volumes/WD_SSD/Solar Eclipse 2024 Data/` and is
left untouched. Paths are relative to that folder. "Migrated" = re-implemented
in this folder's `utils.py` / notebooks; nothing else is needed to reproduce the
paper.

## Migrated (source of this pipeline)

| legacy file | now |
|---|---|
| `categorize_data.ipynb` | `01_frame_inventory` |
| `automated_sheet_for_centers_radius_and_SB.py` (6-column CHT sheet; the one whose output is on disk) | `utils.find_moon_center_radius`, `02_moon_center_cht` |
| `aligning_with_gemini_code/find_sun_center.ipynb` cells 9 (batch) & 11 (mosaic) + `de421.bsp` | `utils.Ephemeris`, `utils.sun_center_from_moon`, `03_sun_center_skyfield` |
| `FINAL_PAPER_ANALYSIS/statistical_analysis_figures.ipynb` cells 3–7 (1-σ filter, metadata table) and 11 (histogram figure) | `utils.build_frame_table`, `04_select_align_stack`, `08_histogram_figure` |
| *(no script on the SSD — C. Gandhi's off-line alignment + mean stack)* | `utils.stack_frames`, `04_select_align_stack` |
| `FINAL_PAPER_ANALYSIS/make_exposure_norm_HDR_from_stacked_exposures_Chaitanya.ipynb` | `05_hdr_exposure_normalization` |
| `FINAL_PAPER_ANALYSIS/make_ldic_hdr_from_stacked_exposures_Chaitanya.ipynb` | `utils.ldic_hdr_stacking`, `06_hdr_ldic` |
| `FINAL_PAPER_ANALYSIS/exposure_norm_hdr_paper_figure_Paras.ipynb`, `ldic_hdr_paper_figure_Paras.ipynb` | `07_paper_figures` (+ `utils` §6–7) |
| `Sun_Linear_Composite/change_reference_frame_all_images.ipynb` cell 2 (Stokes helpers) | `utils.mzp_to_stokes`, `stokes_to_mzp` |
| `Paper-II-Analysis/tests_Claude/pipeline_utils.py` (helper source, copied not imported) | `utils` |
| `sun_moon_centers_output.csv`, `sheet_for_centers_and_radius_and_SB.csv`, `FINAL_PAPER_ANALYSIS/paper_metadata_documentation.csv` | `data/` |

Legacy products used only for validation (`99_validate_against_legacy`):
`FINAL_PAPER_ANALYSIS/Linear_Composites_Instrument_Frame_Chaitanya/`,
`ExpNorm_HDR_Images_Paras/`, `LDIC_HDR_Images_Paras/`, and the three paper PNGs
in `FINAL_PAPER_ANALYSIS/`.

## Superseded — Shivam-data variants (not used in the paper)

Parallel chain built on rect-cropped 5356×4121 stacks. Kept on the SSD only.

* `FINAL_PAPER_ANALYSIS/make_ldic_hdr_from_stacked_exposures_Shivam.ipynb`
* `FINAL_PAPER_ANALYSIS/old_make_ldic_hdr_from_stacked_exposures.ipynb`
* `FINAL_PAPER_ANALYSIS/exposure_norm_hdr_paper_figure_Shivam.ipynb` → `exposure_normalized_hdr_solar_frame.png`
* `FINAL_PAPER_ANALYSIS/ldic_hdr_paper_figure_Shivam.ipynb`
* `FINAL_PAPER_ANALYSIS/HDR_Exposure_Normalization_Shivam/`, `LDIC_HDR_Images_Shivam/`, `Linear_Composites_Instrument_Frame_Shivam/`
* `Paper-II-Analysis/Shivam_data_Calibrated_Aligned/` (81 per-band stacks; Paper II input)

## Experiments off the final chain

* `FINAL_PAPER_ANALYSIS/make_ldic_hdr_from_stacked_exposures_Chaitanya_mgn.ipynb` — LDIC on the `_mgn.fits` composites
* `FINAL_PAPER_ANALYSIS/make_ldic_hdr_from_stacked_exposures_Chaitanya_single_exposure.ipynb` — identical copy of the main LDIC notebook
* `automated_sheet_for_centers_and_radius_and_SB.py` — per-channel CHT variant (18-column sheet)
* `automated_NRGF.py` — NRGF filter batch
* `Sun_Linear_Composite/make_mgn.ipynb` → `Linear_Composite_*_mgn.fits`; `polarization_angle.ipynb`, `change_reference_frame_polarization_angle.ipynb`
* `Sun_Linear_Composite/Solar_Frame_Images/` — April–May 2025 HDR/MGN/figure iterations (`make_hdr_from_stacked_exposures`, `make_mgn_hdr_from_stacked_exposures`, `make_mgn`, `solar_corona_paper_figure`, `make_color_wheel`)
* `Sun_Linear_Composite/Christian_HDR_Algo/` — Julia LDIC port and its python re-implementation (historical reference for step 06)
* `Linear_Composite/` — March 2025 PNG/JPEG composites; `HDR/`, `Calibrated_Lights/HDR/` — early HDR trials

## Dead ends (June 2024 – March 2025)

* pystackreg / NRGF / FNRGF era: `test_pystackreg_wei_hao_data.ipynb`, `test_pystackreg_our_data.ipynb`, `test_NRGF_pystackreg_our_data*.ipynb`, `enhancing an image with FNRGF.ipynb` (root and `generate_combined_polarization_map/`), `test_processing.ipynb`, `test_find_moon_center.ipynb` (prototype of step 02), `NRGF_Filtered_Images/`, `Tangential Data/`, `original_and_tangential_pngs/`, `filtered_*.npy/png`
* phase-correlation / OpenCV registration era: `NEW_/` (`image_registration.ipynb`, `imregpoc.py`, `imreg/`, `opencv-feature-matching.ipynb`, `deconvolution.ipynb`, `mgn_test.ipynb`, `nrgf_*`), `align_two_images_that_share_features/`, `astro_stacker.py` (third-party), `Old_Aligned_Calibrated_Frames/`, `Calibrated_Lights/merge_pos123_aligned.ipynb`, `Position 4 Stacking Siril/`
* `test_flats.ipynb`, `flats/`, `polarization_angle_relation_with_triplets_MZP.ipynb` (stub), `polarized-solar-corona-paper-analysis/` (abandoned one-file consolidation attempt)

## Non-code assets

* `Calibrated_Lights/` — the true input (APP output, 145 frames + `BPM-*.fits`); `DSC_3770-Position4-cal-luminance.fits` is a stray test file
* `D800 TSE Images/`, `GoPro Timelapse TSE24/` — second camera / timelapse (progression figure)
* `finding_celestial_north_with_stellarium.afphoto`, `Aligned_Stacked_Calibrated_frames/stellarium_*` — the manual 168° celestial-North measurement
* `*.key`, PDFs, `FRAME EXPOSURES.csv` (actually an RTF document)
