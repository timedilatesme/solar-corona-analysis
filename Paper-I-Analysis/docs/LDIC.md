# How the Linear Digital Image Composer (LDIC) is used in this work

This is the expanded description of Section 3.3.2 of the paper, written so that
it can be lifted into the manuscript or a referee reply. Code:
[`utils.ldic_hdr_stacking`](../utils.py) (driven by
[`notebooks/07_hdr_ldic.ipynb`](../notebooks/07_hdr_ldic.ipynb)); parameters in
[`config.LDIC_PARAMS`](../config.py).

## 1. The problem LDIC solves

A single exposure cannot hold the corona's dynamic range: the innermost corona
saturates in the 1/3 s frames while the outer corona is lost in noise in the
1/800 s frames. Simply averaging the exposure-normalised frames (Section 3.3.1)
weights every exposure equally, so saturated pixels of the long exposures and
noise-dominated pixels of the short exposures both leak into the result.
LDIC (Druckmüllerová 2014, Ch. 4) instead builds, pixel by pixel, a weighted
average in which **each pixel is taken from the exposures in which it is well
exposed**, after bringing all exposures onto a common linear intensity scale.
The composite remains a *linear* function of the incident irradiance — no
tone mapping is involved — so it can still be used quantitatively.

## 2. The formula

For $N$ stacked exposures $f_i(r,\phi)$ (index $i=0$ = longest exposure,
$i=N-1$ = shortest; polar coordinates about the Sun centre) the composite is

$$
g(r,\phi)=\frac{\sum_{i=0}^{N-1} w_i\!\big(f_i(r,\phi)\big)\,\big[k_i(\phi)\,f_i(r,\phi)+q_i(\phi)\big]}
              {\sum_{i=0}^{N-1} w_i\!\big(f_i(r,\phi)\big)} .
$$

Three ingredients: a **weight** $w_i$ that says how trustworthy a pixel of
exposure $i$ is, a **gain** $k_i(\phi)$ that puts exposure $i$ on the intensity
scale of the exposures already composed, and an **offset** $q_i(\phi)$ for
residual background differences.

## 3. Weights — which pixels are trusted

$w_i$ depends only on the pixel value, expressed as a fraction of the sensor's
dynamic range (the frames are scaled so that the brightest pixel of the whole
bracket is 65535):

| pixel value | weight |
|---|---|
| < 2 % | 0 — noise-dominated |
| 2 % → 15 % | ramps linearly 0 → 1 |
| 15 % → 65 % | 1 — the linear, well-exposed part of the response |
| 65 % → 80 % | ramps linearly 1 → 0 |
| > 80 % | 0 — approaching saturation / non-linear |

Two exceptions keep the ends of the bracket usable: for the **longest** exposure
the low-end rejection is switched off (otherwise the faintest outer corona,
which only that exposure records, would be discarded), and for the **shortest**
exposure the high-end rejection is switched off (otherwise the brightest inner
corona, which only that exposure records unsaturated, would vanish into a dark
"moat" around the Moon). The ramps make the hand-over between exposures smooth,
so no seams appear where one exposure takes over from the next.

![weight function and gains](../figures/ldic_weights_and_gains.png)

*(left: the three variants of the weight function; right: the fitted gains
$k_i(\phi)$, see §5.)*

## 4. The master-luminance rule — specific to polarimetry

LDIC was designed for a single image. We have three polariser channels
(0°, +60°, −60°) that must keep their intensity *ratios* exact, because the
ratios are the polarisation signal. If the weight were evaluated on each channel
separately, the hand-over between exposures would happen at slightly different
pixels in the three channels (their intensities differ by up to ~40 % at a given
pixel), and every hand-over seam would become a ring of false colour.

We therefore evaluate $w_i$ on a **master luminance** image — the pixel-wise
maximum of the three (scaled) channels for that exposure — and apply the *same*
weight map to all three channels. Using the maximum rather than the mean makes
the saturation test conservative: a pixel is rejected as soon as *any* channel
approaches saturation. This is the one modification to the published algorithm
and it is what removes the colour artefacts at the exposure boundaries.

## 5. Gains — putting exposures on one scale

Exposures are processed from the longest to the shortest. The longest is the
reference ($k_0=1$, $q_0=0$). For each following exposure $i$:

1. the running composite of exposures $0\ldots i-1$, $g^{(i-1)}_{\rm avg}$, is
   the target;
2. the image is cut into $n_s = 60$ angular sectors of 6°;
3. in each sector, pixels where both exposure $i$ is trusted ($w_i \ge 0.8$) and
   the running composite already has signal are used for a linear regression of
   $g^{(i-1)}_{\rm avg}$ against $f_i$ **through the origin**: $k_i = \sum f g/\sum f^2$,
   $q_i = 0$. The offset is fixed at zero because the stacked frames are already
   sky-subtracted; letting $q_i$ float over-fits the noise and produces negative
   trenches in the composite;
4. the 60 sector values $k_i(\phi_s)$ are smoothed by least-squares fitting a
   trigonometric polynomial $c_0+\sum_{p=1}^{4}(a_p\cos p\phi + b_p\sin p\phi)$
   (order reduced automatically if fewer than 9 sectors are valid), giving a
   continuous $k_i(\phi)$; a sector is skipped if it has fewer than two usable
   pixels, and the whole exposure is skipped if fewer than 15 sectors are valid.

Nominally $k_i \approx t_0/t_i$ (the exposure-time ratio); fitting it instead of
imposing it absorbs shutter-timing errors, the slight non-linearity of the
sensor near saturation and vignetting differences, and its variation with $\phi$
is a diagnostic (a few per cent here — see `docs/ERROR_BUDGET.md`).

## 6. Composition and normalisation

Every exposure's contribution $w_i\,[k_i f_i + q_i]$ is accumulated together with
its weight, and the composite is the weight-normalised sum. Pixels that no
exposure trusts (only the very centre of the lunar disk) are set to zero. The
composite is stored as float32 FITS in the intensity units of the longest
exposure; the exposure-normalised HDR of Section 3.3.1 is in counts per second,
so the two products differ by a constant factor.

## 7. Parameters used

| parameter | value | in `config.LDIC_PARAMS` |
|---|---|---|
| exposures | 1/800 … 1/3 s (9) | `INV_EXPOSURES` |
| angular sectors | 60 | `num_angular_segments` |
| trig. polynomial order (k, q) | 4, 4 | `trig_poly_order_k/q` |
| weight thresholds | 2 %, 15 %, 65 %, 80 % | `wf_*_percent` |
| regression pixel threshold | $w_i \ge 0.8$ | `regression_weight_threshold` |
| minimum valid sectors | 15 | `min_valid_segments` |
| intercept | fixed 0 | `fit_intercept = False` |
| dynamic-range scale | 65535 | `LDIC_TARGET_MAX` |
| Sun centre | (4128, 2752) | `SUN_CENTER_XY` |

## 8. Pseudocode

```
scale each channel so that max over the bracket = 65535
ref[i] = max(pol1[i], pol2[i], pol3[i])            # master luminance per exposure
for channel in (pol1, pol2, pol3):
    g = 0; W = 0
    for i, f in enumerate(exposures longest -> shortest):
        w = weight(ref[i]; bypass low end if i == 0, high end if i == N-1)
        if i == 0: k(phi) = 1
        else:
            for each of 60 sectors: k_s = sum(f*g/W) / sum(f*f) over trusted pixels
            k(phi) = trig-poly fit of k_s (order <= 4)
        g += w * k(phi) * f ;  W += w
    composite = g / W  (0 where W == 0)
```

## 9. What differs from Druckmüllerová (2014)

* weights evaluated on the master luminance shared by the three polariser
  channels (§4) — new, required for polarimetry;
* $q_i \equiv 0$ (sky already subtracted);
* bypass of the low-end rejection for the longest and of the high-end rejection
  for the shortest exposure (§3);
* adaptive reduction of the trigonometric-polynomial order when few sectors are
  valid;
* a single global scale factor per channel (rather than per-image
  normalisation), so the composite stays linear across the bracket.
