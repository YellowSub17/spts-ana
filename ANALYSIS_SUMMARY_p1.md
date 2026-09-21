# Particle Sizing Analysis — Summary (Part 1)

## 1. Starting point

`thumbnails.h5` contains ~19,000 particle detections across six groups — four
polystyrene (PS) size standards (20/30/40/50 nm) and two proteins (ferritin,
GroEL). Each detection already carried a boolean `flags` field from a prior
focus-filtering process (the "old focus test"), along with shape descriptors
(area, eccentricity, etc.) and intensity statistics per particle.

## 2. Adding shape-based filtering (solidity)

The old focus test wasn't catching everything — a meaningful fraction of
particles it called "focused" were actually bean/kidney-shaped or
double-lobed, not clean single blobs. To catch these, for every thumbnail we:
- background-subtract and threshold the image to isolate the bright blob
  around the peak,
- take its connected-component area, and
- compute the area of its convex hull.

**Solidity** = blob area / convex hull area. A single round or slightly
elongated blob has solidity close to 1; a bean or double-lobed shape
(non-convex) has noticeably lower solidity. Across every group, this cleanly
separated the two: particles from the old test's "focused" set split into a
majority with solidity ≥ 0.80 (clean single blobs) and a substantial minority
(10–20%) below that, which visually confirmed as bean/multi-lobe shapes. We
adopted `flags AND solidity ≥ 0.80` as the combined filter, and now track
three categories throughout: RED (fails the old test), YELLOW (passes old
test, fails solidity), GREEN (passes both — the actual analysis set).

## 3. Intensity → sixth-root → size

For a Rayleigh scatterer, scattered intensity scales with volume², so
intensity^(1/6) is expected to scale linearly with particle diameter. Using
the GREEN (fully filtered) PS particles, we computed the median
intensity^(1/6) for each PS size and plotted it against nominal diameter.
ps20/30/40nm fell on an excellent line (R² ≈ 0.99); ps50nm did not (median
intensity was *lower* than ps40nm's — physically backwards), consistent with
that stock being old/degraded, so it was excluded from the calibration fit.
That fit was then used to read off a "PS-equivalent size" for ferritin
(~17 nm) and GroEL (~9 nm) — with the caveat that this equivalence ignores
refractive-index differences between polystyrene and protein, so it's not a
true physical size for the proteins.

## 4. Why the minimum measurable size can't be pinned down from this

Checking the calibration's reliability by applying it to individual particles
(rather than the group median) showed the problem: each PS group's
per-particle inferred size has a standard deviation of ~15–19 nm — as large
as, or larger than, the 10 nm spacing between the standards themselves. We
confirmed this isn't camera/background noise (signal-to-noise is 15–43× even
for the dimmest, smallest standard), so it's real particle-to-particle
variability — polydispersity, residual focus/z-depth differences, orientation
effects — that the intensity measurement can't distinguish from genuine size
differences.

The practical consequence: a population's *median* size can still be
estimated precisely (thanks to large particle counts), but any individual
particle's intensity-derived size is essentially meaningless, and that
per-particle scatter is wide enough to swamp the size differences between
neighboring standards. That means we don't currently have a reliable way to
say "we can measure particles down to X nm" — the intensity variance among
the standards is too large relative to the size steps being calibrated, so a
defensible minimum-detectable-size number can't be extracted from this data
as it stands.

## 5. Background stability check

To rule out the camera/illumination background as a contributor to that large
particle-intensity variance, we checked `data/bkgrnd.npy` — 100 particle-free
frames taken with the same setup.

- **Temporal stability**: frame-to-frame mean is 101.26 ± 0.11 counts
  (~0.1% relative), with no significant drift across the 100-frame
  acquisition. There is a small periodic wobble (~4-frame period, amplitude
  ~0.3 counts) — likely a minor readout/illumination flicker — but it's far
  too small to explain the ~15–17% intensity variance seen in the PS
  standards.
- **Spatial stability**: the per-pixel mean map shows only a mild
  illumination gradient (~1% of the mean, no sharp features), and the
  per-pixel temporal noise map agrees well with the background noise level
  independently estimated from thumbnail corners in the sizing-precision
  analysis. Only 26 out of 640,000 pixels (0.004%) show unusually high noise,
  scattered randomly rather than clustered — a few hot pixels, not a
  systematic defect.

**Conclusion**: the background is stable at the ~0.1–1% level, both spatially
and temporally. This rules out background instability as the source of the
~15–17% particle-intensity variance — that variance is coming from the
particles/optics themselves (polydispersity, focus residuals, etc.), not the
camera or illumination.

## 6. Tightening the solidity threshold (0.80 → 0.90) made things worse

To check whether the wide intensity spread (section 4) was due to residual
non-convex/bean-shaped contamination slipping through at the 0.80 solidity
cutoff, we tried raising the threshold to 0.90 and re-ran the calibration and
sizing scripts.

Result: it didn't clean up the data — it just discarded good particles and
made the statistics worse.

| | @0.80 (adopted) | @0.90 |
|---|---|---|
| ps20/30/40 sample sizes | 566 / 1839 / 440 | 153 / 553 / 140 |
| calibration R² | 0.989 | 0.944 |
| ferritin PS-equivalent size | 17.1 nm | 16.9 nm |
| GroEL PS-equivalent size | 9.4 nm | 6.1 nm |

Raising the threshold removed 65–70% of the already-filtered particles per PS
group, made the calibration fit *worse* (R² dropped from 0.989 to 0.944), and
shifted GroEL's estimated size by several nm — purely as a side effect of
losing statistical power, not from removing bad data. This indicates 0.80 is
a better operating point than 0.90, not just an arbitrary looser choice: the
particles between 0.80–0.90 solidity are evidently still good single-blob
detections, not contamination, and removing them only adds noise. (This also
means the wide intensity variance from section 4 is not explained by leftover
shape contamination — consistent with section 5's finding that it isn't the
background either.)
