# Particle Sizing Analysis — Summary (Part 2)

Current state of the pipeline and findings. See `scripts/` for the code
(`filter_config.py` holds the shared filter constants referenced below) and
`figures/` for the plots (both regenerable by re-running the scripts).

## 1. Data

`data/thumbnails.h5` holds ~19,000 particle detections across six groups:
polystyrene (PS) size standards at 20/30/40/50 nm, and two proteins (ferritin,
GroEL). Each detection has a 60×60 pixel thumbnail, an (x, y) frame position,
summed intensity, and several shape descriptors.

## 2. Filtering pipeline

A particle is included in the analysis if it passes all three of:

1. **`flags`** — the instrument's own focus flag (already present in
   `thumbnails.h5`).
2. **Solidity ≥ 0.80** — `solidity = blob area / convex hull area`, computed
   per particle from the thumbnail in `compute_focus_shape_metrics.py` and
   written to `data/focus_shape_metrics.h5`. This rejects non-convex shapes
   (bean/kidney/double-lobed blobs) that `flags` alone doesn't catch.
3. **y-position ∈ [30, 167] px** — a crop on the particle's y-coordinate in
   the frame, removing a region with a y-dependent illumination gradient
   (median intensity falls off steadily for y > 167, down to ~75% of its
   plateau value by y ≈ 230). The x-coordinate shows no comparable trend, so
   only y is cropped.

Particles passing all three are referred to as **GREEN** (see
`plot_focus_shape_metrics.py`, which also visualizes the particles that fail
each stage). Within the GREEN set, the top/bottom 2% of summed intensity is
additionally trimmed per group before computing any statistics, to reduce the
influence of aggregates/debris.

Background stability was checked separately (`data/bkgrnd.npy`,
`background_stability.py`): the camera/illumination background is stable to
~0.1–1% both temporally and spatially, and is not a source of the intensity
variance discussed below.

## 3. Sixth-root intensity calibration

For a Rayleigh scatterer, scattered intensity scales with volume², so
intensity^(1/6) is expected to scale linearly with particle diameter. Using
the median intensity^(1/6) of the GREEN, y-cropped particles for each PS
group:

```
intensity^(1/6) = 0.0521 * size_nm + 3.058        R^2 = 0.9997
```

fit across all four PS standards (20/30/40/50 nm) — see
`ps_calibration_sizing.py` and `plot_ps_calibration_sizing.py`.

## 4. Ferritin / GroEL sizing

Reading ferritin's and GroEL's median intensity^(1/6) off that calibration
line gives a **polystyrene-equivalent size**:

| particle | PS-equivalent size | 16–84th percentile range |
|---|---|---|
| ferritin | ~18.8 nm | 11.0–30.9 nm |
| GroEL | ~12.6 nm | 9.9–25.0 nm |

**Caveat:** this calibration is built entirely on polystyrene, which has a
uniform, high refractive index. Scattered intensity depends on
volume² × (optical contrast)², not volume alone, so a particle with different
optical contrast than PS will not map onto this line at its true physical
size. Ferritin carries a dense iron-oxide core (higher contrast than PS per
unit volume); GroEL is a low-density protein-only complex (lower contrast,
closer to buffer). The numbers above should be read as PS-equivalent
scattering sizes, not true physical diameters, until a refractive-index
correction is applied.

## 5. Minimum measurable particle size

Two distinct questions:

- **Is the signal strong enough to detect a particle at all?** Estimating the
  background/shot-noise floor on the summed-intensity measurement from the
  thumbnails' corner pixels (`minimum_particle_size.py`) gives an SNR of
  ~16–60× across the PS groups (smallest for ps20nm). Detection sensitivity
  is not the limiting factor, even for the smallest standard.

- **Given a detected particle, how precisely can its size be determined?**
  Applying the calibration to every individual GREEN particle (rather than
  the group median) gives a per-particle inferred-size standard deviation of
  ~11–14 nm across the PS groups — comparable to the 10 nm spacing between
  standards. A single particle's intensity-derived size is therefore not
  reliable on its own. The population *median* is precise (SEM well under
  1 nm, thanks to large particle counts), but only for a large, clean,
  single-species population — it says nothing about any individual particle,
  and does not help if the population itself is mixed or biased.

**Bottom line:** minimum measurable size is meaningful only as a statement
about a population's median, not about individual particles, and that
statement is only as good as the assumption that the population is a single
clean species imaged under uniform conditions.

## 6. Cross-check against DMA (independent sizing measurement)

`data/AIM148.txt` holds 16 DMA (differential mobility analyzer) scans from the
same day, with a confirmed run log mapping scan number to sample: scans 1–3 =
50nm PS, 4–7 = 40nm PS, 8–9 = 30nm PS, 10–11 = 20nm PS, 12–14 = GroEL 1uM,
15–16 = Ferritin (`plot_dma_data.py`). Within most blocks, the first scan(s)
after switching samples still show carryover from the previous sample (peak
pinned at the DMA's lower limit, ~5.8nm); the later scan(s) are the trustworthy
steady-state measurement. 30nm PS and Ferritin only got 2 scans each and
neither settles into a clean distribution, so their DMA statistics are
unreliable — no meaningful comparison is drawn for those two.

Converting both measurements to the same metric (GSD, geometric std. dev. =
D84/D50) gives a direct, apples-to-apples comparison of measurement spread:

| sample | optical GSD | DMA GSD | ratio (optical/DMA) |
|---|---|---|---|
| 20nm PS | 1.79 | 1.31 | 1.37 |
| 40nm PS | 1.41 | 1.56 | 0.90 |
| 50nm PS | 1.38 | 1.93 | 0.72 |
| GroEL | 1.59 | 1.42 | 1.12 |
| 30nm PS | 1.60 | 2.00 | — (DMA unreliable) |
| Ferritin | 1.68 | 1.50 | — (DMA unreliable) |

For every sample with usable DMA data, the optical and DMA spreads are within
~10–40% of each other — the same order of magnitude, not wildly different.
For 40nm and 50nm PS, the optical measurement is actually *tighter* than the
independently-measured physical polydispersity.

**This reframes the "error is too large" concern.** Since scattered intensity
scales as diameter⁶, even the DMA-confirmed real physical polydispersity
(GSD ~1.3–1.9) is expected to produce a large spread in optical intensity/
inferred size on its own — most of what looked like measurement noise is
consistent with genuine sample polydispersity, amplified by the sixth-power
relationship, not a flaw specific to the optical pipeline. The one exception is
20nm PS, where the optical spread is measurably wider than the DMA spread —
that sample still shows some real excess optical noise beyond what its
physical polydispersity explains, and would be the place to focus further
investigation.
