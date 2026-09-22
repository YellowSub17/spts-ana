# Particle Sizing Analysis — Summary (Part 3)

## How small can the detection floor allow us to measure?

Extrapolating the physically-expected `intensity ∝ D⁶` relationship down from
the two least-censored PS groups (ps40nm, ps50nm — see part 2, section 7) to
the point where a *typical* (median) particle's intensity would equal the
detection floor (~1,150 raw intensity units):

- Using ps40nm as anchor: **~25 nm**
- Using ps50nm as anchor: **~29 nm**

Take these as order-of-magnitude, not precise. The exponent implied by those
two "cleanest" points directly (`ln(median₅₀/median₄₀) / ln(50/40)`) comes out
to **2.7**, not the 6 pure Rayleigh scattering predicts. That's a red flag:
even ps40nm/ps50nm's small sample sizes after all our filtering (n=278, n=78)
mean their medians themselves are noisy, so this extrapolation carries real
uncertainty rather than being a clean physical calculation.

What's confident is the qualitative picture, not the exact number:
**somewhere in the ballpark of 25–30 nm is where the median particle's signal
starts dropping toward the detection floor**, and below that, an increasing
majority of a population would be systematically invisible — consistent with
what was already seen directly (ps20nm already has ~31% of its distribution's
scale eaten by the floor, section 7). This also matches the DMA cross-check
(section 6): even the DMA had a comparable low-size cutoff (~5.8nm) where it,
too, couldn't resolve real signal from instrument floor — different
instrument, same fundamental issue at small sizes.

## How to lower the floor

Roughly in order of impact/feasibility:

1. **Check if the current threshold is set conservatively.** The observed
   floor (~1,000–1,250) sits at 3–4x the statistically-estimated camera noise
   (3σ ≈ 300–310). If the hit-detection threshold in the acquisition pipeline
   was chosen as a round/conservative number rather than tuned to the actual
   noise statistics, there may be free headroom to lower it toward ~2σ or so
   (accepting a slightly higher false-positive rate) without new hardware.
   This is the cheapest thing to check first — worth looking at whatever
   threshold parameter the `spts`/hit-finding code uses.
2. **Shorter illumination wavelength.** Rayleigh scattering cross-section
   scales as λ⁻⁴, so a shorter-wavelength laser gives a large signal boost
   specifically for small particles — probably the single biggest lever if
   switching wavelength is feasible.
3. **More laser power / longer exposure per frame** — directly raises photon
   counts per particle, raising SNR (as long as particle transit speed
   through the beam doesn't smear the signal too much at longer exposures).
4. **Lower-noise camera / better collection optics** (higher quantum
   efficiency, lower read noise, higher NA) — pushes the noise floor itself
   down, which independently lets the detection threshold be lowered.
5. **Reduce background/stray light in the chamber** (baffling, cleaner
   windows, better beam dump) — same effect as #4, attacking the noise side
   rather than the signal side.
6. **Frame-to-frame/trajectory averaging**, if particles are tracked across
   multiple frames rather than single-frame detection — integrating signal
   over a track effectively boosts SNR beyond what a single frame can
   achieve, letting particles that wouldn't clear the single-frame threshold
   be recovered.

## Following up on lever #1: what the acquisition pipeline actually does

Went through `src/spts` (the package used to detect/save the particle data)
and the `spts.conf` used for this analysis to check directly.

**Pipeline stages** (`worker.py`): `1_raw → 2_process → 3_denoise → 4_threshold
→ 5_detect → 6_analyse`. The relevant config for `data01271.cxi`:

```
[process]  floor_cut_level = 0, subtract_constant = 0
[denoise]  method = gauss, sigma = 0.1
[threshold] threshold = 20.0
[analyse]  integration_mode = windows, window_size = 25, circle_window = True
```

`window_size=25` with `circle_window=True` confirms the `is`/`peak_sum` field
is a sum over a circle of radius ~12px (not 25px as earlier scripts assumed —
worth double-checking `minimum_particle_size.py`'s `SUMMATION_RADIUS_PX`
against this). Correcting the noise-floor calculation for this smaller area
actually lowers the estimated statistical noise floor further (3σ ≈ 145,
down from the earlier 3σ ≈ 300–310 that assumed a bigger window) — reinforcing
that pure sensor/shot noise on the window sum is not what's setting the
observed ~1,000–1,250 floor.

**Why the floor is size-independent — mechanism identified.** `threshold=20`
is applied to the *denoised* image as a hard per-pixel cutoff (`threshold.py`:
`image_denoised >= threshold`) — a peak-amplitude test, not a total-flux test.
The `gauss` denoiser with `sigma=0.1` (a frequency-domain Gaussian) converts to
a real-space smoothing kernel of σ≈1.59px — almost identical to the particle
PSF width measured directly from the data (σ≈1.4–1.5px, from
`compute_focus_shape_metrics.py`'s Gaussian fits). Since every PS bead
(20–50nm) is far below the diffraction limit, they all produce the *same-width*
spot on the camera; only the amplitude scales with true brightness. So the
conversion factor between "peak pixel value" and "total window-summed
intensity" (`is`) is roughly constant across particle sizes — a fixed
peak-pixel threshold therefore produces an approximately fixed minimum total
intensity, regardless of particle size. That's the observed floor.

This was checked quantitatively: convolving the ~1.45px PSF with the ~1.59px
denoise kernel reduces peak amplitude by a factor of
`σ_psf²/(σ_psf²+σ_denoise²) ≈ 0.45`, so a denoised threshold of 20 requires a
raw peak pixel of ~44 counts, integrating (for a Gaussian PSF) to a total
window sum of ~44 × 2π×1.45² ≈ 580 — the right order of magnitude versus the
observed ~1,000–1,250 (off by ~2x, reasonable given the real PSF isn't a
perfect Gaussian).

**Headroom on the threshold.** The noise level on the *denoised* image itself
(not raw pixel noise) works out to only ~0.4 counts (a normalized Gaussian
smoothing kernel of σ=1.59px reduces IID pixel noise by `1/(2σ√π)`). The
current `threshold=20` sits at roughly **49σ** above that — extremely
conservative. A typical detection threshold targeting 5–10σ would still give
a very low false-positive rate over a multi-hundred-thousand-pixel,
thousands-of-frames dataset, and would leave a plausible 5–10x reduction in
`threshold` on the table. Since `intensity ∝ D⁶`, `D_min ∝ threshold^(1/6)`,
so a 5x lower threshold would only shrink the minimum detectable size by
~5^(1/6) ≈ 24% — real, but the sixth-power law limits how much any single
lever can buy.

## Background subtraction and flatfield correction — confirmed status

`src/spts/spts/scripts/cxd_to_h5.py` (the raw Hamamatsu CXD → `.cxi` converter,
upstream of the `spts.conf` pipeline above) has built-in background and
flatfield correction:

- `estimate_background()`: averages a dedicated no-injection background
  acquisition (`-b`/`--background-filename`) to get a per-pixel background,
  with 6σ bad-pixel masking. This becomes `entry_1/image_1/bg`.
- `estimate_flatfield()`: a separate "laser on paper" reference
  (`-f`/`--flatfield-filename`) for correcting *spatial illumination
  non-uniformity* across the frame. This becomes `entry_1/image_1/ff`.

Checked directly against `data01271.cxi`:

```
h5ls data01271.cxi/entry_1/image_1
bg                       Dataset {1780, 2048}
bg_fullframe             Dataset {2048, 2048}
```

**`bg`/`bg_fullframe` are present — background subtraction was applied.**
**No `ff`/`ff_fullframe` — flatfield correction was not applied.**

This cleanly separates the two issues this investigation has been chasing:

1. **Detection floor / calibration intercept** (sections above): already
   background-corrected, so the `threshold=20` mechanism above is the actual,
   confirmed explanation — not a missing background-subtraction step.
2. **Y-illumination gradient** (part 2, section on the intensity-vs-position
   check): flatfield correction was never applied to this dataset, which
   directly explains why the gradient is real rather than an analysis
   artifact. Reprocessing with `-f <flatfield_file>` (if a flatfield
   acquisition exists, or one can be taken) would correct the illumination
   profile at the source — a cleaner fix than the post-hoc y-position crop
   used in this analysis, which discards the affected region rather than
   correcting it, and could recover particles currently lost to the crop.
