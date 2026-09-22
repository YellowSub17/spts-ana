"""
Checks whether the PS-calibration intensity variance is small enough to extract a
meaningful particle size, and what that implies for the minimum measurable size.

Two separate questions get conflated if you're not careful, so this script keeps
them apart:

  (a) Is the SIGNAL big enough to detect a particle at all? -- a sensor/shot-noise
      question. Answered here by estimating the background noise floor on the
      summed-intensity ('is') measurement directly from the thumbnails' corner
      pixels (propagated over the ~25px summation window used upstream), and
      comparing it to the median signal of each PS group.

  (b) Given a detected particle, how precisely can we assign it a SIZE? -- an
      entirely different question, answered by running the size calibration
      (see ps_calibration_sizing.py) on every individual GREEN particle (rather
      than just the population median) and looking at the spread of inferred
      sizes.

Finding (see ANALYSIS_SUMMARY_p2.md / run this script): (a) is not the
bottleneck -- even the dimmest standard (ps20nm) sits well above the estimated
noise floor. But (b) is: the per-particle inferred-size standard deviation
(~11-14 nm) is comparable to the 10nm spacing between the PS standards. So a
size read off a SINGLE particle's intensity is not meaningful on its own -- the
population median is still well-determined (SEM well under 1nm, thanks to large
N), but only for a large, clean, single-species population. This bounds "minimum
measurable size" as a statement about populations, not individual particles: you
can characterize the median size of a big, clean population of small particles,
but you cannot reliably size any one of them individually, and you cannot rescue
a systematically-biased or mixed population (e.g. ferritin/GroEL's
refractive-index mismatch -- see ps_calibration_sizing.py) by adding more
particles.

Run compute_focus_shape_metrics.py first to generate data/focus_shape_metrics.h5.
"""

import os

import h5py
import numpy as np
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_focus_shape_metrics as pfsm
from filter_config import Y_ILLUMINATION_MIN, Y_ILLUMINATION_MAX

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"

TRIM_PERCENTILES = (2, 98)
PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm", "ps50nm"]
PS_NOMINAL_SIZES_NM = {"ps20nm": 20, "ps30nm": 30, "ps40nm": 40, "ps50nm": 50}

# detection/summation window radius (px) used upstream for the 'is' intensity field
# (see src/comborun/comborun.py: ComboRun(..., 25) -> the r1 / large window)
SUMMATION_RADIUS_PX = 25
N_BACKGROUND_SAMPLES_PER_GROUP = 500
RANDOM_SEED = 1

GROUP_COLOR = {"ps20nm": "tab:blue", "ps30nm": "tab:orange", "ps40nm": "tab:green",
               "ps50nm": "tab:red"}


def green_intensity(fin, fmetrics, group):
    """Summed intensity ('is') for GREEN particles (passes both filters) within the
    y-illumination band (see filter_config.py), trimmed."""
    cat, _ = pfsm.classify_particles(fin, fmetrics, group)
    ys = fin[group]["ys"][:]
    mask = (cat == pfsm.PASSES_BOTH) & (ys >= Y_ILLUMINATION_MIN) & (ys <= Y_ILLUMINATION_MAX)
    intensity = fin[group]["is"][:][mask]
    lo, hi = np.percentile(intensity, TRIM_PERCENTILES)
    return intensity[(intensity >= lo) & (intensity <= hi)]


def fit_ps_calibration(fin, fmetrics):
    medians = {}
    for group in PS_GROUPS_FOR_CALIBRATION:
        sixth_root = green_intensity(fin, fmetrics, group) ** (1 / 6)
        medians[group] = np.median(sixth_root)
    xs = np.array([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION])
    ys = np.array([medians[g] for g in PS_GROUPS_FOR_CALIBRATION])
    slope, intercept, r, _, _ = stats.linregress(xs, ys)
    return slope, intercept, r**2


def per_particle_sizing_precision(fin, fmetrics, slope, intercept):
    """For each PS group, apply the calibration to every individual particle (not
    just the population median) and report the spread of inferred sizes."""
    results = {}
    print("Per-particle sizing precision (calibration applied to each particle individually):")
    print(f"{'group':<10}{'nominal':>9}{'n':>7}{'mean(nm)':>11}{'std(nm)':>10}{'SEM(nm)':>10}")
    for group in PS_GROUPS_FOR_CALIBRATION:
        sixth_root = green_intensity(fin, fmetrics, group) ** (1 / 6)
        sizes_est = (sixth_root - intercept) / slope
        results[group] = sizes_est
        sem = sizes_est.std() / np.sqrt(len(sizes_est))
        print(f"{group:<10}{PS_NOMINAL_SIZES_NM[group]:>9}{len(sizes_est):>7}"
              f"{sizes_est.mean():>11.2f}{sizes_est.std():>10.2f}{sem:>10.3f}")
    print(
        "-> per-particle std is comparable to, or larger than, the 10nm spacing "
        "between standards: a single particle's inferred size is not reliable on "
        "its own. The tiny SEM values show the population MEDIAN is precise, but "
        "only because of the large particle counts -- it says nothing about any "
        "one particle."
    )
    return results


def estimate_noise_floor(fin):
    """Estimate the background/shot-noise floor on the summed-intensity ('is')
    measurement from thumbnail corner-pixel statistics, propagated over the
    SUMMATION_RADIUS_PX window used upstream."""
    rng = np.random.default_rng(RANDOM_SEED)
    pixel_stds = []
    for group in PS_GROUPS_FOR_CALIBRATION:
        thumbnails = fin[group]["thumbnails"]
        n = min(N_BACKGROUND_SAMPLES_PER_GROUP, thumbnails.shape[0])
        idx = rng.choice(thumbnails.shape[0], size=n, replace=False)
        for k in idx:
            t = thumbnails[k]
            corner = np.concatenate(
                [t[:8, :8].ravel(), t[-8:, -8:].ravel(), t[:8, -8:].ravel(), t[-8:, :8].ravel()]
            )
            pixel_stds.append(np.std(corner))
    pixel_std = np.median(pixel_stds)
    area_px = np.pi * SUMMATION_RADIUS_PX**2
    noise_floor_sum = pixel_std * np.sqrt(area_px)
    return pixel_std, noise_floor_sum


def print_snr_table(fin, fmetrics, noise_floor_sum):
    print()
    print(f"Background noise floor on summed intensity (r={SUMMATION_RADIUS_PX}px window): "
          f"{noise_floor_sum:.1f}  (3-sigma: {3 * noise_floor_sum:.1f})")
    print(f"{'group':<10}{'median intensity':>18}{'SNR (median / 3-sigma noise)':>32}")
    for group in PS_GROUPS_FOR_CALIBRATION:
        intensity = green_intensity(fin, fmetrics, group)
        median_intensity = np.median(intensity)
        snr = median_intensity / (3 * noise_floor_sum)
        print(f"{group:<10}{median_intensity:>18.0f}{snr:>32.1f}")
    print(
        "-> signal is well above the sensor/background noise floor even for the "
        "smallest standard: detection sensitivity is NOT the bottleneck here."
    )


def plot_sizing_precision(fin, fmetrics, per_particle_sizes):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    bins = np.linspace(-40, 100, 100)
    for group in PS_GROUPS_FOR_CALIBRATION:
        sizes_est = per_particle_sizes[group]
        ax.hist(sizes_est, bins=bins, alpha=0.5, density=True, color=GROUP_COLOR[group],
                label=f"{group}: {sizes_est.mean():.1f} +/- {sizes_est.std():.1f} nm (n={len(sizes_est)})")
        ax.axvline(PS_NOMINAL_SIZES_NM[group], color=GROUP_COLOR[group], ls="--", lw=1.5)

    ax.set_xlabel("Per-particle inferred size (nm)")
    ax.set_ylabel("Density")
    ax.set_title(
        "Per-particle sizing precision vs. 10nm standard spacing\n"
        "(dashed lines = nominal size; distributions this wide relative to the\n"
        "10nm spacing mean individual particles can't be reliably sized)"
    )
    ax.legend(fontsize=8)
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "minimum_particle_size_sizing_precision.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nSaved {out_path}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")

    slope, intercept, r2 = fit_ps_calibration(fin, fmetrics)
    print(f"Calibration (ps20/30/40/50): intensity^(1/6) = {slope:.5f}*size_nm + {intercept:.5f}  "
          f"R^2={r2:.4f}\n")

    per_particle_sizes = per_particle_sizing_precision(fin, fmetrics, slope, intercept)

    pixel_std, noise_floor_sum = estimate_noise_floor(fin)
    print(f"\nMedian background pixel std: {pixel_std:.3f}")
    print_snr_table(fin, fmetrics, noise_floor_sum)

    plot_sizing_precision(fin, fmetrics, per_particle_sizes)
