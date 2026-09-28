"""
Directly overlays the DMA size distribution and the optical (light-scattering)
inferred size distribution for each sample with clean DMA data, on the same
diameter axis -- a shape-for-shape comparison rather than just summary
mu/sigma numbers (see fit_dma_gaussian_peaks.py for the Gaussian-fit version
of this comparison).

Shows the five samples fit in fit_dma_gaussian_peaks.py: 20nm/30nm/40nm/50nm
PS and GroEL. Ferritin is excluded -- its DMA block never settles into a
clean single-mode distribution, so there's nothing trustworthy to overlay
against.

Both curves are normalized to a probability density over linear diameter (nm)
so they're comparable in shape regardless of the very different units/
concentration scales of a DMA scan vs. an optical particle count:
  - DMA: dW/dlogDp is converted to a linear-diameter density via
    dW/dD = (dW/dlogDp) / (D * ln10), then normalized to integrate to 1 over
    the fit window (same window as fit_dma_gaussian_peaks.py, chosen to
    exclude the low-diameter contamination artifact).
  - Optical: the per-particle intensity of every GREEN, y-cropped particle in
    the group is converted to an inferred size via the PS calibration
    (ps_calibration_sizing.py: through zero, fitted on ps30/40/50nm), then
    histogrammed as a density.

The dotted vertical line marks the optical floor in size units: the size whose
intensity equals the dimmest GREEN PS particle. Below it the optical
distribution is censored (only the bright tail of a population is detected).

Produces figures/dma_vs_optical_distributions<tag>.png.

Run ps_calibration_sizing.py and fit_dma_gaussian_peaks.py first (this script
reuses their filtering/calibration/fit-window logic by import).
"""

import argparse
import os

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_focus_shape_metrics as pfsm
from filter_config import Y_ILLUMINATION_MIN, Y_ILLUMINATION_MAX
from plot_dma_data import parse_dma_file, STEADY_STATE_SCAN_NUMBERS, DMA_FILE
from fit_dma_gaussian_peaks import DMA_SAMPLES_TO_FIT
from ps_calibration_sizing import PS_GROUPS_FOR_CALIBRATION, PS_NOMINAL_SIZES_NM, fit_through_origin

FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"
THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"

TRIM_PERCENTILES = (2, 98)

PLOT_XLIM = {
    "20nm PS": (0, 45),
    "30nm PS": (0, 65),
    "40nm PS": (0, 75),
    "50nm PS": (0, 95),
    "GroEL 1uM": (0, 35),
}


def green_intensity(fin, fmetrics, group):
    cat, _ = pfsm.classify_particles(fin, fmetrics, group)
    ys = fin[group]["ys"][:]
    mask = (cat == pfsm.PASSES_BOTH) & (ys >= Y_ILLUMINATION_MIN) & (ys <= Y_ILLUMINATION_MAX)
    return fin[group]["is"][:][mask]


def optical_sixth_root(fin, fmetrics, group):
    intensity = green_intensity(fin, fmetrics, group)
    lo, hi = np.percentile(intensity, TRIM_PERCENTILES)
    trimmed = intensity[(intensity >= lo) & (intensity <= hi)]
    return trimmed ** (1 / 6)


def dma_linear_density(diameters, distribution, window):
    """Convert dW/dlogDp to a linear-diameter density, normalized to integrate
    to 1 over `window`. Returns only the (diameter, density) points inside
    `window` -- outside it is dominated by the low-diameter contamination
    artifact (see plot_dma_data.py docstring) and isn't meaningful to show."""
    density = distribution / (diameters * np.log(10))
    lo, hi = window
    sel = (diameters >= lo) & (diameters <= hi)
    area = np.trapezoid(density[sel], diameters[sel])
    return diameters[sel], density[sel] / area


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DMA vs optical size distribution overlays")
    parser.add_argument("--thumbnails", default=THUMBNAILS_H5)
    parser.add_argument("--shape-metrics", default=SHAPE_METRICS_H5,
                        help="focus_shape_metrics h5 computed from the same --thumbnails file")
    parser.add_argument("--tag", default="", help='appended to the output filename, e.g. "_t20_ctm"')
    args = parser.parse_args()

    os.makedirs(FIGURES_DIR, exist_ok=True)

    header, diameters, distributions, footer = parse_dma_file(DMA_FILE)

    fin = h5py.File(args.thumbnails, "r")
    fmetrics = h5py.File(args.shape_metrics, "r")

    # PS calibration: through zero, fitted on ps30/40/50 (see ps_calibration_sizing.py)
    medians = {g: np.median(optical_sixth_root(fin, fmetrics, g)) for g in PS_GROUPS_FOR_CALIBRATION}
    slope, r2 = fit_through_origin([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION],
                                   [medians[g] for g in PS_GROUPS_FOR_CALIBRATION])

    def size_from_sixth_root(v):
        return v / slope

    # optical floor in size units: dimmest GREEN particle across the PS groups
    floor_intensity = min(green_intensity(fin, fmetrics, g).min() for g in PS_NOMINAL_SIZES_NM)
    floor_size = size_from_sixth_root(floor_intensity ** (1 / 6))
    print(f"Calibration: intensity^(1/6) = {slope:.5f} * size_nm (R^2={r2:.4f}); "
          f"optical floor {floor_intensity:.0f} -> {floor_size:.1f} nm")

    samples = list(DMA_SAMPLES_TO_FIT.keys())
    fig, axes = plt.subplots(1, len(samples), figsize=(5 * len(samples), 4.5))

    for ax, sample in zip(axes, samples):
        info = DMA_SAMPLES_TO_FIT[sample]
        group = info["optical_group"]
        window = info["window"]

        # DMA distribution, averaged over steady-state scans, as a linear density
        scan_indices = [n - 1 for n in STEADY_STATE_SCAN_NUMBERS[sample]]
        avg_distribution = distributions[:, scan_indices].mean(axis=1)
        d, dma_density = dma_linear_density(diameters, avg_distribution, window)

        # optical inferred-size distribution, as a density histogram
        sixth_root = optical_sixth_root(fin, fmetrics, group)
        optical_sizes = size_from_sixth_root(sixth_root)
        xlim = PLOT_XLIM[sample]
        bins = np.linspace(*xlim, 60)

        ax.hist(optical_sizes, bins=bins, density=True, color="tab:green", alpha=0.5,
                label=f"optical (n={len(optical_sizes)})")
        ax.plot(d, dma_density, color="tab:blue", lw=1.8, label="DMA")
        ax.axvline(floor_size, color="0.4", ls=":", lw=1.2, label=f"optical floor ({floor_size:.0f} nm)")
        print(f"{sample}: optical median {np.median(optical_sizes):.1f} nm (n={len(optical_sizes)}), "
              f"DMA mode {d[np.argmax(dma_density)]:.1f} nm")

        ax.set_xlim(*xlim)
        ax.set_xlabel("Diameter (nm)")
        ax.set_ylabel("probability density")
        ax.set_title(sample, fontsize=10)
        ax.legend(fontsize=8)

    plt.suptitle("DMA vs. optical (scattering) particle size distributions"
                 + (f" ({args.tag.strip('_').replace('_', ', ')})" if args.tag else ""))
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, f"dma_vs_optical_distributions{args.tag}.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved {out_path}")
