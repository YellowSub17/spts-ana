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
    (ps_calibration_sizing.py), then histogrammed as a density.

Produces figures/dma_vs_optical_distributions.png.

Run ps_calibration_sizing.py and fit_dma_gaussian_peaks.py first (this script
reuses their filtering/calibration/fit-window logic by import).
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
from plot_dma_data import parse_dma_file, STEADY_STATE_SCAN_NUMBERS, DMA_FILE
from fit_dma_gaussian_peaks import DMA_SAMPLES_TO_FIT

FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"
THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"

TRIM_PERCENTILES = (2, 98)
PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm", "ps50nm"]
PS_NOMINAL_SIZES_NM = {"ps20nm": 20, "ps30nm": 30, "ps40nm": 40, "ps50nm": 50}

PLOT_XLIM = {
    "20nm PS": (0, 45),
    "30nm PS": (0, 65),
    "40nm PS": (0, 75),
    "50nm PS": (0, 95),
    "GroEL 1uM": (0, 35),
}


def optical_sixth_root(fin, fmetrics, group):
    cat, _ = pfsm.classify_particles(fin, fmetrics, group)
    ys = fin[group]["ys"][:]
    mask = (cat == pfsm.PASSES_BOTH) & (ys >= Y_ILLUMINATION_MIN) & (ys <= Y_ILLUMINATION_MAX)
    intensity = fin[group]["is"][:][mask]
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
    os.makedirs(FIGURES_DIR, exist_ok=True)

    header, diameters, distributions, footer = parse_dma_file(DMA_FILE)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")

    # PS calibration, fit on all four PS standards
    medians = {g: np.median(optical_sixth_root(fin, fmetrics, g)) for g in PS_GROUPS_FOR_CALIBRATION}
    xs = np.array([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION])
    ys = np.array([medians[g] for g in PS_GROUPS_FOR_CALIBRATION])
    slope, intercept, r, _, _ = stats.linregress(xs, ys)

    def size_from_sixth_root(v):
        return (v - intercept) / slope

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

        ax.set_xlim(*xlim)
        ax.set_xlabel("Diameter (nm)")
        ax.set_ylabel("probability density")
        ax.set_title(sample, fontsize=10)
        ax.legend(fontsize=8)

    plt.suptitle("DMA vs. optical (scattering) particle size distributions")
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "dma_vs_optical_distributions.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved {out_path}")
