"""
Plots for the PS sixth-root intensity calibration (see ps_calibration_sizing.py
for the underlying numbers/fit, including ferritin/GroEL protein sizing, which
is not plotted here -- this figure covers the PS standards only).

Uses only the GREEN particles -- flags==True & solidity>=SOLIDITY_THRESHOLD, i.e.
the particles that pass both filtering stages (see plot_focus_shape_metrics.py for
the RED/YELLOW/GREEN classification) -- with the top/bottom 2% of intensity
additionally trimmed to reduce the influence of aggregates/debris on the median.

Produces two figures:
  1. figures/dist/distribution_intensity_green.png
     Histograms of summed intensity (log-scale) for the GREEN particles, one panel
     per group (ps20/30/40/50nm, ferri, groel).
  2. figures/ps_calibration_size_vs_sixthroot.png
     Median intensity^(1/6) vs. nominal size for the polystyrene spheres, with a
     linear fit through all four PS standards (20/30/40/50nm). Also overlays the
     DMA-fitted Gaussian sigma (see fit_dma_gaussian_peaks.py) as a horizontal
     error bar at each PS sample with clean DMA data, alongside the optical
     spread -- a direct visual comparison of the two measurements.

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
from plot_dma_data import parse_dma_file, DMA_FILE
from fit_dma_gaussian_peaks import DMA_SAMPLES_TO_FIT, compute_all_fits

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"
DIST_DIR = os.path.join(FIGURES_DIR, "dist")

TRIM_PERCENTILES = (2, 98)
PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm", "ps50nm"]
PS_NOMINAL_SIZES_NM = {"ps20nm": 20, "ps30nm": 30, "ps40nm": 40, "ps50nm": 50}


def green_intensity(fin, fmetrics, group):
    """Summed intensity ('is') for GREEN particles (passes both filters) within the
    y-illumination band (see filter_config.py), trimmed."""
    cat, _ = pfsm.classify_particles(fin, fmetrics, group)
    ys = fin[group]["ys"][:]
    mask = (cat == pfsm.PASSES_BOTH) & (ys >= Y_ILLUMINATION_MIN) & (ys <= Y_ILLUMINATION_MAX)
    intensity = fin[group]["is"][:][mask]
    lo, hi = np.percentile(intensity, TRIM_PERCENTILES)
    return intensity[(intensity >= lo) & (intensity <= hi)]


def plot_intensity_distributions(fin, fmetrics, groups):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, group in zip(axes.ravel(), groups):
        intensity = green_intensity(fin, fmetrics, group)
        ax.hist(np.log10(intensity), bins=60, color=pfsm.CATEGORY_COLOR[pfsm.PASSES_BOTH])
        ax.set_title(f"{group} (n={len(intensity)})")
        ax.set_xlabel("log10(summed intensity)")
    plt.suptitle("Summed intensity distributions -- GREEN particles only (passes both filters)")
    plt.tight_layout()
    out_path = os.path.join(DIST_DIR, "distribution_intensity_green.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_size_vs_sixthroot(fin, fmetrics, diameters, distributions):
    # calibration fit on all four PS standards
    medians = {}
    spreads = {}
    for group in PS_GROUPS_FOR_CALIBRATION:
        sixth_root = green_intensity(fin, fmetrics, group) ** (1 / 6)
        medians[group] = np.median(sixth_root)
        spreads[group] = np.percentile(sixth_root, [16, 84])

    xs = np.array([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION])
    ys = np.array([medians[g] for g in PS_GROUPS_FOR_CALIBRATION])
    slope, intercept, r, _, _ = stats.linregress(xs, ys)

    def size_from_sixth_root(v):
        return (v - intercept) / slope

    fig, ax = plt.subplots(figsize=(7, 5.5))
    xfit = np.linspace(0, 55, 100)
    ax.plot(xfit, slope * xfit + intercept, "k--", label=f"PS calibration fit (R^2={r**2:.4f})")

    # PS points used in the fit, with 16-84th percentile error bars
    # (plot_position tracks each group's (x, y) on this plot, for the DMA overlay below)
    plot_position = {}
    for group in PS_GROUPS_FOR_CALIBRATION:
        lo, hi = spreads[group]
        plot_position[group] = (PS_NOMINAL_SIZES_NM[group], medians[group])
        ax.errorbar([PS_NOMINAL_SIZES_NM[group]], [medians[group]],
                    yerr=[[medians[group] - lo], [hi - medians[group]]],
                    fmt="o", color="tab:blue", capsize=4, zorder=3)
    ax.scatter([], [], color="tab:blue", label="PS 20/30/40/50nm (calibration)")

    # DMA-fitted Gaussian sigma, drawn as a horizontal error bar at each fitted
    # PS sample's plot position (see fit_dma_gaussian_peaks.py)
    dma_fits = compute_all_fits(diameters, distributions)
    first = True
    for sample, fit in dma_fits.items():
        group = DMA_SAMPLES_TO_FIT[sample]["optical_group"]
        if group not in plot_position:
            continue
        x_center, y_center = plot_position[group]
        ax.errorbar([x_center], [y_center], xerr=[[fit["sigma"]], [fit["sigma"]]],
                    fmt="none", color="tab:red", capsize=5, lw=2, zorder=5,
                    label="DMA Gaussian sigma" if first else None)

        # where the DMA +/- sigma bounds (evaluated on the fit line) fall on
        # this sample's vertical (intensity) error bar
        x_lo, x_hi = x_center - fit["sigma"], x_center + fit["sigma"]
        y_lo, y_hi = slope * x_lo + intercept, slope * x_hi + intercept
        ax.plot([x_center, x_center], [y_lo, y_hi], marker="_", markersize=16, markeredgewidth=3,
                 linestyle="none", color="tab:green", zorder=6,
                 label="DMA sigma bounds on intensity" if first else None)
        first = False

    ax.set_xlabel("Size (nm)")
    ax.set_ylabel("Median intensity^(1/6)")
    ax.legend(fontsize=8)
    ax.set_title("Size vs. median intensity^(1/6): PS calibration & DMA cross-check", fontsize=11)
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "ps_calibration_size_vs_sixthroot.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved {out_path}")
    print(f"Calibration: intensity^(1/6) = {slope:.5f} * size_nm + {intercept:.5f}  R^2={r**2:.4f}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(DIST_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")
    groups = list(fin.keys())
    header, diameters, distributions, footer = parse_dma_file(DMA_FILE)

    plot_intensity_distributions(fin, fmetrics, groups)
    plot_size_vs_sixthroot(fin, fmetrics, diameters, distributions)

    print("\nAll figures written to", FIGURES_DIR)
