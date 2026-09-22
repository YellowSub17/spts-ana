"""
Per-sample scatter plots of particle position (x, y) vs. intensity^(1/6), for
each of the four PS calibration groups. Complements plot_intensity_vs_position.py
(which pools ps20/30/40nm together into binned medians to find the illumination
band) by showing the raw, unbinned, per-particle data for each group
individually -- useful for checking the y-illumination effect holds up
sample-by-sample rather than only in the pooled/binned view.

Uses GREEN particles (flags & solidity >= SOLIDITY_THRESHOLD, see
filter_config.py) with NO y-position crop applied -- the point here is to see
the position-dependence directly, so cropping it out first would defeat the
purpose.

Produces figures/position_vs_intensity_scatter.png: one row per PS group
(ps20/30/40/50nm), with x-position vs intensity^(1/6) in the left column and
y-position vs intensity^(1/6) in the right column.
"""

import os

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_focus_shape_metrics as pfsm
from filter_config import Y_ILLUMINATION_MIN, Y_ILLUMINATION_MAX

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"

PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm", "ps50nm"]


def green_xy_intensity(fin, fmetrics, group):
    """(x, y, intensity^(1/6)) for GREEN particles, no y-crop applied."""
    cat, _ = pfsm.classify_particles(fin, fmetrics, group)
    green = cat == pfsm.PASSES_BOTH
    xs = fin[group]["xs"][:][green]
    ys = fin[group]["ys"][:][green]
    sixth_root = fin[group]["is"][:][green] ** (1 / 6)
    return xs, ys, sixth_root


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")

    fig, axes = plt.subplots(len(PS_GROUPS_FOR_CALIBRATION), 2,
                              figsize=(11, 3 * len(PS_GROUPS_FOR_CALIBRATION)))

    for row, group in enumerate(PS_GROUPS_FOR_CALIBRATION):
        xs, ys, sixth_root = green_xy_intensity(fin, fmetrics, group)

        ax = axes[row, 0]
        ax.scatter(xs, sixth_root, s=6, alpha=0.4, color="tab:blue")
        ax.set_xlabel("x (px)")
        ax.set_ylabel("intensity^(1/6)")
        ax.set_title(f"{group}: x vs intensity (n={len(xs)})", fontsize=9)

        ax = axes[row, 1]
        ax.scatter(ys, sixth_root, s=6, alpha=0.4, color="tab:orange")
        ax.axvspan(Y_ILLUMINATION_MIN, Y_ILLUMINATION_MAX, color="green", alpha=0.08,
                   label="y-illumination band (used in calibration)")
        ax.set_xlabel("y (px)")
        ax.set_ylabel("intensity^(1/6)")
        ax.set_title(f"{group}: y vs intensity (n={len(ys)})", fontsize=9)
        ax.legend(fontsize=7)

    plt.suptitle("Particle position vs. intensity^(1/6), per PS calibration sample")
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "position_vs_intensity_scatter.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved {out_path}")
