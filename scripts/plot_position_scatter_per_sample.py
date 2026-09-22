"""
Per-sample scatter plots of particle position (x, y) vs. intensity^(1/6), one
figure per PS calibration group. Complements plot_intensity_vs_position.py
(which pools ps20/30/40nm together into binned medians to find the
illumination band) by showing the raw, unbinned, per-particle data for each
group individually.

Plots ALL particles in each group (not just GREEN), colored by the combined
filter (see plot_focus_shape_metrics.py for the RED/YELLOW/GREEN
classification): GREEN = passes both filters (kept), RED = fails either one
(the YELLOW category -- passes flags but fails solidity -- is folded into RED
here, since this plot is about the final keep/reject decision, not which
stage a particle failed at).

No y-illumination crop is applied to the plotted data, and no shaded band is
drawn -- the point of this plot is to see the raw position-dependence
directly.

Produces one figure per group: figures/pos/position_vs_intensity_<group>.png,
each with x-position vs intensity^(1/6) on the left and y-position vs
intensity^(1/6) on the right.
"""

import os

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_focus_shape_metrics as pfsm

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"
POS_DIR = os.path.join(FIGURES_DIR, "pos")

PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm", "ps50nm"]

KEEP_COLOR = pfsm.CATEGORY_COLOR[pfsm.PASSES_BOTH]  # "tab:green"
REJECT_COLOR = pfsm.CATEGORY_COLOR[pfsm.FAILS_OLD_TEST]  # "tab:red"


def xy_intensity_and_keep(fin, fmetrics, group):
    """(x, y, intensity^(1/6), keep) for every particle in group -- keep is
    True for GREEN (passes both filters), False otherwise (RED+YELLOW folded
    together)."""
    cat, _ = pfsm.classify_particles(fin, fmetrics, group)
    keep = cat == pfsm.PASSES_BOTH
    xs = fin[group]["xs"][:]
    ys = fin[group]["ys"][:]
    sixth_root = fin[group]["is"][:] ** (1 / 6)
    return xs, ys, sixth_root, keep


if __name__ == "__main__":
    os.makedirs(POS_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")

    for group in PS_GROUPS_FOR_CALIBRATION:
        xs, ys, sixth_root, keep = xy_intensity_and_keep(fin, fmetrics, group)
        colors = np.where(keep, KEEP_COLOR, REJECT_COLOR)

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

        ax = axes[0]
        ax.scatter(xs, sixth_root, s=6, alpha=0.4, c=colors)
        ax.set_xlabel("x (px)")
        ax.set_ylabel("intensity^(1/6)")
        ax.set_title(f"{group}: x vs intensity", fontsize=10)

        ax = axes[1]
        ax.scatter(ys, sixth_root, s=6, alpha=0.4, c=colors)
        ax.set_xlabel("y (px)")
        ax.set_ylabel("intensity^(1/6)")
        ax.set_title(f"{group}: y vs intensity", fontsize=10)

        ax.scatter([], [], color=KEEP_COLOR, label=f"kept (n={keep.sum()})")
        ax.scatter([], [], color=REJECT_COLOR, label=f"filtered (n={(~keep).sum()})")
        ax.legend(fontsize=8)

        plt.suptitle(f"{group}: particle position vs. intensity^(1/6)")
        plt.tight_layout()
        out_path = os.path.join(POS_DIR, f"position_vs_intensity_{group}.png")
        plt.savefig(out_path, dpi=120)
        plt.close(fig)
        print(f"Saved {out_path}")
