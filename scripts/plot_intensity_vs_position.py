"""
Checks whether particle intensity depends on its (x, y) position in the frame --
a signature of non-uniform laser illumination (e.g. a bright stripe/band) rather
than genuine particle-to-particle variation. Motivated by the DMA cross-check
(see chat / ANALYSIS_SUMMARY_p1.md): the DMA shows the PS stocks are physically
tight (GSD ~1.3-1.6), so the ~15-17% intensity spread seen optically has to be
coming from the optical measurement itself, not real size polydispersity.

Uses GREEN particles (flags & solidity >= SOLIDITY_THRESHOLD, see
filter_config.py) from ps20/30/40nm, pooled together after normalizing each
particle's intensity^(1/6) by its own group's median -- this removes the
particle-size-dependent baseline so a spatial illumination pattern (which should
affect all groups the same way, since it's a property of the optical setup, not
the particle) becomes visible pooled across groups with much better statistics
than any single group alone.

Produces figures/intensity_vs_xy.png:
  - 2D binned heatmap of median normalized intensity^(1/6) across (x, y)
  - 1D profile of median normalized intensity^(1/6) vs x (marginalized over y)
  - 1D profile of median normalized intensity^(1/6) vs y (marginalized over x)

This is a diagnostic/illustrative step -- it does NOT apply any filtering. Run
after confirming a bright x/y band exists, then follow up with the region
actually being cropped out for the calibration.
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

PS_GROUPS = ["ps20nm", "ps30nm", "ps40nm"]
N_XY_BINS = 25


def pooled_normalized_intensity(fin, fmetrics):
    """Pool GREEN particles from ps20/30/40, each particle's intensity^(1/6)
    normalized by its own group's median, along with its (x, y) position."""
    all_x, all_y, all_norm_sixth_root = [], [], []
    for group in PS_GROUPS:
        cat, _ = pfsm.classify_particles(fin, fmetrics, group)
        green = cat == pfsm.PASSES_BOTH
        intensity = fin[group]["is"][:][green]
        xs = fin[group]["xs"][:][green]
        ys = fin[group]["ys"][:][green]

        sixth_root = intensity ** (1 / 6)
        norm = sixth_root / np.median(sixth_root)

        all_x.append(xs)
        all_y.append(ys)
        all_norm_sixth_root.append(norm)

    return np.concatenate(all_x), np.concatenate(all_y), np.concatenate(all_norm_sixth_root)


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")

    x, y, norm_intensity = pooled_normalized_intensity(fin, fmetrics)
    print(f"Pooled n={len(x)} GREEN particles from {PS_GROUPS}")
    print(f"x range: {x.min():.1f}-{x.max():.1f}   y range: {y.min():.1f}-{y.max():.1f}")

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))

    # 2D heatmap
    ax = axes[0]
    x_bins = np.linspace(x.min(), x.max(), N_XY_BINS + 1)
    y_bins = np.linspace(y.min(), y.max(), N_XY_BINS + 1)
    x_idx = np.clip(np.digitize(x, x_bins) - 1, 0, N_XY_BINS - 1)
    y_idx = np.clip(np.digitize(y, y_bins) - 1, 0, N_XY_BINS - 1)
    grid = np.full((N_XY_BINS, N_XY_BINS), np.nan)
    for i in range(N_XY_BINS):
        for j in range(N_XY_BINS):
            sel = (x_idx == i) & (y_idx == j)
            if sel.sum() >= 3:
                grid[j, i] = np.median(norm_intensity[sel])
    im = ax.imshow(grid, origin="lower", aspect="auto", cmap="RdBu_r", vmin=0.7, vmax=1.3,
                    extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]])
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    ax.set_title("Median normalized intensity^(1/6) vs (x, y)")
    plt.colorbar(im, ax=ax, fraction=0.046)

    # 1D profile vs x
    ax = axes[1]
    x_centers = 0.5 * (x_bins[:-1] + x_bins[1:])
    x_medians = [np.median(norm_intensity[x_idx == i]) if (x_idx == i).sum() >= 3 else np.nan
                 for i in range(N_XY_BINS)]
    ax.plot(x_centers, x_medians, marker="o")
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("median normalized intensity^(1/6)")
    ax.set_title("Intensity vs x (marginalized over y)")

    # 1D profile vs y
    ax = axes[2]
    y_centers = 0.5 * (y_bins[:-1] + y_bins[1:])
    y_medians = [np.median(norm_intensity[y_idx == j]) if (y_idx == j).sum() >= 3 else np.nan
                 for j in range(N_XY_BINS)]
    ax.plot(y_centers, y_medians, marker="o", color="tab:orange")
    ax.axhline(1.0, color="gray", ls="--", lw=1)
    ax.set_xlabel("y (px)")
    ax.set_ylabel("median normalized intensity^(1/6)")
    ax.set_title("Intensity vs y (marginalized over x)")

    plt.suptitle("Intensity vs. position, pooled GREEN ps20/30/40nm "
                  "(normalized by each group's median)")
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "intensity_vs_xy.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nSaved {out_path}")
