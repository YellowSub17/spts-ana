"""
Diagnostic plots for the focus/shape filtering work (see compute_focus_shape_metrics.py).

Particles go through two filtering stages, and we color-code by where they land:

  RED    "fails old focus test" -- flags==False. Rejected by the original two-radius
         intensity-ratio focus test (src/comborun/filtering.py: filter_focused).
  YELLOW "passes old test, fails solidity test" -- flags==True but solidity <
         SOLIDITY_THRESHOLD. These are the bean/kidney/double-lobed particles that
         the intensity-ratio test alone doesn't catch (not a single convex blob).
  GREEN  "passes both filters" -- flags==True and solidity >= SOLIDITY_THRESHOLD.
         The final analysis set used downstream for intensity/sizing.

(A small number of flags==True particles have no measurable solidity -- the
thresholded blob was too small/degenerate to form a hull. These can't be confirmed
as a clean single blob either, so they're grouped with YELLOW/failed-solidity.)

Produces, per particle group in thumbnails.h5:
  1. Three example montages -- one per category above -- with colored image borders
     matching the category color, so it's visually obvious which stage a shown
     particle failed (or passed) at.
  2. Histograms of solidity, gauss_r2, eccen, and area, per group, as stacked
     RED/YELLOW/GREEN populations -- shows how well each metric separates the
     categories, and how consistent that separation is across groups/particle types.

Run compute_focus_shape_metrics.py first to generate data/focus_shape_metrics.h5.
Figures are written to the figures/ directory (created if it doesn't exist).
"""

import os

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from filter_config import SOLIDITY_THRESHOLD

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"

N_EXAMPLES = 10
RANDOM_SEED = 0

# category codes
FAILS_OLD_TEST = 0
PASSES_OLD_FAILS_SOLIDITY = 1
PASSES_BOTH = 2

CATEGORY_COLOR = {
    FAILS_OLD_TEST: "tab:red",
    PASSES_OLD_FAILS_SOLIDITY: "gold",
    PASSES_BOTH: "tab:green",
}
CATEGORY_LABEL = {
    FAILS_OLD_TEST: "RED -- fails old focus test (flags==False)",
    PASSES_OLD_FAILS_SOLIDITY: "YELLOW -- passes old test, fails solidity test",
    PASSES_BOTH: "GREEN -- passes both filters (used in analysis)",
}
CATEGORY_FILENAME = {
    FAILS_OLD_TEST: "examples_RED_fails_old_focus_test.png",
    PASSES_OLD_FAILS_SOLIDITY: "examples_YELLOW_fails_solidity_test.png",
    PASSES_BOTH: "examples_GREEN_passes_both_filters.png",
}


def classify_particles(fin, fmetrics, group):
    """Assign each particle in `group` to FAILS_OLD_TEST / PASSES_OLD_FAILS_SOLIDITY / PASSES_BOTH."""
    flags = fin[group]["flags"][:]
    solidity = fmetrics[group]["solidity"][:]
    solidity_valid = ~np.isnan(solidity)

    category = np.full(flags.shape, PASSES_OLD_FAILS_SOLIDITY, dtype=int)
    category[~flags] = FAILS_OLD_TEST
    # flags==True & solidity valid & high enough -> passes both
    category[flags & solidity_valid & (solidity >= SOLIDITY_THRESHOLD)] = PASSES_BOTH
    # flags==True & (solidity < threshold, or not measurable) -> fails solidity stage
    # (default value above already covers this; nothing more to set)
    return category, solidity


def plot_filter_stage_examples(fin, fmetrics, groups):
    """One example montage per category, with colored borders matching RED/YELLOW/GREEN."""
    rng = np.random.default_rng(RANDOM_SEED)

    for category in (FAILS_OLD_TEST, PASSES_OLD_FAILS_SOLIDITY, PASSES_BOTH):
        color = CATEGORY_COLOR[category]
        fig, axes = plt.subplots(len(groups), N_EXAMPLES, figsize=(2 * N_EXAMPLES, 2.1 * len(groups)))

        for i, group in enumerate(groups):
            cat, solidity = classify_particles(fin, fmetrics, group)
            thumbnails = fin[group]["thumbnails"]

            idx_pool = np.where(cat == category)[0]
            n = min(N_EXAMPLES, len(idx_pool))
            idx = rng.choice(idx_pool, size=n, replace=False) if n > 0 else []

            for j in range(N_EXAMPLES):
                ax = axes[i, j]
                if j < len(idx):
                    k = idx[j]
                    ax.imshow(thumbnails[k], cmap="gray")
                    s = solidity[k]
                    ax.set_title(f"s={s:.2f}" if not np.isnan(s) else "s=n/a", fontsize=7)
                ax.set_xticks([])
                ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_edgecolor(color)
                    spine.set_linewidth(3)
            axes[i, 0].set_ylabel(group, fontsize=10, rotation=90)

        plt.suptitle(CATEGORY_LABEL[category], color=color, fontweight="bold")
        plt.tight_layout()
        out_path = os.path.join(FIGURES_DIR, CATEGORY_FILENAME[category])
        plt.savefig(out_path, dpi=110)
        plt.close(fig)
        print(f"Saved {out_path}")


def plot_metric_distributions(fin, fmetrics, groups):
    """Histograms of solidity, gauss_r2, eccen, and area, stacked as RED/YELLOW/GREEN, per group."""
    metrics = [
        ("solidity", "Solidity (blob area / convex hull area)", lambda g: fmetrics[g]["solidity"][:]),
        ("gauss_r2", "2D Gaussian fit R^2", lambda g: fmetrics[g]["gauss_r2"][:]),
        ("eccen", "Eccentricity (from thumbnails.h5)", lambda g: fin[g]["eccen"][:]),
        ("area", "Blob area, log10(px) (from thumbnails.h5)", lambda g: np.log10(np.clip(fin[g]["area"][:], 1, None))),
    ]

    for metric_name, metric_label, getter in metrics:
        fig, axes = plt.subplots(2, 3, figsize=(15, 8))
        for ax, group in zip(axes.ravel(), groups):
            cat, _ = classify_particles(fin, fmetrics, group)
            values = getter(group)
            valid = ~np.isnan(values)

            lo = np.nanpercentile(values[valid], 0.5)
            hi = np.nanpercentile(values[valid], 99.5)
            bins = np.linspace(lo, hi, 60)

            for category in (FAILS_OLD_TEST, PASSES_OLD_FAILS_SOLIDITY, PASSES_BOTH):
                sel = valid & (cat == category)
                short_label = CATEGORY_LABEL[category].split(" -- ")[0]
                ax.hist(values[sel], bins=bins, alpha=0.5,
                        label=f"{short_label} (n={sel.sum()})",
                        color=CATEGORY_COLOR[category], density=True)
            if metric_name == "solidity":
                ax.axvline(SOLIDITY_THRESHOLD, color="k", ls="--", lw=1)
            ax.set_title(group)
            ax.legend(fontsize=6)

        plt.suptitle(f"{metric_label} distributions by group and filter stage")
        plt.tight_layout()
        out_path = os.path.join(FIGURES_DIR, f"distribution_{metric_name}.png")
        plt.savefig(out_path, dpi=110)
        plt.close(fig)
        print(f"Saved {out_path}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")
    groups = list(fin.keys())

    plot_filter_stage_examples(fin, fmetrics, groups)
    plot_metric_distributions(fin, fmetrics, groups)

    print("\nAll figures written to", FIGURES_DIR)
