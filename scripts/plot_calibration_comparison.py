"""
One-figure summary of the PS calibration after switching spts peak centering from
center_of_mass (_com) to center_to_max (_ctm):

  (a) t20: _com vs _ctm. Ordinary (free-intercept) line fits on ps30/40/50nm, drawn
      back to size 0, to show that _ctm passes through the origin on its own while
      _com has a large intercept.
  (b) _ctm: t10 vs t20 -- the detection threshold barely changes the calibration.
  (c) _ctm t20: ps30/40/50nm fitted with a free intercept vs forced through zero.

Points are GREEN medians of intensity^(1/6) (same filtering as
ps_calibration_sizing.py: flags & solidity & y-illumination crop, 2-98% trimmed);
error bars are bootstrapped 95% CIs on the median. ps20nm is drawn as an open
marker and never fitted (a typical 20nm bead is at the detection floor).

Writes figures/calibration_com_vs_ctm_summary.png.
"""

import os

import h5py
import numpy as np
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ps_calibration_sizing import (PS_GROUPS_FOR_CALIBRATION, PS_GROUPS_NOT_FITTED, PS_NOMINAL_SIZES_NM,
                                   fit_through_origin, focused_sixth_root_intensity)

DATA_DIR = "/Users/pat/Documents/work/spts-ana/data"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"
DATASETS = {
    "t20_com": ("thumbnails_t20_com.h5", "focus_shape_metrics_t20.h5"),
    "t20_ctm": ("thumbnails_t20_ctm.h5", "focus_shape_metrics_t20_ctm.h5"),
    "t10_ctm": ("thumbnails_t10_ctm.h5", "focus_shape_metrics_t10_ctm.h5"),
}
N_BOOT = 2000
RANDOM_SEED = 0

# reference palette (dataviz skill), categorical slots in order
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK_2, INK_MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"

ALL_PS = PS_GROUPS_NOT_FITTED + PS_GROUPS_FOR_CALIBRATION


def medians_with_ci(name):
    thumbs, metrics = DATASETS[name]
    rng = np.random.default_rng(RANDOM_SEED)
    out = {}
    with h5py.File(os.path.join(DATA_DIR, thumbs), "r") as fin, h5py.File(os.path.join(DATA_DIR, metrics), "r") as fm:
        for g in ALL_PS:
            v, _ = focused_sixth_root_intensity(fin, fm, g)
            boot = np.median(rng.choice(v, size=(N_BOOT, len(v)), replace=True), axis=1)
            out[g] = (np.median(v), *np.percentile(boot, [2.5, 97.5]), len(v))
    return out


def free_fit(res):
    x = np.array([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION], float)
    y = np.array([res[g][0] for g in PS_GROUPS_FOR_CALIBRATION])
    f = stats.linregress(x, y)
    return f.slope, f.intercept, f.intercept_stderr, f.rvalue**2


def zero_fit(res):
    x = [PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION]
    y = [res[g][0] for g in PS_GROUPS_FOR_CALIBRATION]
    return fit_through_origin(x, y)


def draw_points(ax, res, color, label, dx=0.0):
    for g in ALL_PS:
        med, lo, hi, _ = res[g]
        fitted = g in PS_GROUPS_FOR_CALIBRATION
        ax.errorbar(PS_NOMINAL_SIZES_NM[g] + dx, med, yerr=[[med - lo], [hi - med]], fmt="o", ms=8,
                    color=color, mfc=color if fitted else "white", mec=color, mew=2, elinewidth=2,
                    capsize=3, zorder=4, label=label if g == ALL_PS[-1] else None)


def style(ax, title):
    ax.set_title(title, loc="left", fontsize=11, color=INK)
    ax.set_xlim(-2, 55)
    ax.set_ylim(-1.0, 10)
    ax.set_xlabel("Nominal PS diameter (nm)", color=INK_2)
    ax.axhline(0, color=INK_MUTED, lw=0.8, zorder=1)
    ax.axvline(0, color=INK_MUTED, lw=0.8, zorder=1)
    ax.grid(color=GRID, lw=0.8, zorder=0)
    ax.tick_params(colors=INK_2)
    for s in ax.spines.values():
        s.set_visible(False)


if __name__ == "__main__":
    res = {name: medians_with_ci(name) for name in DATASETS}
    xfit = np.linspace(0, 55, 100)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True)
    fig.patch.set_facecolor("white")

    # (a) com vs ctm, free-intercept fits extended to 0
    ax = axes[0]
    for name, color, lab, dx, text_xy in (("t20_com", ORANGE, "t20, centre of mass", -0.7, (4, 2.1)),
                                          ("t20_ctm", BLUE, "t20, centre to max", 0.7, (7, -0.75))):
        s, b, b_se, r2 = free_fit(res[name])
        ax.plot(xfit, s * xfit + b, color=color, lw=2, ls="--", zorder=3)
        draw_points(ax, res[name], color, lab, dx=dx)
        ax.annotate(f"intercept {b:+.2f} ± {b_se:.2f}", xy=(0, b), xytext=text_xy,
                    color=INK, fontsize=9, arrowprops=dict(arrowstyle="-", color=color, lw=1))
        ax.plot(0, b, "o", ms=6, color=color, zorder=5)
    style(ax, "(a) t20: centre of mass vs centre to max\nfree-intercept fit on ps30–50, drawn to 0 (offset ±0.7 nm)")
    ax.set_ylabel("Median intensity$^{1/6}$", color=INK_2)

    # (b) ctm t10 vs t20
    ax = axes[1]
    for name, color, lab, dx, ls in (("t20_ctm", BLUE, "t20, centre to max", -0.7, "-"),
                                     ("t10_ctm", AQUA, "t10, centre to max", 0.7, (0, (4, 3)))):
        k, r2 = zero_fit(res[name])
        ax.plot(xfit, k * xfit, color=color, lw=2, ls=ls, zorder=3)
        draw_points(ax, res[name], color, f"{lab}: k = {k:.4f}", dx=dx)
    diffs = [100 * (res["t10_ctm"][g][0] / res["t20_ctm"][g][0] - 1) for g in ALL_PS]
    ax.text(54, -0.6, "t10 vs t20 median: " + ", ".join(f"{PS_NOMINAL_SIZES_NM[g]} nm {d:+.1f}%" for g, d in zip(ALL_PS, diffs)),
            ha="right", fontsize=8.5, color=INK_2)
    style(ax, "(b) centre to max: t10 vs t20\nfit through zero on ps30–50 (points offset ±0.7 nm)")

    # (c) ctm t20: free vs through zero
    ax = axes[2]
    s, b, b_se, r2f = free_fit(res["t20_ctm"])
    k, r2z = zero_fit(res["t20_ctm"])
    ax.plot(xfit, s * xfit + b, color=ORANGE, lw=2, ls="--", zorder=3,
            label=f"free:  {s:.4f}·D {b:+.2f}   R² = {r2f:.3f}")
    ax.plot(xfit, k * xfit, color=BLUE, lw=2, zorder=3, label=f"through 0:  {k:.4f}·D   R² = {r2z:.3f}")
    draw_points(ax, res["t20_ctm"], INK_2, "t20, centre to max")
    for g in PS_GROUPS_NOT_FITTED:
        med = res["t20_ctm"][g][0]
        ax.annotate(f"ps20 not fitted\n(line predicts {k * PS_NOMINAL_SIZES_NM[g]:.2f})", xy=(PS_NOMINAL_SIZES_NM[g], med),
                    xytext=(PS_NOMINAL_SIZES_NM[g] - 17, med + 1.6), fontsize=8.5, color=INK_2,
                    arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=1))
    style(ax, "(c) t20 centre to max: ps30–50 fit\nfree intercept vs forced through zero")

    for ax in axes:
        leg = ax.legend(loc="upper left", fontsize=8.5, frameon=False, labelcolor=INK)

    fig.text(0.01, 0.005, "Points: GREEN medians (focus flag & solidity ≥ 0.8 & y 30–167, 2–98% trimmed); "
             "error bars: bootstrapped 95% CI of the median. Open markers (ps20nm) are never fitted.",
             fontsize=8.5, color=INK_MUTED)
    plt.tight_layout(rect=(0, 0.03, 1, 1))
    out = os.path.join(FIGURES_DIR, "calibration_com_vs_ctm_summary.png")
    plt.savefig(out, dpi=150, facecolor="white")
    print(f"Saved {out}")
    for name in DATASETS:
        s, b, b_se, r2f = free_fit(res[name]); k, r2z = zero_fit(res[name])
        print(f"{name}: free {s:.5f}*D {b:+.3f} (se {b_se:.3f}) R2 {r2f:.4f} | through 0 {k:.5f}*D R2 {r2z:.4f} | "
              + "  ".join(f"{g} {res[name][g][0]:.3f} [{res[name][g][1]:.3f},{res[name][g][2]:.3f}] n={res[name][g][3]}" for g in ALL_PS))
