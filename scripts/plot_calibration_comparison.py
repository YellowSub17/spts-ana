"""
Four figures summarising the PS calibration after switching spts peak centering
from center_of_mass (_com) to center_to_max (_ctm):

  (a) t20: _com vs _ctm. Ordinary (free-intercept) line fits on ps30/40/50nm, drawn
      back to size 0, to show that _ctm passes through the origin on its own while
      _com has a large intercept.
  (b) _ctm: t10 vs t20 -- the detection threshold barely changes the calibration.
  (c) _ctm t20: ps30/40/50nm fitted with a free intercept vs forced through zero.
  (d) _ctm t20: through-zero fit on ps30/40/50nm vs including ps20nm in the fit.

Points are GREEN medians of intensity^(1/6) (same filtering as
ps_calibration_sizing.py: flags & solidity & y-illumination crop, 2-98% trimmed);
error bars are bootstrapped 95% CIs on the median. ps20nm is drawn as an open
marker and never fitted (a typical 20nm bead is at the detection floor).

Writes (in figures/):
  (a) calibration_t20_com_vs_ctm.png
  (b) calibration_ctm_t10_vs_t20.png
  (c) calibration_ctm_t20_free_vs_zero.png
  (d) calibration_ctm_t20_with_without_ps20.png
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
    ax.set_ylabel("Median intensity$^{1/6}$", color=INK_2)
    ax.axhline(0, color=INK_MUTED, lw=0.8, zorder=1)
    ax.axvline(0, color=INK_MUTED, lw=0.8, zorder=1)
    ax.grid(color=GRID, lw=0.8, zorder=0)
    ax.tick_params(colors=INK_2)
    for s in ax.spines.values():
        s.set_visible(False)


FOOTNOTE = ("Points: GREEN medians (focus flag & solidity ≥ 0.8 & y 30–167, 2–98% trimmed); "
            "error bars: bootstrapped 95% CI of the median.\nOpen markers (ps20nm) are never fitted: "
            "a typical 20 nm bead sits at the detection floor.")


def new_figure():
    fig, ax = plt.subplots(figsize=(7.5, 6))
    fig.patch.set_facecolor("white")
    return fig, ax


def save(fig, ax, filename, footnote=FOOTNOTE):
    ax.legend(loc="upper left", fontsize=9, frameon=False, labelcolor=INK)
    fig.text(0.01, 0.01, footnote, fontsize=8, color=INK_MUTED)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    out = os.path.join(FIGURES_DIR, filename)
    fig.savefig(out, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    res = {name: medians_with_ci(name) for name in DATASETS}
    xfit = np.linspace(0, 55, 100)

    # (a) com vs ctm, free-intercept fits extended to 0
    fig, ax = new_figure()
    for name, color, lab, dx, text_xy in (("t20_com", ORANGE, "t20, centre of mass", -0.7, (4, 2.1)),
                                          ("t20_ctm", BLUE, "t20, centre to max", 0.7, (7, -0.75))):
        s, b, b_se, r2 = free_fit(res[name])
        ax.plot(xfit, s * xfit + b, color=color, lw=2, ls="--", zorder=3)
        draw_points(ax, res[name], color, f"{lab}: {s:.4f}·D {b:+.2f}", dx=dx)
        ax.annotate(f"intercept {b:+.2f} ± {b_se:.2f}", xy=(0, b), xytext=text_xy,
                    color=INK, fontsize=9, arrowprops=dict(arrowstyle="-", color=color, lw=1))
        ax.plot(0, b, "o", ms=6, color=color, zorder=5)
    style(ax, "t20: centre of mass vs centre to max\nfree-intercept fit on ps30–50, drawn back to 0 (points offset ±0.7 nm)")
    save(fig, ax, "calibration_t20_com_vs_ctm.png")

    # (b) ctm t10 vs t20
    fig, ax = new_figure()
    for name, color, lab, dx, ls in (("t20_ctm", BLUE, "t20, centre to max", -0.7, "-"),
                                     ("t10_ctm", AQUA, "t10, centre to max", 0.7, (0, (4, 3)))):
        k, r2 = zero_fit(res[name])
        ax.plot(xfit, k * xfit, color=color, lw=2, ls=ls, zorder=3)
        draw_points(ax, res[name], color, f"{lab}: k = {k:.4f}", dx=dx)
    diffs = [100 * (res["t10_ctm"][g][0] / res["t20_ctm"][g][0] - 1) for g in ALL_PS]
    ax.text(54, -0.6, "t10 vs t20 median: " + ", ".join(f"{PS_NOMINAL_SIZES_NM[g]} nm {d:+.1f}%" for g, d in zip(ALL_PS, diffs)),
            ha="right", fontsize=8.5, color=INK_2)
    style(ax, "Centre to max: t10 vs t20\nfit through zero on ps30–50 (points offset ±0.7 nm)")
    save(fig, ax, "calibration_ctm_t10_vs_t20.png")

    # (c) ctm t20: free vs through zero
    fig, ax = new_figure()
    s, b, b_se, r2f = free_fit(res["t20_ctm"])
    k, r2z = zero_fit(res["t20_ctm"])
    ax.plot(xfit, s * xfit + b, color=ORANGE, lw=2, ls="--", zorder=3,
            label=f"free intercept:  {s:.4f}·D {b:+.2f}   R² = {r2f:.3f}")
    ax.plot(xfit, k * xfit, color=BLUE, lw=2, zorder=3, label=f"through zero:  {k:.4f}·D   R² = {r2z:.3f}")
    draw_points(ax, res["t20_ctm"], INK_2, "t20, centre to max")
    for g in PS_GROUPS_NOT_FITTED:
        med = res["t20_ctm"][g][0]
        ax.annotate(f"ps20 not fitted\n(line predicts {k * PS_NOMINAL_SIZES_NM[g]:.2f})", xy=(PS_NOMINAL_SIZES_NM[g], med),
                    xytext=(PS_NOMINAL_SIZES_NM[g] - 17, med + 1.6), fontsize=8.5, color=INK_2,
                    arrowprops=dict(arrowstyle="-", color=INK_MUTED, lw=1))
    style(ax, "t20 centre to max: ps30–50 fit\nfree intercept vs forced through zero")
    save(fig, ax, "calibration_ctm_t20_free_vs_zero.png")

    # (d) ctm t20: through-zero fit without vs with ps20
    fig, ax = new_figure()
    t20 = res["t20_ctm"]
    k30, r2_30 = zero_fit(t20)
    x_all = [PS_NOMINAL_SIZES_NM[g] for g in ALL_PS]
    k20, r2_20 = fit_through_origin(x_all, [t20[g][0] for g in ALL_PS])
    ax.plot(xfit, k30 * xfit, color=BLUE, lw=2, zorder=3,
            label=f"ps30–50 (used):  {k30:.4f}·D   R² = {r2_30:.3f}")
    ax.plot(xfit, k20 * xfit, color=ORANGE, lw=2, ls="--", zorder=3,
            label=f"ps20–50:  {k20:.4f}·D   R² = {r2_20:.3f}")
    draw_points(ax, t20, INK_2, "t20, centre to max")
    for g in ALL_PS:
        med, x = t20[g][0], PS_NOMINAL_SIZES_NM[g]
        # right-and-below for 30/40nm; left-and-above for 20nm and 50nm, which would otherwise sit on the lines
        right_side = g not in (ALL_PS[0], ALL_PS[-1])
        ax.text(x + 1.2 if right_side else x - 1.2, med + (-0.55 if right_side else 0.4),
                f"reads {med / k30:.1f} / {med / k20:.1f} nm", fontsize=8, color=INK_2,
                ha="left" if right_side else "right")
    ax.text(54, -0.6, f"including ps20: slope {100 * (k20 / k30 - 1):+.1f}%, so every size reads "
            f"{100 * (1 - k30 / k20):.1f}% smaller", ha="right", fontsize=8.5, color=INK_2)
    style(ax, "t20 centre to max: fit through zero\nwith vs without ps20nm (labels: size read off ps30–50 / ps20–50 line)")
    save(fig, ax, "calibration_ctm_t20_with_without_ps20.png",
         footnote=FOOTNOTE.split("\n")[0] + "\nOpen marker (ps20nm) is fitted only by the dashed ps20–50 line; "
         "a typical 20 nm bead sits at the detection floor.")

    for name in DATASETS:
        s, b, b_se, r2f = free_fit(res[name]); k, r2z = zero_fit(res[name])
        print(f"{name}: free {s:.5f}*D {b:+.3f} (se {b_se:.3f}) R2 {r2f:.4f} | through 0 {k:.5f}*D R2 {r2z:.4f}")
