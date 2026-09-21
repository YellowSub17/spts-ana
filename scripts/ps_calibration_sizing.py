"""
Sixth-root intensity calibration (polystyrene spheres) and sizing of ferritin / GroEL.

Background
----------
Physical basis: for a Rayleigh scatterer, scattered intensity ~ volume^2, so
intensity^(1/6) is linear in particle diameter. Using the polystyrene (PS) size
ladder as standards, we fit that line and then read off an "equivalent size" for
the protein particles from their measured intensity.

Filtering used before any of the intensity stats below: flags (existing focus flag
in thumbnails.h5) AND solidity >= 0.80 (see compute_focus_shape_metrics.py) to
reject unfocused and non-convex/bean-shaped particles. Within that filtered set, the
top/bottom 2% of summed intensity is additionally trimmed per group to reduce the
influence of aggregates/debris and residual detection-threshold noise on the median.

Findings from this analysis (as of this writing):
  - ps20nm, ps30nm, ps40nm give an excellent linear fit on their own: R^2 ~ 0.989.
  - ps50nm sits well below where that line predicts -- its median/mean summed
    intensity is *lower* than ps40nm's, which is unphysical for a bigger PS bead.
    This isn't explained by saturation (only 1-2 particles are near the 16-bit
    ceiling) or by aggregate outliers (trimming the tails barely moves the median).
    Working hypothesis (per lab notes): the ps50nm stock is old/degraded. ps50nm is
    therefore EXCLUDED from the calibration fit; only ps20/30/40 are used.
  - Ferritin comes out around ~17 nm equivalent size, GroEL around ~9 nm -- i.e.
    GroEL scores *smaller* than ferritin on this scale, despite being the physically
    larger complex (~14 nm GroEL barrel vs. ~12-13 nm ferritin shell).

IMPORTANT CAVEAT -- refractive index / optical contrast, not yet corrected for:
    This calibration curve is built entirely on polystyrene, which has a high,
    uniform refractive index. Scattered intensity for a small particle actually
    depends on volume^2 * (optical contrast)^2, where the contrast is set by the
    particle's refractive index relative to the surrounding medium -- not on volume
    alone. Polystyrene beads all share the same contrast, so for them intensity^(1/6)
    is a clean proxy for size. Ferritin and GroEL do NOT share PS's contrast:
      - Ferritin carries a dense iron-oxide mineral core, giving it much higher
        optical contrast per unit volume than a pure-protein particle -- it will
        scatter like a bigger-than-it-really-is PS bead.
      - GroEL is a hollow, low-density protein-only complex, closer in index to the
        surrounding buffer -- it will scatter like a smaller-than-it-really-is PS
        bead.
    So the "nm" values this script reports for ferritin/GroEL are POLYSTYRENE-
    EQUIVALENT scattering sizes, not true physical diameters, and the apparent
    GroEL < ferritin ordering is likely this contrast effect rather than a real
    size measurement. To get true physical sizes, this fit needs a refractive-index
    correction: scale each protein's intensity by (contrast_PS / contrast_protein)^2
    before inverting the calibration line, using literature (or measured) refractive
    index / index-increment values for polystyrene, ferritin, and GroEL. That
    correction is NOT implemented here yet -- treat size_nm_ps_equivalent below as
    provisional until it is.

Output
------
Prints the calibration fit, per-group median sixth-root intensity, and the
PS-equivalent size estimate (with 16-84th percentile spread) for ferritin and GroEL.
Saves a comparison plot to ps_calibration_and_protein_sizing.png in the working
directory.

Run compute_focus_shape_metrics.py first to generate data/focus_shape_metrics.h5.
"""

import os

import h5py
import numpy as np
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from filter_config import SOLIDITY_THRESHOLD

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"
OUTPUT_PLOT = os.path.join(FIGURES_DIR, "ps_calibration_and_protein_sizing.png")

TRIM_PERCENTILES = (2, 98)

PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm"]  # ps50nm excluded, see docstring
PS_NOMINAL_SIZES_NM = {"ps20nm": 20, "ps30nm": 30, "ps40nm": 40, "ps50nm": 50}
PROTEIN_GROUPS = ["ferri", "groel"]


def focused_sixth_root_intensity(fin, fmetrics, group):
    """Sixth-root of summed intensity for particles passing flags & solidity filter,
    with the top/bottom TRIM_PERCENTILES of intensity trimmed."""
    flags = fin[group]["flags"][:]
    solidity = fmetrics[group]["solidity"][:]
    valid = ~np.isnan(solidity)
    mask = flags & valid & (solidity >= SOLIDITY_THRESHOLD)

    intensity = fin[group]["is"][:][mask]
    lo, hi = np.percentile(intensity, TRIM_PERCENTILES)
    trimmed = intensity[(intensity >= lo) & (intensity <= hi)]
    return trimmed ** (1 / 6), mask.sum()


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")

    # --- fit the PS calibration line on ps20/30/40 ---
    medians = {}
    for group in PS_GROUPS_FOR_CALIBRATION:
        sixth_root, n_before_trim = focused_sixth_root_intensity(fin, fmetrics, group)
        medians[group] = np.median(sixth_root)
        print(f"{group}: n={n_before_trim} (after trim {len(sixth_root)})  "
              f"median intensity^(1/6) = {medians[group]:.4f}")

    xs = np.array([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION])
    ys = np.array([medians[g] for g in PS_GROUPS_FOR_CALIBRATION])
    slope, intercept, r, _, _ = stats.linregress(xs, ys)
    print(f"\nCalibration (ps20/30/40 only): "
          f"intensity^(1/6) = {slope:.5f} * size_nm + {intercept:.5f}   R^2={r**2:.4f}")

    def size_from_sixth_root(sixth_root_value):
        return (sixth_root_value - intercept) / slope

    # sanity check: ps50nm, excluded from the fit, for reference
    ps50_sixth_root, _ = focused_sixth_root_intensity(fin, fmetrics, "ps50nm")
    ps50_median = np.median(ps50_sixth_root)
    print(f"\n[reference, excluded from fit] ps50nm: median intensity^(1/6)={ps50_median:.4f} "
          f"-> back-calculated size={size_from_sixth_root(ps50_median):.1f} nm "
          f"(nominal 50; low value supports old/degraded-sample hypothesis)")

    # --- size the proteins against that line ---
    protein_results = {}
    print()
    for group in PROTEIN_GROUPS:
        sixth_root, n_before_trim = focused_sixth_root_intensity(fin, fmetrics, group)
        protein_results[group] = sixth_root
        median_sixth_root = np.median(sixth_root)
        p16, p84 = np.percentile(sixth_root, [16, 84])
        size_est = size_from_sixth_root(median_sixth_root)
        size_lo, size_hi = size_from_sixth_root(p16), size_from_sixth_root(p84)
        print(f"{group}: n={n_before_trim} (after trim {len(sixth_root)})  "
              f"median intensity^(1/6)={median_sixth_root:.3f}  "
              f"-> PS-equivalent size ~ {size_est:.1f} nm "
              f"(16-84th percentile: {size_lo:.1f}-{size_hi:.1f} nm)")

    # --- plot ---
    fig, ax = plt.subplots(figsize=(7, 5.5))
    xfit = np.linspace(0, 55, 100)
    ax.plot(xfit, slope * xfit + intercept, "k--", label=f"PS calibration (R^2={r**2:.3f})")
    ax.scatter(xs, ys, color="tab:blue", s=70, zorder=3, label="PS 20/30/40 (calibration)")
    ax.scatter([50], [ps50_median], color="gray", marker="x", s=70,
               label="PS50 (excluded/degraded)")

    colors = {"ferri": "tab:green", "groel": "tab:purple"}
    for group in PROTEIN_GROUPS:
        sixth_root = protein_results[group]
        median_sixth_root = np.median(sixth_root)
        size_est = size_from_sixth_root(median_sixth_root)
        p16, p84 = np.percentile(sixth_root, [16, 84])
        size_lo, size_hi = size_from_sixth_root(p16), size_from_sixth_root(p84)
        ax.scatter([size_est], [median_sixth_root], color=colors[group], s=90,
                   marker="D", zorder=4, label=f"{group}: ~{size_est:.1f} nm (PS-equivalent)")
        ax.errorbar([size_est], [median_sixth_root],
                    xerr=[[size_est - size_lo], [size_hi - size_est]],
                    color=colors[group], capsize=4)

    ax.set_xlabel("Size (nm) -- PS-equivalent for ferritin/GroEL, see docstring caveat")
    ax.set_ylabel("Intensity^(1/6)")
    ax.legend(fontsize=8)
    ax.set_title("PS sixth-root calibration and protein sizing")
    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=120)
    print(f"\nSaved {OUTPUT_PLOT}")
