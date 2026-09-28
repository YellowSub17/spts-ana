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
reject unfocused and non-convex/bean-shaped particles, AND a crop on y-position
(see filter_config.py / plot_intensity_vs_position.py) to remove a y-dependent
illumination gradient that otherwise biases group medians. Within that filtered
set, the top/bottom 2% of summed intensity is additionally trimmed per group to
reduce the influence of aggregates/debris and residual detection-threshold noise
on the median.

The calibration is fitted through zero on ps30/40/50nm only: Rayleigh scattering
gives intensity ~ D^6, so intensity^(1/6) = k * D with no intercept. ps20nm is
reported (and plotted) but left out of the fit: a typical 20nm bead sits at the
detection floor, so only its bright tail is detected and its median is biased
high. (The earlier 4-point affine fit with R^2 ~ 0.9997 came from spts
center_of_mass peak centering, which made the focus test reject most bright
ps40/ps50 particles -- use spts peak_centering = center_to_max output.)

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
PS-equivalent size estimate (with 16-84th percentile spread) for ferritin and
GroEL. For the plotted version of this calibration (with DMA error bars), see
plot_ps_calibration_sizing.py -> figures/ps_calibration_size_vs_sixthroot.png.

Run compute_focus_shape_metrics.py first to generate data/focus_shape_metrics.h5.
"""

import argparse

import h5py
import numpy as np

from filter_config import SOLIDITY_THRESHOLD, Y_ILLUMINATION_MIN, Y_ILLUMINATION_MAX

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"

TRIM_PERCENTILES = (2, 98)

PS_GROUPS_FOR_CALIBRATION = ["ps30nm", "ps40nm", "ps50nm"]
# shown alongside the calibration but not fitted (censored by the detection floor)
PS_GROUPS_NOT_FITTED = ["ps20nm"]
PS_NOMINAL_SIZES_NM = {"ps20nm": 20, "ps30nm": 30, "ps40nm": 40, "ps50nm": 50}
PROTEIN_GROUPS = ["ferri", "groel"]


def fit_through_origin(sizes_nm, sixth_roots):
    """Least-squares slope k for intensity^(1/6) = k * size_nm (no intercept).
    R^2 is 1 - SS_res/SS_tot about the mean, so it's directly comparable to an
    ordinary fit's R^2 (and can be low or negative if a line through zero fits badly)."""
    sizes_nm, sixth_roots = np.asarray(sizes_nm, float), np.asarray(sixth_roots, float)
    k = np.sum(sizes_nm * sixth_roots) / np.sum(sizes_nm**2)
    ss_res = np.sum((sixth_roots - k * sizes_nm) ** 2)
    ss_tot = np.sum((sixth_roots - sixth_roots.mean()) ** 2)
    return k, 1 - ss_res / ss_tot


def focused_sixth_root_intensity(fin, fmetrics, group, use_flags=True):
    """Sixth-root of summed intensity for particles passing flags & solidity filter
    and the y-illumination-band crop (see filter_config.py), with the top/bottom
    TRIM_PERCENTILES of intensity trimmed. use_flags=False skips the flags filter."""
    flags = fin[group]["flags"][:] if use_flags else np.ones(len(fin[group]["is"]), dtype=bool)
    solidity = fmetrics[group]["solidity"][:]
    ys = fin[group]["ys"][:]
    valid = ~np.isnan(solidity)
    mask = (
        flags
        & valid
        & (solidity >= SOLIDITY_THRESHOLD)
        & (ys >= Y_ILLUMINATION_MIN)
        & (ys <= Y_ILLUMINATION_MAX)
    )

    intensity = fin[group]["is"][:][mask]
    lo, hi = np.percentile(intensity, TRIM_PERCENTILES)
    trimmed = intensity[(intensity >= lo) & (intensity <= hi)]
    return trimmed ** (1 / 6), mask.sum()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PS sixth-root intensity calibration and protein sizing")
    parser.add_argument("--thumbnails", default=THUMBNAILS_H5)
    parser.add_argument("--shape-metrics", default=SHAPE_METRICS_H5,
                        help="focus_shape_metrics h5 computed from the same --thumbnails file")
    parser.add_argument("--no-flags", action="store_true", help="skip the flags focus filter (solidity + y crop only)")
    args = parser.parse_args()

    fin = h5py.File(args.thumbnails, "r")
    fmetrics = h5py.File(args.shape_metrics, "r")
    for group in fin:
        if len(fmetrics[group]["solidity"]) != len(fin[group]["is"]):
            raise SystemExit(f"{args.shape_metrics} doesn't match {args.thumbnails} ({group} row counts differ)")

    # --- fit the PS calibration line through zero on ps30/40/50 (ps20 reported, not fitted) ---
    medians = {}
    for group in PS_GROUPS_NOT_FITTED + PS_GROUPS_FOR_CALIBRATION:
        sixth_root, n_before_trim = focused_sixth_root_intensity(fin, fmetrics, group, use_flags=not args.no_flags)
        medians[group] = np.median(sixth_root)
        print(f"{group}: n={n_before_trim} (after trim {len(sixth_root)})  "
              f"median intensity^(1/6) = {medians[group]:.4f}")

    xs = np.array([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION])
    ys = np.array([medians[g] for g in PS_GROUPS_FOR_CALIBRATION])
    slope, r2 = fit_through_origin(xs, ys)
    print(f"\nCalibration (ps30/40/50, through zero): "
          f"intensity^(1/6) = {slope:.5f} * size_nm   R^2={r2:.4f}")

    def size_from_sixth_root(sixth_root_value):
        return sixth_root_value / slope

    for group in PS_GROUPS_NOT_FITTED:
        predicted = slope * PS_NOMINAL_SIZES_NM[group]
        print(f"{group} (not fitted): measured median {medians[group]:.4f} vs line {predicted:.4f} "
              f"-> reads as {size_from_sixth_root(medians[group]):.1f} nm "
              f"(expected median intensity {predicted**6:.0f})")

    # --- size the proteins against that line ---
    protein_results = {}
    print()
    for group in PROTEIN_GROUPS:
        sixth_root, n_before_trim = focused_sixth_root_intensity(fin, fmetrics, group, use_flags=not args.no_flags)
        protein_results[group] = sixth_root
        median_sixth_root = np.median(sixth_root)
        p16, p84 = np.percentile(sixth_root, [16, 84])
        size_est = size_from_sixth_root(median_sixth_root)
        size_lo, size_hi = size_from_sixth_root(p16), size_from_sixth_root(p84)
        print(f"{group}: n={n_before_trim} (after trim {len(sixth_root)})  "
              f"median intensity^(1/6)={median_sixth_root:.3f}  "
              f"-> PS-equivalent size ~ {size_est:.1f} nm "
              f"(16-84th percentile: {size_lo:.1f}-{size_hi:.1f} nm)")
