"""
Fits a Gaussian (in linear diameter, nm) to the main peak of each clean DMA
scan, and overlays the resulting size std as a horizontal error bar on the PS
calibration curve, next to the existing optical-derived (vertical) spread --
a direct visual comparison of DMA-measured physical size uncertainty against
the optical sizing pipeline.

Fits five samples (see plot_dma_data.py: STEADY_STATE_SCAN_NUMBERS / the
run-log sample mapping): 20nm PS (scan 10), 30nm PS (scans 8-9), 40nm PS
(scans 6-7), 50nm PS (scan 3), GroEL (scans 13-14). Ferritin is skipped --
its DMA block never settles into a clean distribution, even in its last scan.

30nm PS's raw Mode/Median statistics look contaminated too, but the real
population is recoverable: it shows up as a second, broader peak (~28nm)
underneath the low-diameter contamination spike, consistent across both
scans in the block, with a clean falling tail out to ~90nm -- so it's fit
with a window that isolates that peak instead of being skipped.

For samples with more than one steady-state scan, the fit is done on the
average of those scans' distributions.

Produces figures/dma_gaussian_fits.png (the fits themselves, for a sanity
check) and figures/ps_calibration_with_dma_errorbars.png (the calibration
curve with DMA sigma drawn as horizontal error bars).
"""

import os

import h5py
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import plot_focus_shape_metrics as pfsm
from filter_config import Y_ILLUMINATION_MIN, Y_ILLUMINATION_MAX
from plot_dma_data import parse_dma_file, SCAN_SAMPLE_LABELS, STEADY_STATE_SCAN_NUMBERS, DMA_FILE

FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"
THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
SHAPE_METRICS_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"

TRIM_PERCENTILES = (2, 98)
PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm", "ps50nm"]
PS_NOMINAL_SIZES_NM = {"ps20nm": 20, "ps30nm": 30, "ps40nm": 40, "ps50nm": 50}

# sample -> (optical group name in thumbnails.h5, fit window in nm around the peak)
# windows chosen by eye from dma_distributions_grid.png -- wide enough to cover
# the whole peak, narrow enough to exclude the low-diameter contamination rise
# and the noisy high-diameter tail.
DMA_SAMPLES_TO_FIT = {
    "20nm PS": {"optical_group": "ps20nm", "window": (10, 35)},
    "30nm PS": {"optical_group": "ps30nm", "window": (18, 42)},
    "40nm PS": {"optical_group": "ps40nm", "window": (30, 52)},
    "50nm PS": {"optical_group": "ps50nm", "window": (45, 75)},
    "GroEL 1uM": {"optical_group": "groel", "window": (13, 24)},
}
SKIPPED_SAMPLES = ["Ferritin"]  # DMA data not purged, see docstring


def gaussian(x, amplitude, mu, sigma, offset):
    return offset + amplitude * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def fit_peak(diameters, distribution, window):
    lo, hi = window
    sel = (diameters >= lo) & (diameters <= hi)
    x, y = diameters[sel], distribution[sel]

    amplitude0 = y.max() - y.min()
    mu0 = x[np.argmax(y)]
    sigma0 = (hi - lo) / 6
    offset0 = y.min()

    popt, _ = curve_fit(
        gaussian, x, y, p0=[amplitude0, mu0, sigma0, offset0],
        bounds=([0, lo, 0.1, 0], [np.inf, hi, hi - lo, y.max()]),
        maxfev=5000,
    )
    return popt  # amplitude, mu, sigma, offset


def optical_sixth_root(fin, fmetrics, group):
    cat, _ = pfsm.classify_particles(fin, fmetrics, group)
    ys = fin[group]["ys"][:]
    mask = (cat == pfsm.PASSES_BOTH) & (ys >= Y_ILLUMINATION_MIN) & (ys <= Y_ILLUMINATION_MAX)
    intensity = fin[group]["is"][:][mask]
    lo, hi = np.percentile(intensity, TRIM_PERCENTILES)
    trimmed = intensity[(intensity >= lo) & (intensity <= hi)]
    return trimmed ** (1 / 6)


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    header, diameters, distributions, footer = parse_dma_file(DMA_FILE)

    # --- fit each clean sample's averaged steady-state distribution ---
    fits = {}
    for sample, info in DMA_SAMPLES_TO_FIT.items():
        scan_numbers = STEADY_STATE_SCAN_NUMBERS[sample]
        scan_indices = [n - 1 for n in scan_numbers]
        avg_distribution = distributions[:, scan_indices].mean(axis=1)
        amplitude, mu, sigma, offset = fit_peak(diameters, avg_distribution, info["window"])
        fits[sample] = dict(amplitude=amplitude, mu=mu, sigma=sigma, offset=offset,
                             avg_distribution=avg_distribution,
                             scan_numbers=scan_numbers, window=info["window"])
        print(f"{sample}: fit over scans {scan_numbers}, window {info['window']}nm -> "
              f"mu={mu:.2f}nm, sigma={sigma:.2f}nm")

    for sample in SKIPPED_SAMPLES:
        print(f"{sample}: skipped -- DMA scans not purged (see docstring)")

    # --- sanity-check plot: DMA data + Gaussian fit, one panel per sample ---
    fig, axes = plt.subplots(1, len(fits), figsize=(4.5 * len(fits), 4))
    for ax, (sample, fit) in zip(axes, fits.items()):
        ax.plot(diameters, fit["avg_distribution"], color="tab:blue", lw=1,
                label="DMA (avg steady scans)")
        xfit = np.linspace(*fit["window"], 200)
        ax.plot(xfit, gaussian(xfit, fit["amplitude"], fit["mu"], fit["sigma"], fit["offset"]),
                "r--", lw=1.5, label="Gaussian fit")
        ax.axvline(fit["mu"], color="gray", ls=":", lw=1)
        ax.set_xlim(0, max(fit["window"][1] * 1.3, 100))
        ax.set_title(f"{sample}\nmu={fit['mu']:.1f}nm sigma={fit['sigma']:.2f}nm", fontsize=9)
        ax.set_xlabel("Diameter (nm)")
        ax.legend(fontsize=7)
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "dma_gaussian_fits.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nSaved {out_path}")

    # --- PS calibration curve, with DMA sigma drawn as horizontal error bars ---
    fin = h5py.File(THUMBNAILS_H5, "r")
    fmetrics = h5py.File(SHAPE_METRICS_H5, "r")

    medians = {g: np.median(optical_sixth_root(fin, fmetrics, g)) for g in PS_GROUPS_FOR_CALIBRATION}
    xs = np.array([PS_NOMINAL_SIZES_NM[g] for g in PS_GROUPS_FOR_CALIBRATION])
    ys = np.array([medians[g] for g in PS_GROUPS_FOR_CALIBRATION])
    slope, intercept, r, _, _ = stats.linregress(xs, ys)

    def size_from_sixth_root(v):
        return (v - intercept) / slope

    fig, ax = plt.subplots(figsize=(8, 6))
    xfit = np.linspace(0, 90, 100)
    ax.plot(xfit, slope * xfit + intercept, "k--", lw=1, label=f"PS calibration (R^2={r**2:.4f})")

    # optical spread (vertical, y-error) for the four PS standards
    for group in PS_GROUPS_FOR_CALIBRATION:
        sixth_root = optical_sixth_root(fin, fmetrics, group)
        p16, p84 = np.percentile(sixth_root, [16, 84])
        med = medians[group]
        ax.errorbar([PS_NOMINAL_SIZES_NM[group]], [med],
                    yerr=[[med - p16], [p84 - med]],
                    fmt="o", color="tab:blue", capsize=4, zorder=3, lw=1)
    ax.scatter([], [], color="tab:blue", label="optical spread (16-84th %ile, vertical)")

    # GroEL, from the optical calibration inversion
    groel_sixth_root = optical_sixth_root(fin, fmetrics, "groel")
    groel_median = np.median(groel_sixth_root)
    groel_size = size_from_sixth_root(groel_median)
    p16, p84 = np.percentile(groel_sixth_root, [16, 84])
    ax.errorbar([groel_size], [groel_median],
                xerr=[[groel_size - size_from_sixth_root(p16)],
                      [size_from_sixth_root(p84) - groel_size]],
                fmt="D", color="tab:purple", capsize=4, zorder=3, lw=1,
                label="GroEL: optical spread (16-84th %ile, horizontal)")

    # DMA sigma, drawn as a horizontal error bar at the calibration point for
    # each fitted sample (PS points use their nominal size; GroEL uses its
    # optically-inverted size)
    dma_color = "tab:red"
    first = True
    for sample, fit in fits.items():
        group = DMA_SAMPLES_TO_FIT[sample]["optical_group"]
        if group in PS_GROUPS_FOR_CALIBRATION:
            x_center = PS_NOMINAL_SIZES_NM[group]
            y_center = medians[group]
        else:  # groel
            x_center = groel_size
            y_center = groel_median
        ax.errorbar([x_center], [y_center], xerr=[[fit["sigma"]], [fit["sigma"]]],
                    fmt="none", color=dma_color, capsize=5, lw=2, zorder=4,
                    label="DMA Gaussian sigma (horizontal)" if first else None)
        first = False

    ax.set_xlabel("Size (nm)")
    ax.set_ylabel("Median intensity^(1/6)")
    ax.set_title("PS calibration: optical spread (vertical/blue) vs. DMA Gaussian sigma (horizontal/red)")
    ax.legend(fontsize=8)
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "ps_calibration_with_dma_errorbars.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved {out_path}")
