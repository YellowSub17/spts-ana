"""
Fits a Gaussian (in linear diameter, nm) to the main peak of each clean DMA
scan. The resulting mu/sigma per sample is used by
plot_ps_calibration_sizing.py to draw the DMA-measured physical size spread
as a horizontal error bar on the PS calibration curve, next to the
optical-derived spread -- a direct visual comparison of DMA-measured physical
size uncertainty against the optical sizing pipeline.

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

Run directly, this produces figures/dma_gaussian_fits.png (the fits
themselves, for a sanity check). compute_all_fits() is also imported by
plot_ps_calibration_sizing.py to overlay the fitted sigmas elsewhere.
"""

import os

import numpy as np
from scipy.optimize import curve_fit
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from plot_dma_data import parse_dma_file, STEADY_STATE_SCAN_NUMBERS, DMA_FILE

FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"

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


def compute_all_fits(diameters, distributions):
    """Fit each sample in DMA_SAMPLES_TO_FIT, returning
    {sample: {amplitude, mu, sigma, offset, avg_distribution, scan_numbers, window}}."""
    fits = {}
    for sample, info in DMA_SAMPLES_TO_FIT.items():
        scan_numbers = STEADY_STATE_SCAN_NUMBERS[sample]
        scan_indices = [n - 1 for n in scan_numbers]
        avg_distribution = distributions[:, scan_indices].mean(axis=1)
        amplitude, mu, sigma, offset = fit_peak(diameters, avg_distribution, info["window"])
        fits[sample] = dict(amplitude=amplitude, mu=mu, sigma=sigma, offset=offset,
                             avg_distribution=avg_distribution,
                             scan_numbers=scan_numbers, window=info["window"])
    return fits


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    header, diameters, distributions, footer = parse_dma_file(DMA_FILE)

    fits = compute_all_fits(diameters, distributions)
    for sample, fit in fits.items():
        print(f"{sample}: fit over scans {fit['scan_numbers']}, window {fit['window']}nm -> "
              f"mu={fit['mu']:.2f}nm, sigma={fit['sigma']:.2f}nm")
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
