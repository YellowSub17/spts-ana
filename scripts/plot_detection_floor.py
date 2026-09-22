"""
Checks whether the nonzero intercept in the PS sixth-root calibration
(ps_calibration_sizing.py: intensity^(1/6) = 0.052*size_nm + 3.058) is
explained by the sensor/background noise floor, or by something else.

Forcing the calibration fit through the origin makes it dramatically worse
(R^2 drops from 0.9997 to -1.55), and the sensor/illumination noise floor
(~101 raw intensity units, 3-sigma ~303, see minimum_particle_size.py) is far
too small to produce the observed intercept -- reproducing it as a constant
additive background would require B ~ 4500, ~15x the measured noise floor,
and even then the fit is poor (a real background pedestal would need to be
constant across groups; matching each group's intercept-implied background
requires B to scale nearly as steeply as the signal itself, which a constant
background cannot do).

Instead, this checks the ENTIRE unfiltered detection population (every
particle EVER flagged as a hit, before any of our flags/solidity/y-crop
filtering) for a common, size-independent minimum intensity -- the signature
of a fixed hit-detection threshold in the acquisition pipeline. Finding: all
four PS groups bottom out at almost the same absolute intensity (~1000-1250,
roughly 3-4x the measured noise floor) regardless of nominal particle size.
A per-particle SNR-driven detection limit would shift with the particle's
true brightness; a fixed pipeline threshold would not -- this matches the
latter.

That fixed floor is a large fraction of the ps20nm group's own median
intensity (~22%) but a small fraction of ps50nm's (~4%), which is exactly the
size-dependent gradient needed to turn a true zero-intercept
intensity~diameter^6 relationship into an apparent nonzero-intercept affine
fit over the 20-50nm range: the acquisition threshold censors the dim tail of
the true distribution, pulling the median up more for dim (small) particles
than for bright (large) ones.

Produces figures/detection_floor.png: histograms (log-x) of ALL detected
intensities per PS group, with the common floor and the independently
estimated noise floor marked.
"""

import os

import h5py
import numpy as np
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from minimum_particle_size import estimate_noise_floor

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"

PS_GROUPS_FOR_CALIBRATION = ["ps20nm", "ps30nm", "ps40nm", "ps50nm"]
PS_NOMINAL_SIZES_NM = {"ps20nm": 20, "ps30nm": 30, "ps40nm": 40, "ps50nm": 50}


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    fin = h5py.File(THUMBNAILS_H5, "r")

    _, noise_floor_1sigma = estimate_noise_floor(fin)
    noise_floor_3sigma = 3 * noise_floor_1sigma
    print(f"Noise floor: 1-sigma={noise_floor_1sigma:.1f}, 3-sigma={noise_floor_3sigma:.1f}\n")

    print(f"{'group':<10}{'n_total':>9}{'min':>10}{'min/3sigma':>12}{'median':>12}{'floor%median':>14}")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, group in zip(axes.ravel(), PS_GROUPS_FOR_CALIBRATION):
        intensity = fin[group]["is"][:]
        intensity = intensity[intensity > 0]
        median = np.median(intensity)
        min_val = intensity.min()
        print(f"{group:<10}{len(intensity):>9}{min_val:>10.0f}{min_val / noise_floor_3sigma:>12.2f}"
              f"{median:>12.0f}{100 * min_val / median:>13.1f}%")

        ax.hist(np.log10(intensity), bins=80, color="tab:blue", alpha=0.8)
        ax.axvline(np.log10(min_val), color="black", ls="-", lw=1.5,
                   label=f"min detected = {min_val:.0f}")
        ax.axvline(np.log10(noise_floor_3sigma), color="tab:red", ls="--", lw=1.5,
                   label=f"3-sigma noise floor = {noise_floor_3sigma:.0f}")
        ax.set_title(f"{group} (nominal {PS_NOMINAL_SIZES_NM[group]}nm), n={len(intensity)}", fontsize=10)
        ax.set_xlabel("log10(summed intensity)")
        ax.legend(fontsize=7)

    plt.suptitle("All detections (unfiltered) per PS group: common minimum-intensity floor")
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "detection_floor.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nSaved {out_path}")
