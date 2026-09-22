"""
Parses and plots the DMA (differential mobility analyzer) scan data in
data/AIM148.txt.

File format (TSI AIM software export, tab-separated, latin-1 encoded):
  - A block of scalar instrument/scan metadata (one row per field, one column
    per sample/scan), including Date and Start Time per scan.
  - A "Diameter Midpoint" table: one row per log-spaced diameter bin (nm), one
    column per scan, values are dw/dlogDp (number-weighted concentration).
  - A footer block of derived per-scan statistics: Median(nm), Mean(nm),
    Geo. Mean(nm), Mode(nm), Geo. Std. Dev., Total Concentration(#/cm3).

This file has 16 scans, all taken on 2026-03-30, spanning roughly 08:15-15:14.
Sample identity per scan number (from the run log) is in SCAN_SAMPLE_LABELS
below: scans 1-3 = 50nm PS, 4-7 = 40nm PS, 8-9 = 30nm PS, 10-11 = 20nm PS,
12-14 = GroEL 1uM, 15-16 = Ferritin. Within most blocks, the *first* scan(s)
after switching samples show a dominant peak pinned at the DMA's lower
measurement limit (~5.8nm) with the real population barely visible -- almost
certainly carryover/residual aerosol from the previous sample or drying buffer
salt, not yet purged from the lines. The *later* scan(s) in each block are the
more trustworthy steady-state measurement of that sample. 30nm PS and Ferritin
only got 2 scans each and neither one settles into a clean single-mode
distribution in that block -- both are likely still contaminated/not
equilibrated, so treat their DMA stats as unreliable.

Produces figures/dma_distributions_grid.png -- one panel per scan, each
annotated with its sample, start time, and Mode/Median/GSD from the footer
stats, for a closer look at each individual scan's shape (single mode vs.
multiple populations).

Also prints a summary table of Mode/Median/Mean/GSD/Total Concentration per
scan.
"""

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DMA_FILE = "/Users/pat/Documents/work/spts-ana/data/AIM148.txt"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"

PS_NOMINAL_SIZES_NM = [20, 30, 40, 50]

# scan number (1-indexed) -> sample identity, from the run log
SCAN_SAMPLE_LABELS = (
    ["50nm PS"] * 3
    + ["40nm PS"] * 4
    + ["30nm PS"] * 2
    + ["20nm PS"] * 2
    + ["GroEL 1uM"] * 3
    + ["Ferritin"] * 2
)

# the later scan(s) in each block, taken as the steady-state/purged measurement
# for that sample (see docstring -- earlier scans in a block are often still
# showing carryover from the previous sample)
STEADY_STATE_SCAN_NUMBERS = {
    "50nm PS": [3],
    "40nm PS": [6, 7],
    "30nm PS": [8, 9],  # neither looks fully purged -- treat with caution
    "20nm PS": [10],
    "GroEL 1uM": [13, 14],
    "Ferritin": [15, 16],  # neither looks fully purged -- treat with caution
}


def parse_dma_file(path):
    lines = open(path, encoding="latin-1").read().splitlines()

    header = {}
    data_start = data_end = None
    for i, line in enumerate(lines):
        parts = line.split("\t")
        key = parts[0].strip()
        if key == "Diameter Midpoint":
            data_start = i + 1
            continue
        if key == "Scan Up Time(s)":
            data_end = i
            break
        if data_start is None and len(parts) > 1:
            header[key] = parts[1:]

    n_scans = len(header.get("Sample #", []))

    diameters = []
    distributions = []
    for line in lines[data_start:data_end]:
        parts = line.split("\t")
        diameters.append(float(parts[0]))
        distributions.append([float(x) for x in parts[1 : 1 + n_scans]])
    diameters = np.array(diameters)
    distributions = np.array(distributions)  # (n_bins, n_scans)

    footer = {}
    for line in lines[data_end:]:
        parts = line.split("\t")
        key = parts[0].strip()
        if key and len(parts) > 1:
            footer[key] = parts[1 : 1 + n_scans]

    return header, diameters, distributions, footer


def print_summary_table(header, footer, n_scans):
    print(f"{'#':>3}{'sample':>12}{'time':>10}{'mode(nm)':>10}{'median(nm)':>12}"
          f"{'mean(nm)':>10}{'GSD':>7}{'total conc (#/cm3)':>20}")
    for i in range(n_scans):
        time = header["Start Time"][i]
        sample = SCAN_SAMPLE_LABELS[i]
        mode = float(footer["Mode(nm)"][i])
        median = float(footer["Median(nm)"][i])
        mean = float(footer["Mean(nm)"][i])
        gsd = float(footer["Geo. Std. Dev."][i])
        total_conc_key = [k for k in footer if k.startswith("Total Concentration")][0]
        total_conc = float(footer[total_conc_key][i])
        steady = i + 1 in STEADY_STATE_SCAN_NUMBERS.get(sample, [])
        flag = " *" if steady else ""
        print(f"{i + 1:>3}{sample:>12}{time:>10}{mode:>10.1f}{median:>12.1f}{mean:>10.1f}"
              f"{gsd:>7.2f}{total_conc:>20.3e}{flag}")
    print("  (* = steady-state/purged scan, used for the optical comparison)")


def plot_grid(diameters, distributions, header, footer, n_scans):
    ncols = 4
    nrows = int(np.ceil(n_scans / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows))
    for i in range(n_scans):
        ax = axes.ravel()[i]
        ax.plot(diameters, distributions[:, i], color="tab:blue", lw=1.2)
        ax.set_xscale("log")
        for size in PS_NOMINAL_SIZES_NM:
            ax.axvline(size, color="gray", ls="--", lw=0.6, zorder=0)
        mode = float(footer["Mode(nm)"][i])
        median = float(footer["Median(nm)"][i])
        gsd = float(footer["Geo. Std. Dev."][i])
        steady = i + 1 in STEADY_STATE_SCAN_NUMBERS.get(SCAN_SAMPLE_LABELS[i], [])
        star = " *" if steady else ""
        ax.set_title(f"#{i + 1} {SCAN_SAMPLE_LABELS[i]}{star} ({header['Start Time'][i]})\n"
                      f"mode={mode:.1f} median={median:.1f} GSD={gsd:.2f}", fontsize=8)
        ax.tick_params(labelsize=7)
    for j in range(n_scans, nrows * ncols):
        axes.ravel()[j].axis("off")

    plt.suptitle("DMA distributions, per scan")
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "dma_distributions_grid.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    os.makedirs(FIGURES_DIR, exist_ok=True)

    header, diameters, distributions, footer = parse_dma_file(DMA_FILE)
    n_scans = distributions.shape[1]
    print(f"Parsed {n_scans} scans, {len(diameters)} diameter bins "
          f"({diameters.min():.2f}-{diameters.max():.1f} nm)\n")

    print_summary_table(header, footer, n_scans)
    plot_grid(diameters, distributions, header, footer, n_scans)
