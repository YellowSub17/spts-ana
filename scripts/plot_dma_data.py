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

This file has 16 scans, all taken on 2026-03-30, spanning roughly 08:15-15:14 --
a day of DMA runs, presumably of (some subset of) the same particle stocks used
in the optical thumbnails.h5 analysis (PS 20/30/40/50nm, ferritin, GroEL),
though the DMA scans aren't labeled with which stock is which -- only sample
number and time. Match scans to stocks using the run log / acquisition order,
and by eye against the known nominal sizes (dashed reference lines on the
overlay plot at 20/30/40/50 nm).

Produces two figures in figures/:
  1. dma_distributions_overlay.png -- all 16 scans' dw/dlogDp vs. diameter,
     log-x, colored by scan order, with dashed reference lines at the PS
     nominal sizes for comparison to the optical sizing work.
  2. dma_distributions_grid.png -- one panel per scan, each annotated with its
     start time and Mode/Median/GSD from the footer stats, for a closer look
     at each individual scan's shape (single mode vs. multiple populations).

Also prints a summary table of Mode/Median/Mean/GSD/Total Concentration per
scan.
"""

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

DMA_FILE = "/Users/pat/Documents/work/spts-ana/data/AIM148.txt"
FIGURES_DIR = "/Users/pat/Documents/work/spts-ana/figures"

PS_NOMINAL_SIZES_NM = [20, 30, 40, 50]


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
    print(f"{'#':>3}{'time':>10}{'mode(nm)':>10}{'median(nm)':>12}{'mean(nm)':>10}"
          f"{'GSD':>7}{'total conc (#/cm3)':>20}")
    for i in range(n_scans):
        time = header["Start Time"][i]
        mode = float(footer["Mode(nm)"][i])
        median = float(footer["Median(nm)"][i])
        mean = float(footer["Mean(nm)"][i])
        gsd = float(footer["Geo. Std. Dev."][i])
        total_conc_key = [k for k in footer if k.startswith("Total Concentration")][0]
        total_conc = float(footer[total_conc_key][i])
        print(f"{i + 1:>3}{time:>10}{mode:>10.1f}{median:>12.1f}{mean:>10.1f}"
              f"{gsd:>7.2f}{total_conc:>20.3e}")


def plot_overlay(diameters, distributions, header, n_scans):
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = cm.viridis(np.linspace(0, 1, n_scans))
    for i in range(n_scans):
        ax.plot(diameters, distributions[:, i], color=colors[i],
                label=f"#{i + 1} {header['Start Time'][i]}", lw=1.2)

    for size in PS_NOMINAL_SIZES_NM:
        ax.axvline(size, color="gray", ls="--", lw=0.8, zorder=0)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Diameter (nm)")
    ax.set_ylabel("dW/dlogDp (number)")
    ax.set_title("DMA distributions, all scans (dashed lines: PS nominal 20/30/40/50nm)")
    ax.legend(fontsize=6, ncol=2)
    plt.tight_layout()
    out_path = os.path.join(FIGURES_DIR, "dma_distributions_overlay.png")
    plt.savefig(out_path, dpi=120)
    plt.close(fig)
    print(f"\nSaved {out_path}")


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
        ax.set_title(f"#{i + 1} {header['Start Time'][i]}\n"
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
    plot_overlay(diameters, distributions, header, n_scans)
    plot_grid(diameters, distributions, header, footer, n_scans)
