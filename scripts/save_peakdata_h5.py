"""
Build data/thumbnails_t{thresh}.h5 from the per-run spts analyses at a given
detection threshold.

Reads   data/data{run:05}_spts{d}t{thresh}/spts.cxi   for d in (25, 7)
Writes  data/thumbnails_t{thresh}.h5

The d25 analysis supplies the peaks/thumbnails; the d7 analysis is only used
for the focus filter (ComboRun.filter_focused), stored as 'flags'. Output
layout matches the original data/thumbnails.h5 (one group per sample, same
dataset names/dtypes/chunking), so downstream scripts work unchanged if pointed
at the new file. The old thumbnails.h5 mixed thresholds: the ps* groups were
analysed at threshold 20 and groel/ferri at threshold 10. So thumbnails_t20.h5
should reproduce its ps* groups and thumbnails_t10.h5 its groel/ferri groups.

Usage:
    python save_peakdata_h5.py 20
    python save_peakdata_h5.py 10 --groups ps40nm ps50nm
"""
import argparse
import os

import h5py
import hdf5plugin  # noqa: F401  (registers compression filters needed to read spts.cxi)

import comborun
import config


parser = argparse.ArgumentParser(description="Save thumbnails + peak data for a given spts threshold")
parser.add_argument("thresh", type=int, help="spts detection threshold of the analysis to read (e.g. 10 or 20)")
parser.add_argument("--groups", type=str, nargs='+', default=list(config.ps_data_ranges.keys()),
                    choices=list(config.ps_data_ranges.keys()))
args = parser.parse_args()


def generate_full_paths(filenames, d):
    # filenames are 'data01245.cxd' style, as in config.ps_data_ranges
    return [f'{config.DATA_DIR}/{fname[:-4]}_spts{d}t{args.thresh}/spts.cxi' for fname in filenames]


# Fail up front rather than partway through a long run
missing = [p for group in args.groups for d in (25, 7)
           for p in generate_full_paths(config.ps_data_ranges[group], d) if not os.path.exists(p)]
if missing:
    raise SystemExit("Missing spts.cxi files:\n  " + "\n  ".join(missing))


out_file_path = f'{config.DATA_DIR}/thumbnails_t{args.thresh}.h5'

with h5py.File(out_file_path, "w") as f:
    for group in args.groups:
        print(group)
        cbr1 = comborun.ComboRun(generate_full_paths(config.ps_data_ranges[group], 25))
        cbr2 = comborun.ComboRun(generate_full_paths(config.ps_data_ranges[group], 7))
        cbr1.filter_focused(cbr2)
        thumbnails = cbr1.get_thumbnails()
        # --- H5PY WRITE START ---

        grp = f.create_group(group)
        grp.create_dataset("thumbnails", data=thumbnails, compression="gzip", chunks=(1, 60, 60))
        grp.create_dataset("flags", data=cbr1.filter[:])
        grp.create_dataset("xs", data=cbr1.peak_xs[:])
        grp.create_dataset("ys", data=cbr1.peak_ys[:])
        grp.create_dataset("is", data=cbr1.peak_is[:])
        grp.create_dataset("circum", data=cbr1.peak_circum[:])
        grp.create_dataset("max", data=cbr1.peak_max[:])
        grp.create_dataset("mean", data=cbr1.peak_mean[:])
        grp.create_dataset("median", data=cbr1.peak_median[:])
        grp.create_dataset("min", data=cbr1.peak_min[:])
        grp.create_dataset("disloc", data=cbr1.peak_disloc[:])
        grp.create_dataset("area", data=cbr1.peak_area[:])
        grp.create_dataset("eccen", data=cbr1.peak_eccen[:])

        # --- H5PY WRITE END ---

print(f"All groups successfully written to {out_file_path}")
