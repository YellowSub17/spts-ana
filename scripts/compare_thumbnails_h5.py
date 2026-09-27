"""
Check that two thumbnails h5 files hold the same data, group by group.

    python compare_thumbnails_h5.py data/thumbnails.h5 data/thumbnails_t20.h5
    python compare_thumbnails_h5.py data/thumbnails.h5 data/thumbnails_t10.h5 --groups groel ferri
"""
import argparse

import h5py
import numpy as np


parser = argparse.ArgumentParser(description="Compare two thumbnails h5 files")
parser.add_argument("a_path")
parser.add_argument("b_path")
parser.add_argument("--groups", type=str, nargs='+', default=None, help="only compare these groups (default: all)")
args = parser.parse_args()
a_path, b_path = args.a_path, args.b_path
ok = True
with h5py.File(a_path, 'r') as a, h5py.File(b_path, 'r') as b:
    for group in args.groups or sorted(set(a) | set(b)):
        if group not in a or group not in b:
            print(f'{group}: only in {a_path if group in a else b_path}')
            ok = False
            continue
        for key in sorted(set(a[group]) | set(b[group])):
            if key not in a[group] or key not in b[group]:
                print(f'{group}/{key}: only in {a_path if key in a[group] else b_path}')
                ok = False
                continue
            da, db = a[group][key][:], b[group][key][:]
            if da.shape != db.shape:
                print(f'{group}/{key}: shape {da.shape} vs {db.shape}')
                ok = False
            elif not np.array_equal(da, db, equal_nan=da.dtype.kind == 'f'):
                print(f'{group}/{key}: values differ ({np.sum(da != db)} elements)')
                ok = False
print('IDENTICAL' if ok else 'DIFFERENT')
