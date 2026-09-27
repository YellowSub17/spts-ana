"""
Check that two thumbnails h5 files hold the same data, group by group.

    python compare_thumbnails_h5.py data/thumbnails.h5 data/thumbnails_t20.h5
"""
import sys

import h5py
import numpy as np


a_path, b_path = sys.argv[1], sys.argv[2]
ok = True
with h5py.File(a_path, 'r') as a, h5py.File(b_path, 'r') as b:
    for group in sorted(set(a) | set(b)):
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
