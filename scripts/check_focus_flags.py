"""
Sanity check that each particle's 'flags' (focus filter) belongs to that particle.

The focus filter compares a peak's intensity in a small vs a large integration
window (d7 vs d25 spts analyses). A proxy for the same ratio can be computed
straight from each row's thumbnail: intensity within r=3.5 px of the centre over
intensity within r=12.5 px. If flags are attached to the right particles, the
proxy separates flagged from unflagged almost perfectly (AUC ~0.99 on the
original thumbnails.h5). An AUC near 0.5 means the flags are effectively
scrambled relative to the particles.

    python check_focus_flags.py data/thumbnails_t20.h5 data/thumbnails_t10.h5
"""
import sys

import h5py
import numpy as np


R_SMALL, R_BIG = 3.5, 12.5


def auc(score, label):
    # Mann-Whitney AUC: P(score of a flagged particle > score of an unflagged one)
    order = np.argsort(score)
    ranks = np.empty(len(score)); ranks[order] = np.arange(1, len(score) + 1)
    n1 = label.sum(); n0 = len(label) - n1
    return (ranks[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


for path in sys.argv[1:]:
    print(path)
    with h5py.File(path, 'r') as f:
        for group in f:
            thumbnails = f[group]['thumbnails'][:]
            flags = f[group]['flags'][:]
            yy, xx = np.mgrid[:thumbnails.shape[1], :thumbnails.shape[2]]
            r = np.hypot(yy - (thumbnails.shape[1] - 1) / 2, xx - (thumbnails.shape[2] - 1) / 2)

            # peaks too close to the detector edge get all-zero thumbnails
            valid = thumbnails.reshape(len(thumbnails), -1).any(1)
            small = thumbnails[valid][:, r <= R_SMALL].sum(1)
            big = thumbnails[valid][:, r <= R_BIG].sum(1)
            proxy = np.divide(small, big, out=np.zeros_like(small), where=big > 0)

            score = auc(proxy, flags[valid])
            status = 'OK' if score > 0.95 else 'SUSPECT'
            print(f'  {group:8s} n={valid.sum():5d} focused={flags[valid].sum():5d} AUC={score:.3f}  {status}')
