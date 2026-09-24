"""
Recreates the particle-grouping / focus-filtering work from the thumbnails.h5 analysis.

Background
----------
thumbnails.h5 already carries a boolean `flags` dataset per group, produced by the
old two-radius method (src/spts_ana/comborun/filtering.py: filter_focused): a particle is
"focused" if the summed intensity in a 7px window is >= 90% of the summed intensity
in a 25px window (comparing amplitudes, i.e. intensity**(1/6)).

Checking that method against the thumbnails directly showed two things:
  1. The focused/unfocused populations are smooth, heavily overlapping distributions
     with no natural gap near the 0.9 threshold -- especially for the larger
     polystyrene beads (ps30-ps50nm), where a bigger *in-focus* spot is hard to tell
     apart from a smaller *defocused* one using intensity ratio alone.
  2. Even within flags==True, a meaningful fraction of thumbnails are visibly
     non-round: comet tails, double lobes, bean/kidney shapes. These pass the
     intensity-ratio test but are not a single compact spot.

This script adds two per-particle shape metrics, computed directly from the
thumbnail pixels, to catch (2):

  solidity  = (thresholded blob area) / (area of that blob's convex hull)
              Low solidity (<~0.8) means the bright region is non-convex --
              exactly the bean/double-lobe shapes we want to filter out.

  gauss_r2  = R^2 of a 2D elliptical Gaussian fit to the blob.
              Lower for multi-lobed / non-Gaussian spots. Less reliable for the
              larger PS beads (real diffraction pattern is closer to an Airy
              profile than a Gaussian even when in focus), so treat as a
              secondary check, not the primary filter.

Combined focus criterion used downstream:
    focused_v2 = flags & (solidity >= SOLIDITY_THRESHOLD)

Output
------
Writes data/focus_shape_metrics.h5, with one group per particle type (same names,
same row order as thumbnails.h5), each containing:
    solidity, gauss_r2, gauss_sx, gauss_sy, blob_area

Run this after thumbnails.h5 exists. It does not modify thumbnails.h5.
"""

import time

import h5py
import numpy as np
from scipy import ndimage
from scipy.optimize import curve_fit
from scipy.spatial import ConvexHull, QhullError

from filter_config import SOLIDITY_THRESHOLD

THUMBNAILS_H5 = "/Users/pat/Documents/work/spts-ana/data/thumbnails.h5"
OUTPUT_H5 = "/Users/pat/Documents/work/spts-ana/data/focus_shape_metrics.h5"

# SOLIDITY_THRESHOLD is imported from filter_config.py -- the single source of
# truth for every script that filters/classifies by solidity. This script itself
# doesn't apply the threshold (it just computes solidity for every particle); it's
# imported here only so the summary line below reports the value actually in use
# elsewhere, instead of a second, independent copy that could drift out of sync.

THUMB_SIZE = 60
_yy, _xx = np.mgrid[0:THUMB_SIZE, 0:THUMB_SIZE]


def gauss2d(coords, amp, x0, y0, sx, sy, theta, bg):
    x, y = coords
    a = (np.cos(theta) ** 2) / (2 * sx**2) + (np.sin(theta) ** 2) / (2 * sy**2)
    b = -(np.sin(2 * theta)) / (4 * sx**2) + (np.sin(2 * theta)) / (4 * sy**2)
    c = (np.sin(theta) ** 2) / (2 * sx**2) + (np.cos(theta) ** 2) / (2 * sy**2)
    return bg + amp * np.exp(
        -(a * (x - x0) ** 2 + 2 * b * (x - x0) * (y - y0) + c * (y - y0) ** 2)
    )


def analyze(thumb):
    """Compute (solidity, gauss_r2, gauss_sx, gauss_sy, blob_area) for one thumbnail."""
    # background from the corners, away from the particle in the center
    corner = np.concatenate(
        [
            thumb[:6, :6].ravel(),
            thumb[-6:, -6:].ravel(),
            thumb[:6, -6:].ravel(),
            thumb[-6:, :6].ravel(),
        ]
    )
    bg0 = np.median(corner)
    img = thumb - bg0

    peak_val = img.max()
    if peak_val <= 0:
        return np.nan, np.nan, np.nan, np.nan, 0

    # threshold at 25% of peak, keep the connected component containing the peak
    mask = img > 0.25 * peak_val
    lbl, n_components = ndimage.label(mask)
    if n_components == 0:
        return np.nan, np.nan, np.nan, np.nan, 0

    cy, cx = np.unravel_index(np.argmax(img), img.shape)
    comp_id = lbl[cy, cx]
    if comp_id == 0:
        # peak pixel didn't land in a labeled component (shouldn't normally happen);
        # fall back to the largest component
        sizes = ndimage.sum(mask, lbl, range(1, n_components + 1))
        comp_id = np.argmax(sizes) + 1
    comp_mask = lbl == comp_id
    area = int(comp_mask.sum())
    if area < 4:
        return np.nan, np.nan, np.nan, np.nan, area

    # solidity: blob area vs. its convex hull area.
    # Use pixel *corners* (not centers) for the hull so it properly bounds the
    # pixel area -- using centers alone underestimates the hull area and gives
    # solidity > 1.
    ys_pix, xs_pix = np.where(comp_mask)
    corner_offsets = np.array([[-0.5, -0.5], [0.5, -0.5], [-0.5, 0.5], [0.5, 0.5]])
    pts = (
        np.column_stack([xs_pix, ys_pix])[:, None, :].astype(float)
        + corner_offsets[None, :, :]
    ).reshape(-1, 2)
    try:
        hull_area = ConvexHull(pts).volume  # 2D "volume" is area
    except QhullError:
        hull_area = area
    solidity = area / hull_area if hull_area > 0 else np.nan

    # 2D elliptical Gaussian fit, in a window around the blob centroid
    y0_guess, x0_guess = ys_pix.mean(), xs_pix.mean()
    r = 20
    y0i, x0i = int(round(y0_guess)), int(round(x0_guess))
    ylo, yhi = max(0, y0i - r), min(THUMB_SIZE, y0i + r)
    xlo, xhi = max(0, x0i - r), min(THUMB_SIZE, x0i + r)
    sub = img[ylo:yhi, xlo:xhi]
    sub_x, sub_y = _xx[ylo:yhi, xlo:xhi], _yy[ylo:yhi, xlo:xhi]
    p0 = [peak_val, x0_guess, y0_guess, 3.0, 3.0, 0.0, 0.0]
    try:
        popt, _ = curve_fit(
            gauss2d, (sub_x.ravel(), sub_y.ravel()), sub.ravel(), p0=p0, maxfev=1500
        )
        fit = gauss2d((sub_x, sub_y), *popt)
        ss_res = np.sum((sub - fit) ** 2)
        ss_tot = np.sum((sub - sub.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        sx, sy = abs(popt[3]), abs(popt[4])
    except Exception:
        r2, sx, sy = np.nan, np.nan, np.nan

    return solidity, r2, sx, sy, area


if __name__ == "__main__":
    fin = h5py.File(THUMBNAILS_H5, "r")
    fout = h5py.File(OUTPUT_H5, "w")

    for group in fin.keys():
        thumbnails = fin[group]["thumbnails"]
        n = thumbnails.shape[0]

        solidity = np.full(n, np.nan)
        gauss_r2 = np.full(n, np.nan)
        gauss_sx = np.full(n, np.nan)
        gauss_sy = np.full(n, np.nan)
        blob_area = np.full(n, np.nan)

        t0 = time.time()
        for i in range(n):
            solidity[i], gauss_r2[i], gauss_sx[i], gauss_sy[i], blob_area[i] = analyze(
                thumbnails[i]
            )
            if i % 500 == 0:
                print(f"{group}: {i}/{n}  ({time.time() - t0:.1f}s)", end="\r")
        print(f"{group}: {n}/{n} done in {time.time() - t0:.1f}s")

        grp = fout.create_group(group)
        grp.create_dataset("solidity", data=solidity)
        grp.create_dataset("gauss_r2", data=gauss_r2)
        grp.create_dataset("gauss_sx", data=gauss_sx)
        grp.create_dataset("gauss_sy", data=gauss_sy)
        grp.create_dataset("blob_area", data=blob_area)

    fout.close()
    print(f"All done. Wrote {OUTPUT_H5}")
    print(
        f"Combined focus filter used downstream: flags & (solidity >= {SOLIDITY_THRESHOLD})"
    )
