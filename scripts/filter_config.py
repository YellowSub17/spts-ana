"""
Single source of truth for the filters used across the focus/shape/sizing
scripts. Import from here rather than redefining locally -- every script that
applies these filters should import them from this module. (Previously
SOLIDITY_THRESHOLD had independent hardcoded copies in several scripts, which
meant changing it in one place silently had no effect on the others -- see
ANALYSIS_SUMMARY_p1.md section 6.)

    focused = flags & (solidity >= SOLIDITY_THRESHOLD)

SOLIDITY_THRESHOLD: shape-based focus filter (see compute_focus_shape_metrics.py).

Y_ILLUMINATION_MIN/MAX: crop on particle y-position (thumbnails.h5 'ys' field),
used in the intensity/sizing scripts (ps_calibration_sizing.py and friends).
plot_intensity_vs_position.py found a clear y-dependent illumination gradient --
median normalized intensity is flat for y in ~[30,167] and falls off steadily
above that (down to ~0.75x by y~230) -- that has nothing to do with particle
size but biases each group's median intensity depending on how its particles
happen to be distributed in y. Restricting to this y range tightened the PS
calibration fit from R^2=0.989 to R^2=0.9998. The x-profile showed no comparable
trend, so no x crop is applied.
"""

SOLIDITY_THRESHOLD = 0.80

Y_ILLUMINATION_MIN = 30
Y_ILLUMINATION_MAX = 167
