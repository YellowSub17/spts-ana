"""
Single source of truth for the combined focus filter used across all the
focus/shape/sizing scripts:

    focused = flags & (solidity >= SOLIDITY_THRESHOLD)

Import SOLIDITY_THRESHOLD from here rather than redefining it locally -- every
script in this directory that filters or classifies particles by solidity should
import it from this module. (Previously each script had its own hardcoded copy,
which meant changing the threshold in one script silently had no effect on the
others -- see ANALYSIS_SUMMARY_p1.md section 6.)
"""

SOLIDITY_THRESHOLD = 0.80
