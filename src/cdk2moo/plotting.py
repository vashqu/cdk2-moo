"""Shared figure style and small plotting helpers.

Will live here:

- one matplotlib/seaborn style applied across every figure
- a fixed colour assignment per molecule set (GA arms, known actives, random
  ChEMBL, decoys) so colours mean the same thing in every panel
- save helper writing to ``config.FIGURES_DIR`` at a fixed dpi

Axis limits are never trimmed to flatter a trend. If a figure looks bad, the
figure is telling you something.
"""
