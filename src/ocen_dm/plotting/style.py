"""House style for static figures: palette, ink, surfaces, mark sizes.

The categorical slots and the sequential ramp are the validated reference palette
of the data-viz method (first three slots pass every colour-vision check for
scatter plots; past three, fold to "other"). Colours identify entities, never
rank; text always wears ink, never a series colour.
"""

from __future__ import annotations

from matplotlib import colors as mcolors
from matplotlib import pyplot as plt

#: categorical slots, fixed order -- never cycled
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")
#: what each slot means across every figure in this project
COLOR_KUZMA2025 = SERIES[0]     # Gaia + Pristine periphery members
COLOR_KUZMA2026 = SERIES[1]     # spectroscopic members
COLOR_FIMBULTHUL = SERIES[2]    # Fimbulthul stream candidates
COLOR_FIELD = "#b9b8b3"         # non-members / background: neutral, recessive

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e6e5e1"

#: single-hue sequential ramp (blue, light -> dark) for densities
_SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
        "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]
SEQUENTIAL = mcolors.LinearSegmentedColormap.from_list("ocen_blue", _SEQ)

#: diverging: blue <-> red through a neutral grey midpoint, for signed velocities
DIVERGING = mcolors.LinearSegmentedColormap.from_list(
    "ocen_div", ["#0d366b", "#2a78d6", "#86b6ef", "#f0efec", "#f0a3a2", "#e34948", "#8f1d1c"])
#: series-4 slot used only where a fourth line is unavoidable
SERIES_EXTRA = "#eda100"

MARKER_PT = 5.0        # >= 8 px at 150 dpi
LINE_PT = 1.5          # ~2 px
ERR_PT = 1.0


def apply() -> None:
    """Set the matplotlib defaults used by every figure."""
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "text.color": INK,
        "axes.labelcolor": INK,
        "axes.edgecolor": INK_SECONDARY,
        "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "medium",
        "legend.frameon": False,
        "legend.fontsize": 9,
        "lines.linewidth": LINE_PT,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })
