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


ARCSEC_PER_RAD = 206264.806


def add_pc_axis(ax, distance_kpc: float, per_unit_arcsec: float = 1.0, label: str | None = None):
    """Add a secondary top x axis in pc to an axis whose x is angular.

    ``per_unit_arcsec`` is the number of arcsec in one unit of the parent axis
    (1 for arcsec, 60 for arcmin, 3600 for degrees).
    """
    k = per_unit_arcsec * distance_kpc * 1e3 / ARCSEC_PER_RAD
    sec = ax.secondary_xaxis("top", functions=(lambda a: a * k, lambda p: p / k))
    sec.set_xlabel(label if label is not None else f"r  [pc]  (D = {distance_kpc:.2f} kpc)",
                   color=INK_SECONDARY, fontsize=9)
    sec.tick_params(colors=INK_SECONDARY, labelsize=8)
    return sec


def add_arcsec_axis(ax, distance_kpc: float, per_unit_pc: float = 1.0, label: str | None = None,
                    unit: str = "arcsec"):
    """Add a secondary top x axis in arcsec (or arcmin) to an axis whose x is in pc."""
    per_arcsec = 60.0 if unit == "arcmin" else 1.0
    k = per_unit_pc * ARCSEC_PER_RAD / (distance_kpc * 1e3 * per_arcsec)
    sec = ax.secondary_xaxis("top", functions=(lambda p: p * k, lambda a: a / k))
    sec.set_xlabel(label if label is not None else f"r  [{unit}]  (D = {distance_kpc:.2f} kpc)",
                   color=INK_SECONDARY, fontsize=9)
    sec.tick_params(colors=INK_SECONDARY, labelsize=8)
    return sec


# --------------------------------------------------------------------- datasets ---
# One symbol and one colour per dataset, used identically in every figure. The
# component is carried by the marker fill, never by a different symbol, so a radial
# and a tangential measurement of the same stars are visibly the same measurement:
#   full  = the scalar quantity (combined PM dispersion, or a line-of-sight one)
#   left  = radial component
#   right = tangential component
DATASET_STYLE = {
    "hst":       {"color": SERIES[1], "marker": "o", "label": "HST (oMEGACat)"},
    "muse":      {"color": SERIES_EXTRA, "marker": "s", "label": "MUSE, line of sight"},
    "gaia_dr2":  {"color": INK, "marker": "P", "label": "Gaia DR2 (published)"},
    "gaia_edr3": {"color": SERIES[0], "marker": "D", "label": "Gaia EDR3 (our measurement)"},
    "pristine":  {"color": SERIES[2], "marker": "v", "label": "Pristine periphery"},
    "spectro":   {"color": INK_SECONDARY, "marker": "*", "label": "periphery spectroscopy"},
}

#: dataset key of the likelihood -> (style key, component)
DATASET_KEY_MAP = {
    "hst_pm_radial": ("hst", "radial"),
    "hst_pm_tangential": ("hst", "tangential"),
    "hst_pm_combined": ("hst", "combined"),
    "muse_los_dispersion": ("muse", "combined"),
    "gaia_dr2_pm": ("gaia_dr2", "combined"),
    "gaia_edr3_pm": ("gaia_edr3", "combined"),
    "gaia_edr3_ours": ("gaia_edr3", "combined"),
    "gaia_edr3_ours_radial": ("gaia_edr3", "radial"),
    "gaia_edr3_ours_tangential": ("gaia_edr3", "tangential"),
    "hst_pm_radial_ours": ("hst", "radial"),
    "hst_pm_tangential_ours": ("hst", "tangential"),
}

_FILL = {"combined": "full", "radial": "left", "tangential": "right"}


def dataset_style(key: str, component: str = "combined", size: float = 6.5,
                  fitted: bool = True) -> dict:
    """Marker keyword arguments for one dataset and component.

    ``fitted=False`` only lightens the mark, for data shown alongside the primary ones. The
    label is unchanged: which datasets the likelihood receives is stated once, in the
    coverage panel, not repeated on every legend entry.
    """
    if key in DATASET_KEY_MAP:
        key, component = DATASET_KEY_MAP[key]
    spec = DATASET_STYLE[key]
    kw = {"color": spec["color"], "marker": spec["marker"], "markersize": size,
          "fillstyle": _FILL[component], "markerfacecoloralt": SURFACE,
          "markeredgewidth": 1.2, "alpha": 1.0 if fitted else 0.55}
    if not fitted:
        kw["markeredgewidth"] = 1.0
    return kw


def dataset_label(key: str, component: str = "combined", fitted: bool = True) -> str:
    """Legend text matching :func:`dataset_style`."""
    if key in DATASET_KEY_MAP:
        key, component = DATASET_KEY_MAP[key]
    text = DATASET_STYLE[key]["label"]
    if component != "combined":
        text += ", %s" % component
    return text
