"""AlloQuant manuscript figure palette — LOCKED from Fig1.odg (Flat-UI system).

Colours extracted from the Fig 1 embedded SVG (dominant hues) + Draw text boxes.
Usage rules (see PALETTE.md): encode MAGNITUDE (MAC) with the sequential blue
RIGIDITY ramp; reserve RED = tense/rigid pole and GREEN = active/output as
semantic accents that are ALWAYS direct-labelled (never rely on red-vs-green
alone — it fails CVD). Keep categorical sets to CVD-distinguishable hues + labels.
"""
from matplotlib.colors import LinearSegmentedColormap

# --- core ink / neutrals ---
INK      = "#2C3E50"   # primary text, axes, borders
INK2     = "#34495E"   # secondary ink
MUTED    = "#7F8C8D"   # tertiary / annotation grey
GRID     = "#BDC3C7"   # gridlines
FILL_LT  = "#ECF0F1"   # light panel fill
WHITE    = "#FFFFFF"

# --- semantic hues (Flat-UI, matching Fig 1) ---
BLUE, BLUE_D     = "#3498DB", "#2980B9"   # structure / N-lobe / primary data
RED,  RED_D      = "#E74C3C", "#C0392B"   # tense / rigid pole / C-lobe features
GREEN, GREEN_D   = "#27AE60", "#1E8449"   # active / output / positive
AMBER            = "#F39C12"              # nucleotide / intermediate
PURPLE           = "#8E44AD"              # special landmark / 4th category

# --- sequential rigidity ramp: fluid (light) -> rigid (navy). MAGNITUDE encoding. ---
RIGIDITY_STOPS = ["#EAF2FB", "#AED4F0", "#5DAEE3", "#2E86C1", "#21618C", "#1B2F45"]
RIGIDITY = LinearSegmentedColormap.from_list("rigidity", RIGIDITY_STOPS)

# 3-series categorical for canonical features (CVD-distinguishable; always labelled)
FEATURE = {"αC-In": BLUE, "DFG-in": AMBER, "BLAminus": PURPLE}

FONT_STACK = ["Liberation Sans", "Arial", "Noto Sans", "DejaVu Sans"]

# Main-figure type scale (pt) — nothing below ~7.5 at final print size
TS = dict(base=9.5, title=10.5, axlabel=10.5, tick=8.5, legend=8.5,
          panel=14, anno=7.5, anno_big=10, value=8.0, small=7.5)

def apply_rc(mpl):
    mpl.rcParams.update({
        "font.family": FONT_STACK, "font.size": TS["base"],
        "axes.titlesize": TS["title"], "axes.labelsize": TS["axlabel"],
        "xtick.labelsize": TS["tick"], "ytick.labelsize": TS["tick"], "legend.fontsize": TS["legend"],
        "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK,
        "axes.linewidth": 0.9, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
        "figure.dpi": 150, "savefig.dpi": 300, "svg.fonttype": "none",
    })
