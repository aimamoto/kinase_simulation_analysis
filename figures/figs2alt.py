#!/usr/bin/env python3
"""Figure S2 alt — monomeric CSK apo vs holo side-by-side composite.

ChimeraX renders in struct/figs2alt_apo.png and struct/figs2alt_holo.png.
Run those via ChimeraX --nogui before this script.
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import numpy as np
import os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import palette as P

STRUCT = os.environ.get("ALLOQUANT_STRUCT", "struct/")  # ChimeraX-rendered structural panels
OUTDIR = os.environ.get("ALLOQUANT_FIGDIR", ".")  # rendered-figure output dir
APO_RENDER  = os.path.join(STRUCT, "figs2alt_apo.png")
HOLO_RENDER = os.path.join(STRUCT, "figs2alt_holo.png")
BASE        = os.path.join(OUTDIR, "FigS2alt")

LABEL_APO  = "Apo (−ATP)"
LABEL_HOLO = "Holo (+ATP)"

def autocrop(img, threshold=0.97, margin=6):
    """Trim near-white borders from a float [0,1] RGB(A) image."""
    gray = img[:, :, :3].min(axis=2)
    mask = gray < threshold
    rows = np.where(mask.any(axis=1))[0]
    cols = np.where(mask.any(axis=0))[0]
    if not len(rows) or not len(cols):
        return img
    r0 = max(rows[0]  - margin, 0);  r1 = min(rows[-1]  + margin + 1, img.shape[0])
    c0 = max(cols[0]  - margin, 0);  c1 = min(cols[-1]  + margin + 1, img.shape[1])
    return img[r0:r1, c0:c1]

apo_img  = autocrop(mpimg.imread(APO_RENDER))
holo_img = autocrop(mpimg.imread(HOLO_RENDER))

fig = plt.figure(figsize=(6.5, 4.28), facecolor="white", dpi=220)
gs  = GridSpec(2, 2, figure=fig,
               left=0.02, right=0.98, top=0.95, bottom=0.10,
               hspace=0.18, wspace=0.04,
               height_ratios=[3.2, 1.0])

ax_apo  = fig.add_subplot(gs[0, 0]); ax_apo.set_axis_off()
ax_holo = fig.add_subplot(gs[0, 1]); ax_holo.set_axis_off()
ax_apo.imshow(apo_img)
ax_holo.imshow(holo_img)
for ax, label in [(ax_apo, LABEL_APO), (ax_holo, LABEL_HOLO)]:
    ax.text(0.5, 1.0, label, transform=ax.transAxes,
            ha="center", va="bottom", fontsize=12, fontweight="bold", color=P.INK)

ax_apo.text(0.5, -0.02,
    "αC-In, active",
    transform=ax_apo.transAxes, ha="center", va="top",
    fontsize=9.5, fontstyle="italic", color=P.INK)
ax_holo.text(0.5, -0.02,
    "αC-Out, inactive",
    transform=ax_holo.transAxes, ha="center", va="top",
    fontsize=9.5, fontstyle="italic", color=P.INK)

ax_tab = fig.add_subplot(gs[1, :])
ax_tab.set_axis_off()

col_labels = ["", LABEL_APO, LABEL_HOLO]
rows = [
    ["αC-In",    "225 / 225  (100%)",  "1 / 225  (0.4%)"],
    ["BLAminus", "216 / 225  (96%)",   "60 / 225  (27%)"],
]

tbl = ax_tab.table(
    cellText=rows,
    colLabels=col_labels,
    loc="center",
    cellLoc="center",
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(9.5)
tbl.scale(1, 1.45)

for col in range(3):
    cell = tbl[0, col]
    cell.set_facecolor(P.FILL_LT)
    cell.set_text_props(fontweight="bold", color=P.INK)

tbl[1, 2].set_text_props(color=P.RED, fontweight="bold")

for row in (1, 2):
    tbl[row, 0].set_text_props(fontweight="bold")

fig.text(0.5, 0.03,
    "Fisher's exact $p$ (αC-In) = 4×10⁻¹³²"
    "  |  $n$ = 225 per condition, 45 independent seeds",
    ha="center", va="bottom", fontsize=8.5, fontstyle="italic", color=P.INK2)

for ext, kw in [("png", {"dpi": 220}), ("pdf", {})]:
    out = f"{BASE}.{ext}"
    fig.savefig(out, facecolor="white", bbox_inches="tight", **kw)
    print("Wrote", out)
