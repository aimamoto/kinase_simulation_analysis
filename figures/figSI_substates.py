#!/usr/bin/env python3
"""SI figure: pT161 (active) vs triple-P (inhibited) each coexist as a rigid (State 2)
and floppy (State 1) conformer; inhibition repopulates toward floppy (Fisher p=0.022)."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P

D = os.environ.get("ALLOQUANT_STRUCT", "struct/")  # ChimeraX-rendered structural panels
# grid[row][col]; rows = State1(floppy), State2(rigid); cols = pT161 active, triple-P inhibited
cells = [
    [("fig3ens_pt161_S1_w.png", "49%"), ("fig3ens_triplep_S1_w.png", "66%")],
    [("fig3ens_pt161_S2_w.png", "51%"), ("fig3ens_triplep_S2_w.png", "34%")],
]
col_titles = ["pT161  (active)", "triple-P  (inhibited)"]
row_titles = ["State 1 — floppy\nintrinsic MAC 0.14", "State 2 — rigid\nintrinsic MAC 0.25"]

fig, axes = plt.subplots(2, 2, figsize=(8.2, 9.2), facecolor="white")
for i in range(2):
    for j in range(2):
        ax = axes[i][j]; fn, pct = cells[i][j]
        ax.imshow(mpimg.imread(D + fn)); ax.set_axis_off()
        # population badge
        ax.text(0.035, 0.955, pct, transform=ax.transAxes, color="white",
                fontsize=15, fontweight="bold", va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.28", fc=P.INK, ec=P.GRID, lw=0.8, alpha=0.85))
        if i == 0:
            ax.set_title(col_titles[j], color=P.INK, fontsize=13, fontweight="bold", pad=10)
        if j == 0:
            ax.text(-0.045, 0.5, row_titles[i], transform=ax.transAxes, color=P.INK,
                    fontsize=12, fontweight="bold", rotation=90, va="center", ha="center")

fig.suptitle("Active and inhibited CDK1 coexist as rigid / floppy conformers",
             color=P.INK, fontsize=15, fontweight="bold", y=0.985)
fig.text(0.5, 0.028,
         "Ensemble Cα ropes (10 models) + centroid representative. The active kinase\n"
         "sits at a ~50/50 rigid–floppy balance; inhibitory pT14/pY15 significantly\n"
         "repopulate it toward the floppy State 1 (Fisher exact p = 0.022).",
         color=P.MUTED, fontsize=11.5, ha="center", va="bottom", linespacing=1.4)
fig.subplots_adjust(left=0.075, right=0.99, top=0.88, bottom=0.135, wspace=0.02, hspace=0.06)
fig.savefig(D + "../figSI_substates.png", dpi=200, facecolor="white")
fig.savefig(D + "../figSI_substates.pdf", facecolor="white")
print("wrote figSI_substates.{png,pdf}")
