#!/usr/bin/env python
"""Fig S4 (Model 2 dynamics) — Panel A: rigid loader / breathing receiver (RMSF);
Panel B: SRC C-tail (Y530) tucked, CSK pocket accessible. Combines former S4+S5."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from PIL import Image
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)
S=os.environ.get("ALLOQUANT_STRUCT", "struct/")  # ChimeraX-rendered structural panels
OUTDIR = os.environ.get("ALLOQUANT_FIGDIR", ".")  # rendered-figure output dir
def crop(fn,pad=8,thr=246):
    im=np.array(Image.open(fn).convert("RGB")); mask=(im<thr).any(2); ys,xs=np.where(mask)
    y0,y1,x0,x1=ys.min(),ys.max(),xs.min(),xs.max()
    return im[max(0,y0-pad):y1+pad,max(0,x0-pad):x1+pad]

fig=plt.figure(figsize=(7.6,4.3))
# --- Panel A: RMSF ---
axA=fig.add_axes([0.015,0.06,0.44,0.86]); axA.imshow(crop(S+"m2_rmsf_a.png")); axA.axis("off")
axA.text(-0.01,1.0,"A",transform=axA.transAxes,fontsize=P.TS["panel"],fontweight="bold",va="top",ha="right")
axA.text(0.5,1.01,"Rigid SRC loader, breathing CSK receiver",transform=axA.transAxes,fontsize=P.TS["small"],fontweight="bold",va="bottom",ha="center")
axA.text(0.60,0.90,"SRC:\nrigid",transform=axA.transAxes,fontsize=P.TS["small"],ha="left",va="center",color="#2166ac")
axA.text(0.02,0.14,"CSK:\nmobile A-loop",transform=axA.transAxes,fontsize=P.TS["small"],ha="left",va="center",color="#B03030")
cmap=LinearSegmentedColormap.from_list("rmsf",["#2166ac","#f7f7f7","#b2182b"])
cax=fig.add_axes([0.455,0.24,0.016,0.50]); cb=plt.colorbar(ScalarMappable(Normalize(0,4),cmap),cax=cax)
cb.set_label("Cα RMSF (Å)",fontsize=P.TS["small"]); cb.ax.tick_params(labelsize=P.TS["small"])
# --- Panel B: C-tail ---
# NOTE (2026-07-19): figS4_model2_dynamics.pdf was MANUALLY adjusted by the user (label
# positions) and is the MASTER. Re-running this script OVERWRITES that manual layout.
# The label coords below approximate the user's manual placement as a fallback.
axB=fig.add_axes([0.525,0.06,0.46,0.86]); axB.imshow(crop(S+"m2_ctail_v2_a.png")); axB.axis("off")
axB.text(-0.01,1.0,"B",transform=axB.transAxes,fontsize=P.TS["panel"],fontweight="bold",va="top",ha="right")
axB.text(0.5,1.01,"Substrate tail (Y530) tucked, ≈32 Å from the CSK pocket",transform=axB.transAxes,fontsize=P.TS["small"],fontweight="bold",va="bottom",ha="center")
axB.annotate("SRC C-tail (Y530 ●)",xy=(0.56,0.62),xytext=(0.52,0.94),
             xycoords="axes fraction",textcoords="axes fraction",
             fontsize=P.TS["small"],ha="center",va="center",color="#2E4053",
             arrowprops=dict(arrowstyle="->",color="#2E4053",lw=0.9))
axB.annotate("CSK catalytic\npocket (ATP)",xy=(0.31,0.55),xytext=(0.03,0.86),
             xycoords="axes fraction",textcoords="axes fraction",
             fontsize=P.TS["small"],ha="left",va="center",color=P.INK,
             arrowprops=dict(arrowstyle="->",color=P.INK,lw=0.9))
for ext in ("png","pdf"):
    fig.savefig(f"{OUTDIR}/figS4_model2_dynamics.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
print("saved figS4_model2_dynamics.{png,pdf}")
