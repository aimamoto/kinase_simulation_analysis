#!/usr/bin/env python
"""Figure 6 — Model 2 CSK-SRC structural compound figure. Paired rows are h-concatenated
(controlled gap) so left/right panels sit tight. Auto-crops struct/*.png.

Representative models (see the matching struct/*.cxc for the exact open commands); all from
condition a-csk-wtcat-holo_b-src-wtcat-py159-holo, coordinates bit-identical between the
260606 and 260718 datasets:
  D  seed-1759538_sample-2  SRC State 8 (SB 2.98 A) x CSK State 5 (SB 17.67 A)  [v7r3]
  E  seed-1593656_sample-0  CSK State 3, fully active (SB 2.81 A)               [v7r3]
  E  seed-10345_sample-0    CSK State 9, incomplete   (SB 16.38 A)              [v7r3]
Quoted K-E distances are AlloQuant SB_Dist = min(NZ-OE1, NZ-OE2) from the v7r3 Phase7 tables."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from PIL import Image
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)
S=os.environ.get("ALLOQUANT_STRUCT", "struct/")  # ChimeraX-rendered structural panels
OUTDIR = os.environ.get("ALLOQUANT_FIGDIR", ".")  # rendered-figure output dir

def crop(fn,pad=10,thr=246):
    im=np.array(Image.open(fn).convert("RGB")); mask=(im<thr).any(2); ys,xs=np.where(mask)
    y0,y1,x0,x1=ys.min(),ys.max(),xs.min(),xs.max()
    return im[max(0,y0-pad):y1+pad,max(0,x0-pad):x1+pad]

def hconcat(fns,gap=16,H=760):
    ims=[]
    for fn in fns:
        im=Image.fromarray(crop(S+fn)); w=int(round(im.width*H/im.height))
        ims.append(np.array(im.resize((w,H))))
    widths=[im.shape[1] for im in ims]; tot=sum(widths)+gap*(len(ims)-1)
    canvas=np.full((H,tot,3),255,np.uint8); x=0; lefts=[]; centers=[]
    for im in ims:
        canvas[:,x:x+im.shape[1]]=im; lefts.append(x/tot); centers.append((x+im.shape[1]/2)/tot); x+=im.shape[1]+gap
    return canvas,lefts,centers

def letter(ax,l,x,y=1.0):
    ax.text(x,y,l,transform=ax.transAxes,fontsize=P.TS["panel"],fontweight="bold",va="top",ha="right")
def ttl(ax,t,x,y=1.008):
    ax.text(x,y,t,transform=ax.transAxes,fontsize=P.TS["small"],fontweight="bold",va="bottom",ha="center")
def cap(ax,t,x,color=P.INK):
    ax.text(x,-0.015,t,transform=ax.transAxes,fontsize=P.TS["small"],va="top",ha="center",color=color)

GRN="#1E8449"; RED="#B03030"
fig=plt.figure(figsize=(7.5,9.4))
gs=fig.add_gridspec(4,1,height_ratios=[1.05,0.92,0.92,0.92],hspace=0.22,
                    left=0.055,right=0.985,top=0.955,bottom=0.085)

# A — overview (two views, 180° apart)
axA=fig.add_subplot(gs[0]); img,lf,ct=hconcat(["m2_ov_v3_SRCleft1.png","m2_ov_A2.png"]); axA.imshow(img); axA.axis("off")
letter(axA,"A",lf[0]); ttl(axA,"The CSK–SRC kinase–kinase complex (two views, 180° apart)",0.5)
cap(axA,"interface / pY419 face",ct[0]); cap(axA,"opposite face (C-tail, Y530)",ct[1])

# B|C — interface + spine
axBC=fig.add_subplot(gs[1]); img,lf,ct=hconcat(["m2_panelB_a.png","m2_panelC3_b.png"]); axBC.imshow(img); axBC.axis("off")
letter(axBC,"B",lf[0]); letter(axBC,"C",lf[1])
ttl(axBC,"Basic–acidic docking interface",ct[0]); ttl(axBC,"The coupled R-spine",ct[1])

# D — encounter (SRC active | CSK poised)
axD=fig.add_subplot(gs[2]); img,lf,ct=hconcat(["m2Dv2_SRCactive.png","m2Dv2_CSKpoised.png"]); axD.imshow(img); axD.axis("off")
letter(axD,"D",lf[0]); ttl(axD,"Single-structure handshake: active SRC docked on poised CSK",0.5)
cap(axD,"SRC active: K–E 3.0 Å",ct[0],GRN); cap(axD,"CSK poised: K–E 17.7 Å",ct[1],RED)

# E — CSK activation (State 8 | State 9)
axE=fig.add_subplot(gs[3]); img,lf,ct=hconcat(["m2Ev2_state8.png","m2Ev2_state9.png"]); axE.imshow(img); axE.axis("off")
letter(axE,"E",lf[0]); ttl(axE,"CSK partial activation (≈28% reach the active State 3)",0.5)
cap(axE,"State 3 active: K–E 2.8 Å",ct[0],GRN); cap(axE,"State 9: K–E 16.4 Å",ct[1],RED)

# legend — shortened terms
leg=[("CSK","#cfd3d6"),("SRC","#c9b89a"),("interface","#e67e22"),("αC-helix","#3498DB"),
     ("β3(K)–αC(E)","#27AE60"),("DFG","#D35400"),("HRD-Asp","#E74C3C"),("A-loop","#FF7F50"),
     ("R-spine","#8E44AD"),("ATP","#DBA500"),("pY419 (SRC)","#F1C40F"),("Y530 (SRC)","#E0119B")]
handles=[Patch(facecolor=c,edgecolor=P.INK,linewidth=0.4,label=l) for l,c in leg]
fig.legend(handles=handles,loc="lower center",ncol=6,frameon=False,fontsize=P.TS["small"],
           handlelength=1.1,columnspacing=1.3,handletextpad=0.4,bbox_to_anchor=(0.5,0.005))

for ext in ("png","pdf"):
    fig.savefig(f"{OUTDIR}/Figure6_CSK_SRC_structure.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
print("saved Figure6_CSK_SRC_structure.{png,pdf}")
