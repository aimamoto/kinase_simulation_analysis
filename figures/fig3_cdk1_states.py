#!/usr/bin/env python3
"""Fig 3 — CDK1 activation cycle as an ordered structural series (WHITE bg).
Five ensemble panels (10-model Ca ropes + annotated representative) in biological
order, laid out 3-over-2, with transition arrows and global-MAC labels."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
import numpy as np, sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P

D = os.environ.get("ALLOQUANT_STRUCT", "struct/")  # ChimeraX-rendered structural panels
OUTDIR = os.environ.get("ALLOQUANT_FIGDIR", ".")  # rendered-figure output dir
# global (condition-level) MAC, AlloQuant v7r3: plots_and_stats_CDK1_GMM/Phase5_Global_Network_Density.csv
# Macro-state composition is deliberately NOT annotated here -- it lives in Fig. 2C.
panels = [
    ("fig3ens_apo_w.png",     "apo",           0.3652, None),
    ("fig3ens_ccnb1_w.png",   "+ CCNB1",       0.2167, "+ CCNB1"),
    ("fig3ens_atp_w.png",     "+ CCNB1 + ATP", 0.1227, "+ ATP"),
    ("fig3ens_triplep_w.png", "triple-P",      0.1663, "+ pT14/pY15/pT161"),
    ("fig3ens_pt161_w.png",   "pT161",         0.1770, "CDC25 − pT14/pY15"),
]
# K–E (β3-Lys–αC-Glu) salt bridge: bootstrap median and 95% CI (Å), n=50/condition
# source: CDK1/af3_output/260517_CDK1-CCNB1/master_kinase_analysis_results_v7r3.csv, 2000 resamples
ke_data = {
    "apo":           (10.14,  9.59, 10.62),
    "+ CCNB1":       ( 2.74,  2.71,  2.76),
    "+ CCNB1 + ATP": ( 2.72,  2.71,  2.73),
    "triple-P":      ( 3.02,  3.01,  3.04),
    "pT161":         ( 2.99,  2.96,  3.00),
}
imgs = [mpimg.imread(D + p[0]) for p in panels]
# common crop box (union of non-WHITE regions) so scale/position match across panels
ys=[]; xs=[]
for im in imgs:
    lum = im[..., :3].min(axis=2)
    yy, xx = np.where(lum < 0.94); ys.append((yy.min(),yy.max())); xs.append((xx.min(),xx.max()))
m=30
y0=max(0, min(a for a,_ in ys)-m); y1=min(imgs[0].shape[0], max(b for _,b in ys)+m)
x0=max(0, min(a for a,_ in xs)-m); x1=min(imgs[0].shape[1], max(b for _,b in xs)+m)
imgs=[im[y0:y1, x0:x1] for im in imgs]

def maccol(v):
    # MAC on the rigidity ramp, kept in the darker band that stays visible on white
    t = 0.45 + 0.52*max(0.0, min(1.0, (v-0.12)/(0.40-0.12)))
    return P.RIGIDITY(t)

def kecol(dist):
    return P.GREEN if dist < 4.0 else P.AMBER

FA = plt.matplotlib.patches.FancyArrowPatch
INK=P.INK; SUB=P.MUTED
# Vertical layout is specified in INCHES and converted to figure fractions, so PANEL_H (the drawn
# size of each structure) is held fixed while the dead space between the rows is tuned. Panels are
# height-limited by imshow's equal aspect, so keeping PANEL_H constant keeps the printed structures
# the same size at a given {width=} in the docx.
FW      = 11.5   # figure width -- unchanged
PANEL_H = 2.62   # drawn height of one panel box (matches the original 0.285 x 9.2in layout)
TOP_M   = 0.35   # margin above the top panel row
LBL_H   = 0.85   # space under a panel row for three label lines
ROW_GAP = 0.35   # white gap between the top row's labels and the bottom row's panels
BOT_M   = 0.20   # margin below the bottom row's labels
FH      = TOP_M + PANEL_H + LBL_H + ROW_GAP + PANEL_H + LBL_H + BOT_M
fig=plt.figure(figsize=(FW, FH), facecolor="white")
L,Rm,gap = 0.02,0.02,0.018
pw=(1-L-Rm-2*gap)/3
xpos=[L+i*(pw+gap) for i in range(3)]
ph=PANEL_H/FH
top_pb=(BOT_M+LBL_H+PANEL_H+ROW_GAP+LBL_H)/FH; bot_pb=(BOT_M+LBL_H)/FH
letters=["A","B","C","D","E"]
place=[(0,top_pb,False),(1,top_pb,False),(2,top_pb,False),(1,bot_pb,False),(2,bot_pb,False)]

def labels(xc, pb, above, cond, mac, ke):
    y = pb+ph if above else pb
    s = 1 if above else -1
    med, lo, hi = ke
    ke_str = f"K–E  {med:.2f} [{lo:.2f}–{hi:.2f}] Å"
    order = [
        (cond, 15, INK, "bold"),
        (f"global MAC {mac:.2f}", 12, maccol(mac), "bold"),
        (ke_str, 10, kecol(med), "bold"),
    ]
    offs = [0.50/FH, 0.20/FH, 0.05/FH] if above else [0.28/FH, 0.48/FH, 0.68/FH]
    for (txt,fs,col,fw),d in zip(order,offs):
        fig.text(xc, y+s*d, txt, color=col, fontsize=fs, fontweight=fw, ha="center")

def harrow(xg, y, label):
    fig.text(xg, y+0.2576/FH, label, color=INK, fontsize=8.5, ha="center", va="bottom", linespacing=1.3)
    fig.add_artist(FA((xg-0.012, y), (xg+0.012, y), transform=fig.transFigure,
                      arrowstyle='-|>', mutation_scale=20, color=INK, lw=2.2))

for i,(im,(fn,cond,mac,trans)) in enumerate(zip(imgs,panels)):
    col,pb,above=place[i]; x=xpos[col]
    ax=fig.add_axes([x, pb, pw, ph]); ax.set_axis_off(); ax.imshow(im)
    labels(x+pw/2, pb, above, cond, mac, ke_data[cond])
    fig.text(x+0.004, pb+ph-0.004, letters[i], color=INK, fontsize=15, fontweight="bold",
             ha="left", va="top")

top_mid=top_pb+ph/2; bot_mid=bot_pb+ph/2
harrow((xpos[0]+pw+xpos[1])/2, top_mid, "+ CCNB1")
harrow((xpos[1]+pw+xpos[2])/2, top_mid, "+ ATP")
harrow((xpos[0]+pw+xpos[1])/2, bot_mid, "CAK · Wee1/Myt1\n+ pT14/pY15/pT161")
harrow((xpos[1]+pw+xpos[2])/2, bot_mid, "CDC25\n− pT14/pY15")

fig.savefig(f"{OUTDIR}/Figure3_CDK1_states.png", dpi=200, facecolor="white")
fig.savefig(f"{OUTDIR}/Figure3_CDK1_states.pdf", facecolor="white")
print("wrote Figure3_CDK1_states.{png,pdf} (white)  crop y",y0,y1,"x",x0,x1)
