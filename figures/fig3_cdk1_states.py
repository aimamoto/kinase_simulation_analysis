#!/usr/bin/env python3
"""Fig 3 — CDK1 activation cycle as an ordered structural series (WHITE bg).
Five ensemble panels (10-model Ca ropes + annotated representative) in biological
order, laid out 3-over-2, with transition arrows and global-MAC labels."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
import numpy as np, sys
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P

D = os.environ.get("ALLOQUANT_STRUCT", "struct/")  # ChimeraX-rendered structural panels
panels = [
    ("fig3ens_apo_w.png",     "apo",           "State 5",               0.381, None),
    ("fig3ens_ccnb1_w.png",   "+ CCNB1",       "State 8",               0.219, "+ CCNB1"),
    ("fig3ens_atp_w.png",     "+ CCNB1 + ATP", "State 9",               0.125, "+ ATP"),
    ("fig3ens_triplep_w.png", "triple-P",      "States 1/2 · inhibited", 0.165, "+ pT14/pY15/pT161"),
    ("fig3ens_pt161_w.png",   "pT161",         "States 1/2 · active",    0.174, "CDC25 − pT14/pY15"),
]
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

FA = plt.matplotlib.patches.FancyArrowPatch
INK=P.INK; SUB=P.MUTED
fig=plt.figure(figsize=(11.5, 9.2), facecolor="white")
L,Rm,gap = 0.02,0.02,0.018
pw=(1-L-Rm-2*gap)/3
xpos=[L+i*(pw+gap) for i in range(3)]
ph=0.285
top_pb=0.66; bot_pb=0.155
letters=["A","B","C","D","E"]
place=[(0,top_pb,False),(1,top_pb,False),(2,top_pb,False),(1,bot_pb,False),(2,bot_pb,False)]

def labels(xc, pb, above, cond, state, mac):
    y = pb+ph if above else pb
    s = 1 if above else -1
    order = [(cond,15,INK,"bold"),(state,10.5,SUB,"normal"),(f"MAC {mac:.2f}",12,maccol(mac),"bold")]
    offs = [0.068,0.041,0.014] if above else [0.030,0.057,0.086]
    for (txt,fs,col,fw),d in zip(order,offs):
        fig.text(xc, y+s*d, txt, color=col, fontsize=fs, fontweight=fw, ha="center")

def harrow(xg, y, label):
    fig.text(xg, y+0.028, label, color=INK, fontsize=8.5, ha="center", va="bottom", linespacing=1.3)
    fig.add_artist(FA((xg-0.012, y), (xg+0.012, y), transform=fig.transFigure,
                      arrowstyle='-|>', mutation_scale=20, color=INK, lw=2.2))

for i,(im,(fn,cond,state,mac,trans)) in enumerate(zip(imgs,panels)):
    col,pb,above=place[i]; x=xpos[col]
    ax=fig.add_axes([x, pb, pw, ph]); ax.set_axis_off(); ax.imshow(im)
    labels(x+pw/2, pb, above, cond, state, mac)
    fig.text(x+0.004, pb+ph-0.004, letters[i], color=INK, fontsize=15, fontweight="bold",
             ha="left", va="top")

top_mid=top_pb+ph/2; bot_mid=bot_pb+ph/2
harrow((xpos[0]+pw+xpos[1])/2, top_mid, "+ CCNB1")
harrow((xpos[1]+pw+xpos[2])/2, top_mid, "+ ATP")
harrow((xpos[0]+pw+xpos[1])/2, bot_mid, "+ pT14/pY15\n/pT161")
harrow((xpos[1]+pw+xpos[2])/2, bot_mid, "CDC25\n− pT14/pY15")

fig.savefig("Figure3_CDK1_states.png", dpi=200, facecolor="white")
fig.savefig("Figure3_CDK1_states.pdf", facecolor="white")
print("wrote Figure3_CDK1_states.{png,pdf} (white)  crop y",y0,y1,"x",x0,x1)
