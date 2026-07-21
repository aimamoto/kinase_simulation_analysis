#!/usr/bin/env python3
"""Fig S2 (WHITE bg): pT14 phosphate clashes with the ATP β/γ-phosphate.
Ensemble (P-loop trajectories) + tightest-clasher close-up + distance histogram."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
import numpy as np, glob, gemmi, sys
import scipy.stats as st
from PIL import Image
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)
BASE=os.environ.get("ALLOQUANT_CDK1", "CDK1-CCNB1_output")  # AlloQuant CDK1-CCNB1 output dir
S=os.environ.get("ALLOQUANT_STRUCT", "struct/")  # ChimeraX-rendered structural panels

def dists(cond, r14name, r14atoms):
    out=[]
    for f in sorted(glob.glob(f"{BASE}/{cond}/seed-*/model.cif")):
        m=gemmi.read_structure(f)[0]; grp=[]; atpP=[]
        for ch in m:
            for res in ch:
                if res.seqid.num==14 and res.name==r14name:
                    grp=[a.pos for a in res if a.name in r14atoms]
                if res.name=="ATP":
                    atpP+=[a.pos for a in res if a.name in ("PB","PG")]
        if grp and atpP:
            out.append(min(o.dist(q) for o in grp for q in atpP))
    return np.array(out)
d_tp = dists("a-cdk1-pt14-py15-pt161_b-ccnb1-166_1atp","TPO",("O1P","O2P","O3P","OP1","OP2","OP3"))
d_ac = dists("a-cdk1-pt161_b-ccnb1-166_1atp","THR",("OG1",))
pct_tp=(d_tp<4).mean()*100; pct_ac=(d_ac<4).mean()*100
med_tp,med_ac=np.median(d_tp),np.median(d_ac)
q1_tp,q3_tp=np.percentile(d_tp,[25,75]); q1_ac,q3_ac=np.percentile(d_ac,[25,75])
mn_tp,mx_tp=d_tp.min(),d_tp.max()
n_tp,n_ac=len(d_tp),len(d_ac); n_tp_lt4=int((d_tp<4).sum()); n_ac_lt4=int((d_ac<4).sum())
U_mw,p_mw=st.mannwhitneyu(d_tp,d_ac,alternative="two-sided")
OR_f,p_f=st.fisher_exact([[n_tp_lt4,n_tp-n_tp_lt4],[n_ac_lt4,n_ac-n_ac_lt4]])
def _p(p):
    mm,e=f"{p:.1e}".split("e"); e=int(e)
    return f"{mm}×10{str(e).translate(str.maketrans('-0123456789','⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))}"
p_mw_s=_p(p_mw); p_f_s=_p(p_f)
def crop(fn,pad=10,thr=246):
    im=np.array(Image.open(fn).convert("RGB")); mask=(im<thr).any(2); ys,xs=np.where(mask)
    return im[max(0,ys.min()-pad):ys.max()+pad, max(0,xs.min()-pad):xs.max()+pad]

fig=plt.figure(figsize=(7.6,8.6),facecolor="white")
fig.text(0.5,0.975,"A P-loop phosphate clashes with the ATP γ-phosphate",color=P.INK,fontsize=14,fontweight="bold",ha="center")
fig.text(0.5,0.951,"Inhibitory triple-phosphorylated CDK1 — active-site ensemble and pT14→ATP approach trajectory",color=P.MUTED,fontsize=9,ha="center")
# A: ensemble (P-loop trajectories)
axE=fig.add_axes([0.005,0.50,0.60,0.44]); axE.imshow(crop(S+"cage_white_p1.png")); axE.axis("off")
axE.text(0.02,0.98,"A",transform=axE.transAxes,fontsize=P.TS["panel"],fontweight="bold",va="top")
axE.text(0.5,0.0,"11-model ensemble (P-loop trajectories)",transform=axE.transAxes,fontsize=P.TS["small"],ha="center",va="top",color=P.INK)
# B: pT14 trajectory (min->max distance), single fixed ATP; gradient = distance
axC=fig.add_axes([0.605,0.53,0.39,0.41]); axC.imshow(crop(S+"cage_traj_a.png")); axC.axis("off")
axC.text(0.02,0.98,"B",transform=axC.transAxes,fontsize=P.TS["panel"],fontweight="bold",va="top")
axC.text(0.5,-0.01,"pT14 phosphate trajectory vs. fixed ATP (7 models)",transform=axC.transAxes,fontsize=P.TS["small"],ha="center",va="top",color=P.INK)
from matplotlib.colors import Normalize as _Nrm
from matplotlib.cm import ScalarMappable as _SM
caxB=fig.add_axes([0.655,0.485,0.29,0.013]); cbB=plt.colorbar(_SM(_Nrm(2.5,9.5),plt.cm.coolwarm_r),cax=caxB,orientation="horizontal")
cbB.set_label("pT14 → ATP β/γ-P  (Å)",fontsize=P.TS["small"]); cbB.set_ticks([2.5,5,7.5,9.5]); cbB.ax.tick_params(labelsize=P.TS["small"])
# color key
key=[("pT14 phosphate (actor)","magenta"),("pY15 (bystander)","deepskyblue"),("ATP","gold"),
     ("Mg²⁺","#1E8449"),("P-loop trajectories","0.55")]
for i,(lab,c) in enumerate(key):
    yy=0.45-i*0.028
    fig.add_artist(plt.Line2D([0.05],[yy],marker='o',ms=10,color=c,mec=P.INK,mew=0.5,transform=fig.transFigure))
    fig.text(0.075,yy,lab,color=P.INK,fontsize=11,va="center")
# histogram (white)
axi=fig.add_axes([0.60,0.065,0.375,0.25],facecolor="white")
bins=np.linspace(2,12,26)
axi.hist(d_tp,bins=bins,color="magenta",alpha=0.75); axi.hist(d_ac,bins=bins,color="0.6",alpha=0.55)
axi.axvline(4,color=P.INK,ls="--",lw=1); axi.axvline(med_tp,color="magenta",ls=":",lw=1.4); axi.axvline(med_ac,color="0.4",ls=":",lw=1.4)
axi.set_xlabel("min dist to ATP β/γ-phosphate (Å)",color=P.INK,fontsize=9); axi.set_ylabel("models",color=P.INK,fontsize=9)
axi.tick_params(colors=P.INK,labelsize=8)
from matplotlib.patches import Patch; from matplotlib.lines import Line2D
axi.legend(handles=[Patch(fc="magenta",alpha=0.75,label=f"triple-P (pT14): {pct_tp:.0f}% <4Å"),
                    Patch(fc="0.6",alpha=0.55,label=f"active (T14–OH): {pct_ac:.0f}% <4Å"),
                    Line2D([0],[0],color=P.INK,ls="--",lw=1,label="4 Å cutoff"),
                    Line2D([0],[0],color=P.INK,ls=":",lw=1.4,label="median")],
           fontsize=6.6,facecolor="white",edgecolor=P.GRID,labelcolor=P.INK,loc="upper right",handlelength=1.4)
# text block
fig.text(0.05,0.28,
    f"pT14 apposes the ATP β/γ-phosphates in {pct_tp:.0f}% of inhibited\n"
    f"models (<4 Å): median {med_tp:.1f} Å, IQR {q1_tp:.1f}–{q3_tp:.1f} Å,\n"
    f"range {mn_tp:.1f}–{mx_tp:.1f} Å (n={n_tp}). Unphosphorylated T14:\n"
    f"{pct_ac:.0f}% <4 Å, median {med_ac:.1f} Å (n={n_ac}). The distributions\n"
    f"differ significantly (Mann–Whitney p = {p_mw_s}; Fisher\n"
    f"exact on <4 Å, p = {p_f_s}). Mg²⁺ stays coordinated;\n"
    "pY15 points away.",
    color=P.INK,fontsize=10,ha="left",va="top",linespacing=1.5)
fig.savefig("figSI_cage.png",dpi=200,facecolor="white")
fig.savefig("figSI_cage.pdf",facecolor="white")
print(f"wrote figSI_cage white  triple-P {pct_tp:.0f}%<4A n={n_tp}; active {pct_ac:.0f}%<4A n={n_ac}")
