#!/usr/bin/env python
"""Figure 4 — Model 2 CSK-SRC handshake. Valid (thermodynamically favored, biological) conditions only.
A: pooled (condition-level) MAC only. CSK rises on priming; SRC stays bulk-flat, which is real.
The per-state overlay added for register F7 was REMOVED 2026-08-24. MAC is strongly upward-biased at
small n (0.08 at n=100, 0.22 at n=15), so per-state values at n=15-45 and pooled values at n=100 do
not share a linear axis: a shared axis asserts a comparison the statistic does not license, whatever
the legend says. Retaining it made the paper's central figure open on an apparent contradiction.
The composition control (which is what F7 was really after, and what T9/T11/T12 rest on) now lives
in Supplementary Note 10 and Table 7, tested against a size-matched null rather than read off a bar
height: addons/mac_nmatched_control.py. B: pY419 flips SRC state, not stiffness.
C: single-structure state handshake — two heatmaps (unprimed vs primed), spanning the bottom row.
The former panel C (CSK categorical activation by condition) was retired 2026-08-01: it plotted the
ensemble averages the Results now argue understate the mechanism. Superseded by Fig. S3B, with the
state-resolved version in Fig. S3D (figSI_activation.py).
Data: 260718_rerun_260606_v7_csk-src Module 2 outputs (AlloQuant v7r3)."""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)

B = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
OUTDIR = os.environ.get("ALLOQUANT_FIGDIR", ".")  # rendered-figure output dir
def short(x): return x.replace("\n", "/").replace("wtcat-", "").split("[")[0].strip().rstrip("/")

csk = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv")
src = pd.read_csv(f"{B}/plots_and_stats_SRC_GMM/Phase7_Complete_Structural_Metadata.csv")
csk["cond"] = csk["Condition_reviewed"].map(short); src["cond"] = src["Condition_reviewed"].map(short)
macC = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase5_Global_Network_Density.csv"); macC["c"]=macC["Condition"].map(short)
macS = pd.read_csv(f"{B}/plots_and_stats_SRC_GMM/Phase5_Global_Network_Density.csv"); macS["c"]=macS["Condition"].map(short)
riS = pd.read_csv(f"{B}/plots_and_stats_SRC_GMM/Phase6b_State_Intrinsic_Rigidity_MAC.csv").set_index("Macro_State")["Intrinsic_MAC"]
riC = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase6b_State_Intrinsic_Rigidity_MAC.csv").set_index("Macro_State")["Intrinsic_MAC"]
stC = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase5_Global_Network_Density_Stats.csv")
stC["c1"] = stC["Condition_1"].map(short); stC["c2"] = stC["Condition_2"].map(short)
def macv(df, key):
    r = df[df["c"] == key]; return float(r["Global_Coupling_Score"].iloc[0]) if len(r) else np.nan
def padj(df, a, b):
    r = df[((df["c1"] == a) & (df["c2"] == b)) | ((df["c1"] == b) & (df["c2"] == a))]
    return float(r["p.adj"].iloc[0])
def psci(p):  # 3.279e-08 -> "3.3x10^-8" with unicode superscripts
    e = int(np.floor(np.log10(p)))
    sup = str(e).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{p/10**e:.1f}×10{sup}"
allmac = pd.concat([riC, riS]); norm = plt.Normalize(allmac.min(), allmac.max())
def rc(m): return P.RIGIDITY(norm(m))

# VALID biological conditions (exclude apo/src-holo overwrite artifact + src-py159-apo non-biological)
CONDS = [("csk-apo/src-apo","apo/apo"), ("csk-holo/src-apo","CSK-ATP"),
         ("csk-holo/src-holo","both-ATP"), ("csk-holo/src-py159-holo","primed\n(pY419)")]
Acsk = [macv(macC, k) for k,_ in CONDS]; Asrc = [macv(macS, k) for k,_ in CONDS]

fig = plt.figure(figsize=(6.7, 5.9))
gs = fig.add_gridspec(2, 2, width_ratios=[1.42, 1], height_ratios=[1, 1.12],
                      wspace=0.34, hspace=0.60, left=0.100, right=0.975, top=0.90, bottom=0.095)
axA, axB = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])
gsD = gs[1,:].subgridspec(2, 2, height_ratios=[1, 0.055], hspace=0.62, wspace=0.10)
axD1, axD2 = fig.add_subplot(gsD[0,0]), fig.add_subplot(gsD[0,1])
axDCB = fig.add_subplot(gsD[1, :])

# ---- A: pooled (condition-level) MAC ----
# One claim per panel: CSK rises on priming, SRC does not. The per-state overlay is gone (see the
# module docstring); its numbers are Supplementary Table 7, its test Supplementary Note 10.
x = np.arange(len(CONDS)); w = 0.38
axA.bar(x-w/2, Acsk, w, color=P.BLUE,   edgecolor=P.INK, lw=0.6, label="CSK")
axA.bar(x+w/2, Asrc, w, color=P.PURPLE, edgecolor=P.INK, lw=0.6, label="SRC")
axA.set_xticks(x); axA.set_xticklabels([l for _,l in CONDS], fontsize=P.TS["small"])
axA.set_ylabel("Pooled network rigidity (MAC)", fontsize=P.TS["small"]); axA.set_ylim(0, 0.30)
axA.legend(frameon=False, fontsize=P.TS["small"], loc="upper left", ncol=2,
           handletextpad=0.25, borderpad=0.0, columnspacing=0.8,
           bbox_to_anchor=(0.01, 1.00))
# the pY419-isolating contrast keeps its marker; the legend defines the test and the comparator
axA.text(3-w/2, Acsk[3]+0.007, "*", ha="center", va="bottom", fontsize=P.TS["base"], color=P.BLUE_D)
axA.plot([-0.4, 3.4], [np.mean(Asrc)]*2, ls=":", color=P.PURPLE, lw=0.9, alpha=0.7)
axA.text(-0.34,1.045,"A",transform=axA.transAxes,fontsize=P.TS["panel"],fontweight="bold")
# NB: the title must not say "stiffens" -- that is the reading T9/T11/T12 retracted. It describes
# the pooled statistic and nothing about coupling within states.
axA.text(0.5,1.045,"Pooled rigidity rises for CSK only",transform=axA.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",ha="center")

# ---- B: SRC state redistribution ----
Bconds=[("csk-holo/src-holo","both-ATP"),("csk-holo/src-py159-holo","+pY419")]
src_states=["State 1","State 2","State 3","State 4","State 5","State 7","State 8"]
for i,(k,_) in enumerate(Bconds):
    v=src[src["cond"]==k]["Macro_State"].value_counts(normalize=True)*100; bot=0
    for st in src_states:
        h=v.get(st,0)
        if h<=0: continue
        axB.bar(i,h,0.6,bottom=bot,color=rc(riS.get(st,allmac.min())),edgecolor=P.INK,lw=0.4)
        if h>=7: axB.text(i,bot+h/2,st.replace("State ","S"),ha="center",va="center",
                          fontsize=P.TS["small"],color="white" if norm(riS.get(st,0))>0.55 else P.INK)
        bot+=h
axB.set_xticks([0,1]); axB.set_xticklabels([l for _,l in Bconds],fontsize=P.TS["small"])
axB.set_ylabel("SRC metastable state (%)",fontsize=P.TS["small"]); axB.set_ylim(0,100)
axB.text(-0.30,1.045,"B",transform=axB.transAxes,fontsize=P.TS["panel"],fontweight="bold")
axB.text(0.5,1.045,"pY419 flips SRC state, not stiffness",transform=axB.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",ha="center")

# ---- C (was D): two heatmaps (unprimed vs primed), shared SRC x CSK grid ----
SRC_ROWS=["State 4","State 7","State 8"]; CSK_COLS=["State 3","State 5","State 7","State 8","State 9"]
def contin(cond):
    m=csk[csk["cond"]==cond][["Simulation_ID","Macro_State"]].rename(columns={"Macro_State":"CSK"}).merge(
      src[src["cond"]==cond][["Simulation_ID","Macro_State"]].rename(columns={"Macro_State":"SRC"}),on="Simulation_ID")
    return pd.crosstab(m["SRC"],m["CSK"]).reindex(index=SRC_ROWS,columns=CSK_COLS).fillna(0)
Dun, Dpr = contin("csk-holo/src-holo"), contin("csk-holo/src-py159-holo")
vmax=max(Dun.values.max(), Dpr.values.max())
im=None
for ax, M, ttl, stat in [(axD1,Dun,"unprimed","V = 0.26, n.s."),(axD2,Dpr,"primed","V = 0.44, p < 0.001 \u2020")]:
    im=ax.imshow(M.values, cmap=P.RIGIDITY, aspect="auto", vmin=0, vmax=vmax)
    ax.set_xticks(range(len(CSK_COLS))); ax.set_xticklabels([c.replace("State ","C") for c in CSK_COLS],fontsize=P.TS["small"])
    ax.set_yticks(range(len(SRC_ROWS)))
    ax.set_yticklabels([s.replace("State ","S") for s in SRC_ROWS] if ax is axD1 else [], fontsize=P.TS["small"])
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if M.values[i,j]>0:
                ax.text(j,i,int(M.values[i,j]),ha="center",va="center",fontsize=P.TS["small"],
                        color="white" if M.values[i,j]>vmax*0.55 else P.INK)
    ax.set_title(f"{ttl}\n{stat}",fontsize=P.TS["small"])
    ax.set_xlabel("CSK state",fontsize=P.TS["small"])
axD1.set_ylabel("SRC state",fontsize=P.TS["small"])
# shared horizontal colorbar spanning both heatmaps
cb=fig.colorbar(im, cax=axDCB, orientation="horizontal")
cb.set_label("models per cell", fontsize=P.TS["small"])
cb.ax.tick_params(labelsize=P.TS["small"])
cb.outline.set_linewidth(0.5)
axD1.text(-0.27,1.22,"C",transform=axD1.transAxes,fontsize=P.TS["panel"],fontweight="bold")
axD1.text(1.05,1.22,"Single-structure state handshake",transform=axD1.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",ha="center")

for ext in ("png","pdf"):
    fig.savefig(f"{OUTDIR}/Figure4_CSK_SRC_handshake.{ext}",dpi=200 if ext=="png" else None,bbox_inches="tight")
print("saved Figure4_CSK_SRC_handshake.{png,pdf}")
print("A CSK:",[round(v,3) for v in Acsk],"| A SRC:",[round(v,3) for v in Asrc])
# reported for the legend footnote; no longer annotated on the panel itself (register F7)
print("A both-ATP vs primed p.adj:", psci(padj(stC, "csk-holo/src-holo", "csk-holo/src-py159-holo")))
print("D unprimed:\n",Dun.astype(int)); print("D primed:\n",Dpr.astype(int))
