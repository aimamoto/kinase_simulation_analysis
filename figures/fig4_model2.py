#!/usr/bin/env python
"""Figure 4 — Model 2 CSK-SRC handshake. Valid (thermodynamically favored, biological) conditions only.
A: reciprocal bulk MAC (CSK stiffens on pY419; SRC bulk-flat). B: pY419 flips SRC state, not stiffness.
C: CSK partial activation. D: single-structure state handshake — two heatmaps (unprimed vs primed).
Data: 260718_rerun_260606_v7_csk-src Module 2 outputs."""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)

B = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
def short(x): return x.replace("\n", "/").replace("wtcat-", "").split("[")[0].strip().rstrip("/")

csk = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv")
src = pd.read_csv(f"{B}/plots_and_stats_SRC_GMM/Phase7_Complete_Structural_Metadata.csv")
csk["cond"] = csk["Condition_reviewed"].map(short); src["cond"] = src["Condition_reviewed"].map(short)
macC = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase5_Global_Network_Density.csv"); macC["c"]=macC["Condition"].map(short)
macS = pd.read_csv(f"{B}/plots_and_stats_SRC_GMM/Phase5_Global_Network_Density.csv"); macS["c"]=macS["Condition"].map(short)
riS = pd.read_csv(f"{B}/plots_and_stats_SRC_GMM/Phase6b_State_Intrinsic_Rigidity_MAC.csv").set_index("Macro_State")["Intrinsic_MAC"]
riC = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase6b_State_Intrinsic_Rigidity_MAC.csv").set_index("Macro_State")["Intrinsic_MAC"]
def macv(df, key):
    r = df[df["c"] == key]; return float(r["Global_Coupling_Score"].iloc[0]) if len(r) else np.nan
allmac = pd.concat([riC, riS]); norm = plt.Normalize(allmac.min(), allmac.max())
def rc(m): return P.RIGIDITY(norm(m))

# VALID biological conditions (exclude apo/src-holo overwrite artifact + src-py159-apo non-biological)
CONDS = [("csk-apo/src-apo","apo/apo"), ("csk-holo/src-apo","CSK-ATP"),
         ("csk-holo/src-holo","both-ATP"), ("csk-holo/src-py159-holo","primed\n(pY419)")]
Acsk = [macv(macC, k) for k,_ in CONDS]; Asrc = [macv(macS, k) for k,_ in CONDS]

fig = plt.figure(figsize=(7.4, 6.6))
gs = fig.add_gridspec(2, 2, wspace=0.32, hspace=0.52, left=0.09, right=0.975, top=0.90, bottom=0.10)
axA, axB, axC = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1]), fig.add_subplot(gs[1,0])
gsD = gs[1,1].subgridspec(2, 2, height_ratios=[1, 0.055], hspace=0.85, wspace=0.12)
axD1, axD2 = fig.add_subplot(gsD[0,0]), fig.add_subplot(gsD[0,1])
axDCB = fig.add_subplot(gsD[1, :])

# ---- A: reciprocal bulk MAC ----
x = np.arange(len(CONDS)); w = 0.38
axA.bar(x-w/2, Acsk, w, color=P.BLUE,   edgecolor=P.INK, lw=0.6, label="CSK (receiver)")
axA.bar(x+w/2, Asrc, w, color=P.PURPLE, edgecolor=P.INK, lw=0.6, label="SRC (loader)")
axA.set_xticks(x); axA.set_xticklabels([l for _,l in CONDS], fontsize=P.TS["small"])
axA.set_ylabel("Global network rigidity (MAC)", fontsize=P.TS["small"]); axA.set_ylim(0,0.30)
axA.legend(frameon=False, fontsize=P.TS["small"], loc="upper left")
axA.annotate("CSK stiffens\n***", xy=(3-w/2, Acsk[3]), xytext=(2.05,0.265), fontsize=P.TS["small"],
             color=P.BLUE_D, ha="center", arrowprops=dict(arrowstyle="->", color=P.BLUE_D, lw=1.1))
axA.plot([-0.4, 3.4], [np.mean(Asrc)]*2, ls=":", color=P.PURPLE, lw=0.9, alpha=0.7)
axA.text(3.4, np.mean(Asrc)+0.006, "SRC bulk-flat ~0.12", ha="right", va="bottom",
         fontsize=P.TS["small"], color=P.PURPLE)
axA.text(-0.34,1.045,"A",transform=axA.transAxes,fontsize=P.TS["panel"],fontweight="bold")
axA.text(0.5,1.045,"Bulk network rigidity of both partners",transform=axA.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",ha="center")

# ---- B: SRC state redistribution ----
Bconds=[("csk-holo/src-holo","both-ATP"),("csk-holo/src-py159-holo","+pY419")]
src_states=["State 1","State 2","State 3","State 5","State 6","State 8"]
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
axB.set_ylabel("SRC macro-state (%)",fontsize=P.TS["small"]); axB.set_ylim(0,100)
axB.text(-0.30,1.045,"B",transform=axB.transAxes,fontsize=P.TS["panel"],fontweight="bold")
axB.text(0.5,1.045,"pY419 flips SRC state, not stiffness",transform=axB.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",ha="center")

# ---- C: CSK partial activation ----
def cat(cond):
    s=csk[csk["cond"]==cond]
    return [ (s["C_Helix"].astype(str).str.contains("In")).mean()*100,
             (s["Spatial"].astype(str).str.contains("in")).mean()*100,
             (s["Dihedral"].astype(str).str.contains("BLAminus")).mean()*100 ]
Cconds=[("csk-apo/src-apo","apo/apo"),("csk-holo/src-apo","CSK-ATP"),
        ("csk-holo/src-holo","both-ATP"),("csk-holo/src-py159-holo","primed")]
feats=["αC-In","DFG-in","BLAminus"]; xc=np.arange(len(feats)); wc=0.2
ramp=[P.GRID,"#AED4F0",P.BLUE,P.BLUE_D]
for i,(k,lab) in enumerate(Cconds):
    axC.bar(xc+(i-1.5)*wc, cat(k), wc, color=ramp[i], edgecolor=P.INK, lw=0.5, label=lab)
axC.set_xticks(xc); axC.set_xticklabels(feats,fontsize=P.TS["small"]); axC.set_ylim(0,100)
axC.set_ylabel("CSK models (%)",fontsize=P.TS["small"])
axC.legend(frameon=False,fontsize=P.TS["small"],loc="upper right",ncol=1,handlelength=1)
axC.text(-0.30,1.045,"C",transform=axC.transAxes,fontsize=P.TS["panel"],fontweight="bold")
axC.text(0.5,1.045,"CSK partial activation",transform=axC.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",ha="center")

# ---- D: two heatmaps (unprimed vs primed), shared SRC x CSK grid ----
SRC_ROWS=["State 8","State 6","State 5"]; CSK_COLS=["State 5","State 6","State 7","State 8","State 9"]
def contin(cond):
    m=csk[csk["cond"]==cond][["Simulation_ID","Macro_State"]].rename(columns={"Macro_State":"CSK"}).merge(
      src[src["cond"]==cond][["Simulation_ID","Macro_State"]].rename(columns={"Macro_State":"SRC"}),on="Simulation_ID")
    return pd.crosstab(m["SRC"],m["CSK"]).reindex(index=SRC_ROWS,columns=CSK_COLS).fillna(0)
Dun, Dpr = contin("csk-holo/src-holo"), contin("csk-holo/src-py159-holo")
vmax=max(Dun.values.max(), Dpr.values.max())
im=None
for ax, M, ttl, stat in [(axD1,Dun,"unprimed","V = 0.30, n.s."),(axD2,Dpr,"primed","V = 0.53, p = 0.01")]:
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
axD1.text(-0.55,1.42,"D",transform=axD1.transAxes,fontsize=P.TS["panel"],fontweight="bold")
axD1.text(1.12,1.42,"Single-structure state handshake",transform=axD1.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",ha="center")

for ext in ("png","pdf"):
    fig.savefig(f"Figure4_CSK_SRC_handshake.{ext}",dpi=200 if ext=="png" else None,bbox_inches="tight")
print("saved Figure4_CSK_SRC_handshake.{png,pdf}")
print("A CSK:",[round(v,3) for v in Acsk],"| A SRC:",[round(v,3) for v in Asrc])
print("D unprimed:\n",Dun.astype(int)); print("D primed:\n",Dpr.astype(int))
