#!/usr/bin/env python
"""Model 2 paired coupling figure: regulatory-spine conformational mirroring across the CSK-SRC interface.
Data: 260718_rerun_260606_v7_csk-src Phase7 tables (per-model, ipTM/PAE populated).
Panel A = within-condition partial-rho spectrum (specificity). Panel B = condition-centered exemplar scatters.
"""
import sys, numpy as np, pandas as pd
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)

BASE=os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
c=pd.read_csv(f"{BASE}/plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv")
s=pd.read_csv(f"{BASE}/plots_and_stats_SRC_GMM/Phase7_Complete_Structural_Metadata.csv")

# element -> (display label, module)  module: 'spine' regulatory-spine cluster, 'cat' catalytic/active-site
ELEMENTS=[
 ("K105_E121_Dist","αC-loop (K105–E121)","spine"),
 ("V104_RS2_Dist","R-spine (V104–RS2)","spine"),
 ("Shell_M118_M120_Dist","Reg. shell (M118–M120)","spine"),
 ("Spine_Bridge_Dist","Spine bridge","spine"),
 ("D220_HRD_Dist","HRD backbone (D220)","spine"),
 ("K105_N99_Dist","K105–N99","spine"),
 ("Cleft_Gape_Dist","Substrate cleft","cat"),
 ("I150_HRD_Dist","αE floor (I150–HRD)","cat"),
 ("SB_Dist","Reg. salt bridge","cat"),
 ("Y156_N99_Dist","Deep scaffold (Y156)","cat"),
 ("HRD_ATP_Dist","HRD–ATP","cat"),
 ("aCb4_aE_Dist","αC-β4/αE","cat"),
]
cols=[e[0] for e in ELEMENTS]
cc=c[["Simulation_ID","Condition_reviewed","ipTM"]+cols].rename(columns={x:"C_"+x for x in cols}).rename(columns={"Condition_reviewed":"cond"})
ss=s[["Simulation_ID"]+cols].rename(columns={x:"S_"+x for x in cols})
p=cc.merge(ss,on="Simulation_ID")
conds=[g for _,g in p.groupby("cond") if len(g)>=25]

def sp(a,b):
    d=pd.concat([a,b],axis=1).dropna()
    return (stats.spearmanr(d.iloc[:,0],d.iloc[:,1])[0],len(d)) if len(d)>6 else (np.nan,0)
def metaz(vals):
    zs=[];ws=[]
    for r,n in vals:
        if np.isnan(r) or n<8: continue
        zs.append(np.arctanh(np.clip(r,-.999,.999)));ws.append(n-3)
    zs=np.array(zs);ws=np.array(ws);zb=(zs*ws).sum()/ws.sum();se=1/np.sqrt(ws.sum())
    return np.tanh(zb),2*(1-stats.norm.cdf(abs(zb/se)))
def part(g,x,y,z):
    rxy,n=sp(g[x],g[y]);rxz,_=sp(g[x],g[z]);ryz,_=sp(g[y],g[z])
    if any(np.isnan([rxy,rxz,ryz])): return (np.nan,0)
    den=np.sqrt((1-rxz**2)*(1-ryz**2)); return ((rxy-rxz*ryz)/den,n) if den>1e-9 else (np.nan,0)

res=[]
for col,lab,mod in ELEMENTS:
    rho,pv=metaz([part(g,"S_"+col,"C_"+col,"ipTM") for g in conds])
    res.append(dict(col=col,lab=lab,mod=mod,rho=rho,p=pv))
res=pd.DataFrame(res).sort_values("rho",ascending=True).reset_index(drop=True)

# ---- figure ----
fig=plt.figure(figsize=(7.2,3.9))
gs=fig.add_gridspec(2,2,width_ratios=[1.35,1.0],height_ratios=[1,1],wspace=0.42,hspace=0.55,left=0.30,right=0.965,top=0.90,bottom=0.135)
axA=fig.add_subplot(gs[:,0]); axB1=fig.add_subplot(gs[0,1]); axB2=fig.add_subplot(gs[1,1])

# Panel A: coupling spectrum
ypos=np.arange(len(res))
for i,r in res.iterrows():
    col=P.BLUE_D if r["mod"]=="spine" else P.GRID
    axA.barh(i,r["rho"],color=col,edgecolor=P.INK,linewidth=0.6,height=0.72,zorder=3)
    star="***" if r["p"]<1e-3 else ("**" if r["p"]<1e-2 else ("*" if r["p"]<0.05 else "n.s."))
    axA.text(r["rho"]+0.006,i,star,va="center",ha="left",fontsize=P.TS["small"],color=P.INK)
axA.set_yticks(ypos); axA.set_yticklabels(res["lab"],fontsize=P.TS["tick"])
for tick,mod in zip(axA.get_yticklabels(),res["mod"]):
    tick.set_color(P.INK if mod=="spine" else P.MUTED)
    if mod=="spine": tick.set_fontweight("bold")
axA.axvline(0,color=P.INK,linewidth=0.8)
axA.set_xlim(-0.03,0.34); axA.set_xlabel("Cross-interface coupling  (partial Spearman ρ)",fontsize=P.TS["small"])
axA.text(-0.34,1.045,"A",transform=axA.transAxes,fontsize=P.TS["panel"],fontweight="bold",va="bottom",ha="left")
axA.text(0.5,1.045,"Cross-interface coupling spectrum",transform=axA.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",va="bottom",ha="center")
axA.grid(axis="y",visible=False)
# legend chips
axA.scatter([],[],marker="s",s=42,color=P.BLUE_D,label="Regulatory-spine module")
axA.scatter([],[],marker="s",s=42,color=P.GRID,label="Catalytic / active-site")
axA.legend(loc="lower right",frameon=False,fontsize=P.TS["small"],handletextpad=0.4)

# Panel B: condition-centered exemplar scatters (two strongest)
def centered(col):
    x=p["S_"+col]-p.groupby("cond")["S_"+col].transform("median")
    y=p["C_"+col]-p.groupby("cond")["C_"+col].transform("median")
    return x,y
def draw(ax,col,lab):
    x,y=centered(col); d=pd.concat([x,y],axis=1).dropna(); x,y=d.iloc[:,0],d.iloc[:,1]
    rho,pv=metaz([sp(g["S_"+col],g["C_"+col]) for g in conds])
    ax.axhline(0,color=P.GRID,lw=0.6,zorder=1); ax.axvline(0,color=P.GRID,lw=0.6,zorder=1)
    ax.scatter(x,y,s=8,color=P.BLUE,alpha=0.5,edgecolor="none",zorder=2)
    b,a=np.polyfit(x,y,1); xs=np.array([x.min(),x.max()])
    ax.plot(xs,a+b*xs,color=P.RED_D,lw=1.6,zorder=3)
    ax.text(0.04,0.92,lab,transform=ax.transAxes,fontsize=P.TS["small"],fontweight="bold",va="top")
    pstr = "p<10⁻¹⁵" if pv<1e-15 else f"p={pv:.0e}"
    ax.text(0.04,0.80,f"ρ={rho:.2f}, {pstr}",transform=ax.transAxes,fontsize=P.TS["small"],va="top",color=P.RED_D)
    ax.tick_params(labelsize=P.TS["small"])
draw(axB1,"K105_E121_Dist","K105–E121 αC-loop network")
axB1.text(-0.28,1.045,"B",transform=axB1.transAxes,fontsize=P.TS["panel"],fontweight="bold",va="bottom",ha="left")
axB1.text(0.5,1.045,"Within-condition co-variation",transform=axB1.transAxes,fontsize=P.TS["axlabel"],fontweight="bold",va="bottom",ha="center")
draw(axB2,"V104_RS2_Dist","R-spine")
axB2.set_xlabel("SRC element  (Å, condition-centered)",fontsize=P.TS["small"])
for ax in (axB1,axB2): ax.set_ylabel("CSK  (Å, cent.)",fontsize=P.TS["small"])

for ext in ("png","pdf"):
    fig.savefig(f"Figure5_CSK_SRC_coupling.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
print("saved Figure5_CSK_SRC_coupling.{png,pdf}")
print(res[["lab","rho","p"]].to_string(index=False))
