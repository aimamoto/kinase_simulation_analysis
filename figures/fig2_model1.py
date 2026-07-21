"""Figure 2 — Model 1: CCNB1 thaws CDK1, phosphorylation re-tensions it.
Reads deposited v7r2 CDK1 outputs directly (reproducible). LOCKED TEMPLATE for Figs 3-5.
Panel order: A MAC path | B state rigidity | C state flow | D categorical (x follows A)."""
import os, csv, re, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
TS = P.TS
CDK = os.environ.get("ALLOQUANT_CDK1", "CDK1-CCNB1_output")  # AlloQuant CDK1-CCNB1 output dir
GMM = f"{CDK}/plots_and_stats_CDK1_GMM"
P.apply_rc(matplotlib)
def norm(s): return s.replace("\n", " ").strip()

# ---- data ----------------------------------------------------------------
macs = {norm(r["Condition"]): float(r["Global_Coupling_Score"])
        for r in csv.DictReader(open(f"{GMM}/Phase5_Global_Network_Density.csv"))}
smac = {r["Macro_State"]: float(r["Intrinsic_MAC"])
        for r in csv.DictReader(open(f"{GMM}/Phase6b_State_Intrinsic_Rigidity_MAC.csv"))}
states_sorted = sorted(smac, key=lambda s: -smac[s])
comp = collections.defaultdict(collections.Counter)
for r in csv.DictReader(open(f"{GMM}/Phase6_State_Assignments.csv")):
    comp[norm(r["Condition_reviewed"])][r["Macro_State"]] += 1
feat = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
for r in csv.DictReader(open(f"{CDK}/master_kinase_analysis_results_v7r2.csv")):
    if "CDK1" not in r["Type"]: continue
    c = re.match(r"a-([^/]+)", r["Directory"]).group(1)
    feat[c]["CH"][r["C_Helix"]] += 1; feat[c]["SP"][r["Spatial"]] += 1; feat[c]["DI"][r["Dihedral"]] += 1
def pct(c, grp, key):
    d = feat[c][grp]; n = sum(d.values()); return 100*d.get(key,0)/n if n else 0

# biological order, shared by A/C/D; short labels
PATH = [  # (Phase5 key, Phase6 key, master dir key, bar label, flow label)
    ("cdk1 [CDK1 = Apo]",                        "cdk1 [0atp]",                       "cdk1_0atp",                             "apo",            "apo"),
    ("cdk1 ccnb1 [CDK1 = Apo]",                  "cdk1 ccnb1 [0atp]",                 "cdk1_b-ccnb1-166_0atp",                 "·CCNB1",         "·CCNB1"),
    ("cdk1 ccnb1 [CDK1 = Holo]",                 "cdk1 ccnb1 [1atp]",                 "cdk1_b-ccnb1-166_1atp",                 "·CCNB1+ATP",     "·CCNB1+ATP"),
    ("cdk1-pt14-py15-pt161 ccnb1 [CDK1 = Holo]", "cdk1-pt14-py15-pt161 ccnb1 [1atp]", "cdk1-pt14-py15-pt161_b-ccnb1-166_1atp", "triple-P",       "triple-P"),
    ("cdk1-pt161 ccnb1 [CDK1 = Holo]",           "cdk1-pt161 ccnb1 [1atp]",           "cdk1-pt161_b-ccnb1-166_1atp",           "pT161 (active)", "pT161"),
]
CTRL = ("cdk1 [CDK1 = Holo]", "cdk1 [1atp]", "cdk1_1atp", "+ATP (no cyclin)", "+ATP mono")
norm_mac = lambda v: (v-0.10)/(0.40-0.10)

# ---- figure --------------------------------------------------------------
fig = plt.figure(figsize=(7.5, 7.2))
gs = fig.add_gridspec(2, 2, hspace=0.85, wspace=0.30, left=0.13, right=0.985, top=0.9, bottom=0.17)
def tag(ax, s):
    ax.text(-0.20, 1.13, s, transform=ax.transAxes, fontsize=TS["panel"], fontweight="bold", va="top", ha="left", color=P.INK)

# Panel A — MAC along the biological activation path -----------------------
axA = fig.add_subplot(gs[0, 0]); tag(axA, "A")
xs = list(range(len(PATH))); hs = [macs[k[0]] for k in PATH]
axA.bar(xs, hs, width=0.66, color=[P.RIGIDITY(norm_mac(v)) for v in hs], edgecolor=P.INK, linewidth=0.6, zorder=3)
for xi, v in zip(xs, hs):
    axA.text(xi, v+0.008, f"{v:.3f}", ha="center", va="bottom", fontsize=TS["value"]-0.5, color=P.INK, zorder=4)
cx, cv = 6.0, macs[CTRL[0]]
axA.bar([cx], [cv], width=0.66, color=P.FILL_LT, edgecolor=P.MUTED, linewidth=0.9, hatch="////", zorder=3)
axA.text(cx, cv+0.008, f"{cv:.3f}", ha="center", va="bottom", fontsize=TS["value"]-0.5, color=P.MUTED, zorder=4)
# step arrows + p.adj (explicit label placement to avoid overlaps)
STEPS = [(0,1,"****","1.7×10⁻⁶",P.RED_D,"dn",0.6,0.46),(1,2,"****","1.1×10⁻⁷",P.RED_D,"dn",1.55,0.32),
         (2,3,"*","0.013",P.GREEN_D,"up",2.5,0.255),(3,4,"n.s.","0.29",P.MUTED,"flat",3.6,0.235)]
for i,j,sig,pv,col,d,lx,ly in STEPS:
    axA.annotate("", xy=(j-0.24, hs[j]+0.015), xytext=(i+0.24, hs[i]+0.015),
                 arrowprops=dict(arrowstyle="-|>", color=col, lw=1.5,
                                 connectionstyle=f"arc3,rad={-0.3 if d!='flat' else -0.6}"), zorder=5)
    axA.text(lx, ly, f"{sig}\n{pv}", ha="center", va="bottom", fontsize=TS["small"]-0.8,
             color=col, fontweight="bold", zorder=6, linespacing=1.0)
axA.set_xticks(xs+[cx])
_tA = axA.set_xticklabels([k[3] for k in PATH]+[CTRL[3]], rotation=45, ha="right", rotation_mode="anchor")
_tA[-1].set_color(P.MUTED)
axA.set_ylabel("Global network rigidity (MAC)"); axA.set_ylim(0, 0.55); axA.set_xlim(-0.7, 6.7)
axA.set_title("Rigidity along the biological activation path", pad=6)
axA.text(0.75, 0.51, "thaw", color=P.RED_D, fontsize=TS["anno"], ha="center", style="italic")
axA.text(2.9, 0.31, "activation", color=P.GREEN_D, fontsize=TS["anno"]-0.5, ha="center", style="italic")

# Panel B — meta-stable state rigidity -------------------------------------
axB = fig.add_subplot(gs[0, 1]); tag(axB, "B")
xsb = np.arange(len(states_sorted)); valsb = [smac[s] for s in states_sorted]
axB.bar(xsb, valsb, width=0.7, color=[P.RIGIDITY(norm_mac(v)) for v in valsb], edgecolor=P.INK, linewidth=0.6)
axB.set_xticks(xsb); axB.set_xticklabels([s.replace("State ","S") for s in states_sorted])
axB.set_ylabel("Intrinsic state rigidity (MAC)"); axB.set_ylim(0, 0.46)
axB.set_title("Meta-stable state rigidity", pad=6)
notes = {"State 5":"trapped\ninactive","State 9":"fluid\ninterm.","State 2":"tense","State 1":"relaxed"}
for i, s in enumerate(states_sorted):
    if s in notes:
        axB.text(i, valsb[i]+0.009, notes[s], ha="center", va="bottom", fontsize=TS["small"]-1, color=P.MUTED, style="italic", linespacing=1.0)

# Panel C — GMM state flow (biological order, matches A) --------------------
axC = fig.add_subplot(gs[1, 0]); tag(axC, "C")
flow_keys = [k[1] for k in PATH] + [CTRL[1]]; flow_labs = [k[4] for k in PATH] + [CTRL[4]]
y = np.arange(len(flow_keys))[::-1]; left = np.zeros(len(flow_keys))
for st in states_sorted:
    fr = np.array([100*comp[k].get(st,0)/max(sum(comp[k].values()),1) for k in flow_keys])
    axC.barh(y, fr, left=left, height=0.72, color=P.RIGIDITY(norm_mac(smac[st])), edgecolor=P.WHITE, linewidth=0.8)
    for yi, f0, l0 in zip(y, fr, left):
        if f0 >= 24:
            axC.text(l0+f0/2, yi, st.replace("State ","S"), ha="center", va="center", fontsize=TS["small"], color=P.WHITE, fontweight="bold")
    left += fr
axC.set_yticks(y); axC.set_yticklabels(flow_labs); axC.set_xlim(0,100)
axC.set_xlabel("Macro-state population (%)"); axC.set_title("GMM state flow along activation", pad=6)
axC.grid(axis="y", visible=False)

# Panel D — categorical activation signatures (x follows A) -----------------
axD = fig.add_subplot(gs[1, 1]); tag(axD, "D")
w = 0.26; series = [("αC-In","CH","In"),("DFG-in","SP","DFGin"),("BLAminus","DI","BLAminus")]
for jj,(name,grp,key) in enumerate(series):
    axD.bar([x+(jj-1)*w for x in range(len(PATH))], [pct(k[2],grp,key) for k in PATH], width=w,
            color=P.FEATURE[name], edgecolor=P.INK, linewidth=0.4, label=name)
    axD.bar([6+(jj-1)*w], [pct(CTRL[2],grp,key)], width=w, color=P.FEATURE[name], edgecolor=P.MUTED, linewidth=0.4, hatch="////", alpha=0.85)
axD.set_xticks(list(range(len(PATH)))+[6])
_tD = axD.set_xticklabels([k[3] for k in PATH]+[CTRL[3]], rotation=45, ha="right", rotation_mode="anchor")
_tD[-1].set_color(P.MUTED)
axD.set_ylabel("Canonical active feature (%)"); axD.set_ylim(0, 118); axD.set_xlim(-0.7, 6.7)
axD.set_title("Categorical activation signatures", pad=6)
axD.legend(frameon=False, ncol=1, loc="upper right", bbox_to_anchor=(1.0, 1.0),
           handlelength=1.0, labelspacing=0.3, fontsize=TS["legend"]-1)

fig.suptitle("Figure 2  |  Model 1: CCNB1 thaws the trapped CDK1 core; pT161 re-tensions it",
             x=0.02, y=0.985, ha="left", fontsize=12, fontweight="bold", color=P.INK)
out = os.environ.get("ALLOQUANT_FIGDIR", ".")
fig.savefig(f"{out}/Figure2_Model1.pdf"); fig.savefig(f"{out}/Figure2_Model1.png", dpi=200)
print("wrote Figure2_Model1")
