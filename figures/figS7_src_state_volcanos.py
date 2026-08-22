#!/usr/bin/env python3
"""Figure S7 — SRC intra-state pairwise volcanos (State 4 vs 7, 4 vs 8, 7 vs 8).
State 4 = ATP-loaded unphosphorylated; States 7/8 = pY419 substates.
Data: Phase 8 Wilcoxon + BH correction from plots_and_stats_SRC_GMM."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P

try:
    import pandas as pd
except ImportError:
    raise SystemExit("Run under conda env 'main' (pandas required)")

SRCB = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
VOLC = f"{SRCB}/plots_and_stats_SRC_GMM/Phase8_Volcanos/"

PAIRS = [
    ("4", "7", "State 4 vs State 7"),
    ("4", "8", "State 4 vs State 8"),
    ("7", "8", "State 7 vs State 8"),
]

# Human-readable metric labels (concise for the plot; column names in parentheses for dict lookup)
METRIC_LABELS = {
    "Interface_C_Lobe_Donor_Dist": "C-lobe\ninterface",
    "Interface_N_Lobe_Rec_Dist":   "N-lobe face\n(unengaged)†",
    "Cleft_Gape_Dist":             "Cleft gape",
    "D1_Dist":                     "D1 anchor",
    "D2_Dist":                     "D2 anchor",
    "SB_Dist":                     "K–E salt bridge\n(SB_Dist)",
    "aCb4_aE_Dist":                "αC-β4/αE",
    "Spine_Bridge_Dist":           "Spine bridge",
    "V104_RS2_Dist":               "V-RS2",
    "I150_HRD_Dist":               "I-HRD",
    "Shell_M118_M120_Dist":        "Shell\nM118-M120",
    "Y156_N99_Dist":               "Y156-N99",
    "K105_E107_Dist":              "K105–E107\n(K105_E107_Dist)",
    "K105_E121_Dist":              "K105–E121\n(K105_E121_Dist)",
    "K105_N99_Dist":               "K-N99",
    "D220_HRD_Dist":               "D220-HRD",
    "Mg_Hijack_Dist":              "Mg hijack",
    "HRD_ATP_Dist":                "HRD-ATP",
    "DFG_ATP_Dist":                "DFG-ATP",
    "PLoop_ATP_Dist":              "P-loop–ATP",
}

# Per-panel label placement overrides: (panel_idx, metric) → (ha, dx)
# "left"  = label to left of dot  (ha="right", dx negative)
# "right" = label to right of dot (ha="left",  dx positive)
LABEL_OVERRIDE = {
    # Panel A
    (0, "Interface_C_Lobe_Donor_Dist"): ("right", -0.12),  # avoid stacking with K-N99
    (0, "V104_RS2_Dist"):               ("left",  +0.12),   # blue dots → right side
    (0, "aCb4_aE_Dist"):                ("left",  +0.12),   # blue dots → right side
    (0, "K105_N99_Dist"):               ("right", -0.12),   # avoid crowding right cluster
    # Panel B
    (1, "Interface_C_Lobe_Donor_Dist"): ("right", -0.12),  # avoid stacking with N-lobe face
    (1, "aCb4_aE_Dist"):                ("left",  +0.12),   # blue dot → right side
    (1, "K105_N99_Dist"):               ("right", -0.12),   # avoid crowding right cluster
}

LABEL_THR = 1.0   # |d| threshold for labelling

fig, axes = plt.subplots(1, 3, figsize=(13, 4.8), facecolor="white")
fig.subplots_adjust(left=0.06, right=0.98, top=0.87, bottom=0.18, wspace=0.36)

SIG_THR   = -np.log10(0.05)

for pidx, (ax, (s1, s2, title)) in enumerate(zip(axes, PAIRS)):
    fn = VOLC + "Stats_State_%s_vs_State_%s.csv" % (s1, s2)
    df = pd.read_csv(fn)

    x   = df["Signed_EffSize"].values
    y   = -np.log10(df["p.adj"].values)
    sig = df["p.adj.signif"].values != "ns"
    exp = sig & (df["Direction"].values == "Expanded")
    com = sig & (df["Direction"].values == "Compressed")
    ns  = ~sig

    # scatter
    ax.scatter(x[ns],  y[ns],  s=26, color=P.MUTED, alpha=0.35, zorder=2, linewidths=0)
    ax.scatter(x[exp], y[exp], s=40, color=P.RED,    zorder=3, linewidths=0)
    ax.scatter(x[com], y[com], s=40, color=P.BLUE_D, zorder=3, linewidths=0)

    # threshold lines
    ax.axhline(SIG_THR, color=P.MUTED, lw=0.8, ls="--", zorder=1)
    ax.axvline(0,       color=P.INK,   lw=0.6, zorder=1)
    ax.axvline( LABEL_THR, color=P.MUTED, lw=0.7, ls=(0,(3,4)), zorder=1)
    ax.axvline(-LABEL_THR, color=P.MUTED, lw=0.7, ls=(0,(3,4)), zorder=1)

    # label large-effect significant points
    for _, row in df[sig & (df["Signed_EffSize"].abs() >= LABEL_THR)].iterrows():
        metric = row["Metric"]
        lbl = METRIC_LABELS.get(metric, metric.replace("_Dist","").replace("_"," "))
        xi = row["Signed_EffSize"]
        yi = -np.log10(row["p.adj"])
        if (pidx, metric) in LABEL_OVERRIDE:
            ha, dx = LABEL_OVERRIDE[(pidx, metric)]
        elif xi < 0 and xi < -1.8:
            ha, dx = "left", 0.10
        elif xi > 0:
            ha, dx = "left", 0.10
        else:
            ha, dx = "right", -0.10
        ax.annotate(lbl, xy=(xi, yi), xytext=(xi + dx, yi),
                    fontsize=7.5, color=P.INK, ha=ha, va="center",
                    arrowprops=dict(arrowstyle="-", color=P.MUTED, lw=0.5,
                                    shrinkA=3, shrinkB=2))

    ax.set_xlabel("Cohen's $d$  (State B − State A)", fontsize=8.5, color=P.INK)
    if ax is axes[0]:
        ax.set_ylabel("−log₁₀($p$.adj, BH)", fontsize=8.5, color=P.INK)
    ax.set_title(title, fontsize=10, fontweight="bold", color=P.INK, pad=5)
    letter = "ABC"[list(axes).index(ax)]
    ax.text(-0.08, 1.08, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", color=P.INK, va="top", ha="left")
    ax.tick_params(labelsize=8, color=P.MUTED)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(P.MUTED)

# shared legend
from matplotlib.patches import Patch
handles = [
    Patch(color=P.RED,    label="Expanded in State B"),
    Patch(color=P.BLUE_D, label="Compressed in State B"),
    Patch(color=P.MUTED,  alpha=0.5, label="n.s."),
]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.00),
           ncol=3, fontsize=8, frameon=False, handlelength=1.2)

fig.text(0.5, 0.005,
         "† N-lobe face distance spans the unengaged side of the dimer (landmark >8 Å); "
         "shown for completeness.",
         ha="center", fontsize=6.5, color=P.MUTED, style="italic")

BASE = os.path.join(os.environ.get("ALLOQUANT_FIGDIR", "."),
                    "FigS7_SRC_state_volcanos")
for ext, kw in [("png", {"dpi": 200}), ("pdf", {})]:
    out = f"{BASE}.{ext}"
    fig.savefig(out, facecolor="white", bbox_inches="tight", **kw)
    print("Wrote", out)
