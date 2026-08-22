#!/usr/bin/env python
"""Figure S3 — categorical activation status of both partners, by condition and by metastable state.
A: SRC by condition (flat). B: CSK by condition (climbs). C: SRC by state (flat — the pY419
redistribution of Fig. 4B moves SRC between categorically identical states). D: CSK by state
(all-or-none — the state partition behind the ensemble values; cf. Supplementary Table 4).
Conditions: the four biological ones in Fig. 4A order, then the three excluded designs (hatched).
DFG-in is omitted: 100% in every condition and every state (stated in the legend).
Data: 260718_rerun_260606_v7_csk-src Module 2 outputs (AlloQuant v7r3)."""
import sys, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)

B = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
def short(x): return x.replace("\n", "/").replace("wtcat-", "").split("[")[0].strip().rstrip("/")

csk = pd.read_csv(f"{B}/plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv")
src = pd.read_csv(f"{B}/plots_and_stats_SRC_GMM/Phase7_Complete_Structural_Metadata.csv")
csk["cond"] = csk["Condition_reviewed"].map(short); src["cond"] = src["Condition_reviewed"].map(short)

# Fig. 4A biological order first, then the excluded designs (Suppl. Table 2 "CSK/SRC" labels)
BIO = [("csk-apo/src-apo", "apo/apo"), ("csk-holo/src-apo", "CSK-ATP"),
       ("csk-holo/src-holo", "both-ATP"), ("csk-holo/src-py159-holo", "primed")]
EXC = [("csk-apo/src-py159-apo", "apo/pY"), ("csk-apo/src-py159-holo", "apo/pY+ATP"),
       ("csk-holo/src-py159-apo", "ATP/pY")]
CONDS = BIO + EXC
GAP = 0.85                                   # visual break between the two blocks
XC = np.array([0, 1, 2, 3] + [4 + GAP, 5 + GAP, 6 + GAP], dtype=float)

FEATS = [("αC-In", "C_Helix", "In"), ("Active BLAminus", "Dihedral", "BLAminus")]
COL = [P.FEATURE["αC-In"], P.FEATURE["BLAminus"]]
W = 0.36

def pct(sub, col, token):
    return (sub[col].astype(str).str.contains(token)).mean() * 100 if len(sub) else np.nan

def by_cond(df, key, col, token):
    return pct(df[df["cond"] == key], col, token)

def states_of(df):                            # numeric order, matching Supplementary Table 4
    return sorted(df["Macro_State"].unique(), key=lambda s: int(str(s).split()[-1]))

fig = plt.figure(figsize=(7.5, 6.3))
gs = fig.add_gridspec(2, 2, wspace=0.22, hspace=0.62, left=0.085, right=0.985, top=0.858, bottom=0.100)
axes = [[fig.add_subplot(gs[r, c]) for c in range(2)] for r in range(2)]

def style(ax, letter, title, ylab, xs, labels, ns, rot=0, nrow=-0.14, ttop=1.06):
    ax.set_ylim(0, 108); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel(ylab, fontsize=P.TS["small"], labelpad=2)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, fontsize=7.0, rotation=rot,
                       ha="right" if rot else "center", rotation_mode="anchor" if rot else None)
    ax.set_xlim(xs[0] - 0.62, xs[-1] + 0.62)
    for x, n in zip(xs, ns):
        ax.text(x, nrow, f"{n}", ha="center", va="top", fontsize=6.5, color=P.MUTED,
                transform=ax.get_xaxis_transform(), clip_on=False)
    ax.text(-0.145, ttop, letter, transform=ax.transAxes, fontsize=P.TS["panel"], fontweight="bold")
    ax.text(0.5, ttop, title, transform=ax.transAxes, fontsize=9.0,
            fontweight="bold", ha="center", va="bottom")

def zero_tag(ax, x, v):
    if v == 0:
        ax.text(x, 1.2, "0", ha="center", va="bottom", fontsize=6.8, color=P.MUTED)

# ---- row 1: by condition ----
for ci, (df, who) in enumerate([(src, "SRC"), (csk, "CSK")]):
    ax = axes[0][ci]
    ns = [len(df[df["cond"] == k]) for k, _ in CONDS]
    for fi, (fl, col, tok) in enumerate(FEATS):
        vals = [by_cond(df, k, col, tok) for k, _ in CONDS]
        for j, (x, v) in enumerate(zip(XC, vals)):
            ax.bar(x + (fi - 0.5) * W, v, W, color=COL[fi], edgecolor=P.INK, lw=0.5,
                   hatch="//" if j >= len(BIO) else None,
                   label=fl if (j == 0 and fi in (0, 1)) else None)
            zero_tag(ax, x + (fi - 0.5) * W, v)
    sep = (XC[len(BIO) - 1] + XC[len(BIO)]) / 2
    ax.axvline(sep, color=P.GRID, lw=0.9, ls=(0, (3, 3)), zorder=1)
    ax.text((XC[0] + XC[len(BIO) - 1]) / 2, 1.015, "biological", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=6.8, color=P.MUTED, style="italic")
    ax.text((XC[len(BIO)] + XC[-1]) / 2, 1.015, "excluded\n(non-biological / artifact)",
            transform=ax.get_xaxis_transform(), ha="center", va="bottom", linespacing=1.25,
            fontsize=6.8, color=P.MUTED, style="italic")
    style(ax, "AB"[ci], f"{who} {'loader' if who=='SRC' else 'receiver'} — by condition",
          f"{who} models (%)", XC, [l for _, l in CONDS], ns, rot=32, nrow=-0.30, ttop=1.15)

# ---- row 2: by metastable state ----
for ci, (df, who) in enumerate([(src, "SRC"), (csk, "CSK")]):
    ax = axes[1][ci]
    sts = states_of(df); xs = np.arange(len(sts), dtype=float)
    ns = [int((df["Macro_State"] == s).sum()) for s in sts]
    for fi, (fl, col, tok) in enumerate(FEATS):
        vals = [pct(df[df["Macro_State"] == s], col, tok) for s in sts]
        for x, v in zip(xs, vals):
            ax.bar(x + (fi - 0.5) * W, v, W, color=COL[fi], edgecolor=P.INK, lw=0.5)
            zero_tag(ax, x + (fi - 0.5) * W, v)
    style(ax, "CD"[ci], f"{who} {'loader' if who=='SRC' else 'receiver'} — by metastable state",
          f"{who} models (%)", xs, [str(s).replace("State ", "S") for s in sts], ns)

h = [plt.Rectangle((0, 0), 1, 1, fc=COL[i], ec=P.INK, lw=0.5) for i in range(2)]
h.append(plt.Rectangle((0, 0), 1, 1, fc=P.WHITE, ec=P.INK, lw=0.5, hatch="//"))
fig.legend(h, [FEATS[0][0], FEATS[1][0], "excluded design"], frameon=False,
           fontsize=P.TS["small"], ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.0))
fig.text(0.5, 0.012, "grey numerals beneath each group give n models", ha="center",
         fontsize=6.8, color=P.MUTED)

out = os.path.join(os.environ.get("ALLOQUANT_FIGDIR", "."), "figS3_activation")
fig.savefig(out + ".png", dpi=300, facecolor="white")
fig.savefig(out + ".pdf", facecolor="white")
print("wrote", out + ".{png,pdf}")
