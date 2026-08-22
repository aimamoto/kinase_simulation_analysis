#!/usr/bin/env python
"""Supplementary figure — what SRC pY419 does to the CSK–SRC docking interface.

A: engagement decreases on priming (buried surface area, <4 Å contacts) while charged
   contacts increase — the interface loosens sterically and tightens ionically.
B: the change is a peripheral rim register-shift; the basic–acidic core clamp is invariant.
C: buried area is a property of SRC's metastable state, not of the phosphomark — within
   State 4 the two conditions are indistinguishable.
D: pY419 never touches CSK (13.6 Å from the nearest CSK atom), so the effect is allosteric.

Data: analysis_iface_py419.json (see analysis_iface_py419.py) joined to the AlloQuant v7r3
Module 2 state assignments. n = 100 models per condition.
"""
import sys, json, collections
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import palette as P
P.apply_rc(matplotlib)

OUTDIR = os.environ.get("ALLOQUANT_FIGDIR", ".")
# written by addons/iface_py419.py -- run that first
IFACE_JSON = os.environ.get("ALLOQUANT_IFACE_JSON", "analysis_iface_py419.json")
B = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
UP = "a-csk-wtcat-holo_b-src-wtcat-holo"
PR = "a-csk-wtcat-holo_b-src-wtcat-py159-holo"
C_UP, C_PR = P.MUTED, P.RED_D          # unprimed = neutral, primed = tense pole (direct-labelled)

d = json.load(open(IFACE_JSON))
state = {(r["Directory"].split("/")[0], r["Directory"].split("/")[1]): r["Macro_State"]
         for _, r in pd.read_csv(
             f"{B}/plots_and_stats_SRC_GMM/Phase7_Complete_Structural_Metadata.csv").iterrows()}


def col(cond, key):
    return np.array([r[key] for r in d[cond]], dtype=float)


def pstr(p):
    if p < 1e-10:
        return "$p$ < 1×10$^{-10}$"
    e = int(np.floor(np.log10(p)))
    return f"$p$ = {p/10**e:.0f}×10$^{{{e}}}$"


fig = plt.figure(figsize=(7.5, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.34,
                      left=0.095, right=0.975, top=0.925, bottom=0.085)

# ---------------------------------------------------------------- A: engagement
axA = fig.add_subplot(gs[0, 0])
METS = [("bsa", "buried area\n(Å$^2$)", 1.0), ("n_contacts_4", "contacts\n< 4 Å", 1.0),
        ("n_saltbridge", "charged pairs\n< 4 Å", 1.0)]
for i, (k, lab, _) in enumerate(METS):
    a, b = col(UP, k), col(PR, k)
    # z-score against the unprimed mean so three different units share one axis
    mu, sd = a.mean(), a.std(ddof=1)
    for j, (v, c) in enumerate(((a, C_UP), (b, C_PR))):
        z = (v - mu) / sd
        x = i + (j - 0.5) * 0.42
        parts = axA.violinplot([z], positions=[x], widths=0.36, showextrema=False)
        for pc in parts["bodies"]:
            pc.set_facecolor(c); pc.set_alpha(0.55); pc.set_edgecolor(c); pc.set_linewidth(0.8)
        axA.plot([x - 0.13, x + 0.13], [z.mean()] * 2, color=c, lw=1.6, solid_capstyle="butt")
    p = mannwhitneyu(a, b).pvalue
    axA.text(i, 4.2, pstr(p), ha="center", va="bottom", fontsize=P.TS["small"], color=P.INK)
axA.axhline(0, color=P.GRID, lw=0.8, zorder=0)
axA.set_xticks(range(len(METS)))
axA.set_xticklabels([m[1] for m in METS])
axA.set_ylabel("z-score vs unprimed mean")
axA.set_ylim(-4.8, 5.9)
axA.set_title("A  Priming loosens the interface sterically\n     and tightens it ionically",
              loc="left", fontweight="bold")
axA.plot([], [], color=C_UP, lw=6, alpha=0.6, label="unprimed (both-ATP)")
axA.plot([], [], color=C_PR, lw=6, alpha=0.6, label="primed (SRC pY419)")
axA.legend(loc="lower left", frameon=False, handlelength=1.1, borderpad=0.1)

# ---------------------------------------------------------------- B: rim register shift
axB = fig.add_subplot(gs[0, 1])


def freq(cond):
    n = len(d[cond]); f = collections.Counter()
    for r in d[cond]:
        f.update(tuple(p) for p in r["pairs"])
    return {k: v / n for k, v in f.items()}


fu, fp = freq(UP), freq(PR)
NAMES = {  # local numbering -> canonical label, from the primed representative
    (91, 243): "D91–R243", (207, 264): "L207–T264", (197, 254): "F197–Y254",
    (203, 265): "P203–S265", (204, 265): "R204–S265", (201, 265): "P201–S265"}
sel = sorted(NAMES, key=lambda k: fp.get(k, 0) - fu.get(k, 0))
y = np.arange(len(sel))
axB.barh(y - 0.19, [fu.get(k, 0) for k in sel], height=0.36, color=C_UP, alpha=0.75,
         label="unprimed", edgecolor="none")
axB.barh(y + 0.19, [fp.get(k, 0) for k in sel], height=0.36, color=C_PR, alpha=0.75,
         label="primed", edgecolor="none")
axB.set_yticks(y); axB.set_yticklabels([NAMES[k] for k in sel])
axB.set_xlabel("contact frequency across the ensemble")
axB.set_xlim(0, 1.0)
axB.set_ylim(len(sel) - 0.4, -1.35)          # headroom for the legend, no bar overlap
axB.legend(loc="upper right", frameon=False, handlelength=1.1, ncol=2,
           columnspacing=1.0, borderpad=0.0, bbox_to_anchor=(1.0, 1.005))
axB.set_title("B  Only the rim re-registers;\n     the core clamp is invariant",
              loc="left", fontweight="bold")
inv = sum(1 for k in set(fu) | set(fp) if fu.get(k, 0) >= 0.99 and fp.get(k, 0) >= 0.99)
axB.text(0.0, -0.235, f"a further {inv} residue pairs sit at frequency 1.00 in both "
         f"conditions\n(the basic–acidic core clamp)", transform=axB.transAxes,
         ha="left", va="top", fontsize=P.TS["small"], color=P.MUTED, linespacing=1.35)

# ---------------------------------------------------------------- C: state, not phosphomark
axC = fig.add_subplot(gs[1, 0])
rec = [(c, state.get((c, r["seed"])), r["bsa"]) for c, rows in d.items() for r in rows
       if state.get((c, r["seed"]))]
df = pd.DataFrame(rec, columns=["cond", "st", "bsa"])
big = [s for s in sorted(df.st.unique(), key=lambda s: int(s.split()[-1]))
       if (df.st == s).sum() >= 25]
for i, s in enumerate(big):
    v = df.bsa[df.st == s].values
    parts = axC.violinplot([v], positions=[i], widths=0.7, showextrema=False)
    for pc in parts["bodies"]:
        pc.set_facecolor(P.BLUE); pc.set_alpha(0.45); pc.set_edgecolor(P.BLUE_D); pc.set_linewidth(0.8)
    axC.plot([i - 0.22, i + 0.22], [v.mean()] * 2, color=P.BLUE_D, lw=1.7)
axC.set_xticks(range(len(big)))
axC.set_xticklabels([s.replace("State ", "S") for s in big])
axC.set_xlabel("SRC metastable state (all conditions pooled)")
axC.set_ylabel("buried surface area (Å$^2$)")
axC.set_title("C  Buried area tracks SRC's state,\n     not the phosphomark", loc="left",
              fontweight="bold")
s4u = df.bsa[(df.st == "State 4") & (~df.cond.str.contains("py159"))]
s4p = df.bsa[(df.st == "State 4") & (df.cond.str.contains("py159"))]
i4 = big.index("State 4")
axC.annotate(f"within State 4: unprimed {s4u.mean():.0f}\nvs primed {s4p.mean():.0f} Å$^2$, "
             f"n.s. ($p$ = {mannwhitneyu(s4u, s4p).pvalue:.2f})",
             xy=(i4, df.bsa[df.st == "State 4"].min()), xytext=(0.02, 0.045),
             textcoords="axes fraction", fontsize=P.TS["small"], color=P.INK,
             ha="left", va="bottom", linespacing=1.35,
             arrowprops=dict(arrowstyle="-", color=P.MUTED, lw=0.8,
                             connectionstyle="angle,angleA=0,angleB=90,rad=0"))
axC.set_ylim(1030, 1375)

# ---------------------------------------------------------------- D: the mark never touches
axD = fig.add_subplot(gs[1, 1])
res = np.concatenate([col(c, "ptr_min_dist") for c in d if "py159" in c])
pho = np.concatenate([col(c, "phos_min_dist") for c in d if "py159" in c])
mind = np.concatenate([col(c, "min_dist") for c in d])
axD.hist(res, bins=22, color=P.MUTED, alpha=0.8, edgecolor="none")
axD.hist(pho, bins=22, color=P.PURPLE, alpha=0.9, edgecolor="none")
axD.axvline(mind.mean(), color=P.INK, lw=1.2, ls="--")
axD.set_xlabel("distance to the nearest CSK atom (Å)")
axD.set_ylabel("primed models")
axD.set_xlim(0, 24)
axD.set_title("D  pY419 never contacts CSK", loc="left", fontweight="bold")
axD.text(mind.mean() + 0.5, axD.get_ylim()[1] * 0.97,
         f"closest CSK–SRC\ncontact, {mind.mean():.1f} Å", ha="left", va="top",
         fontsize=P.TS["small"], color=P.INK, linespacing=1.35)
axD.text(res.mean() + 0.7, axD.get_ylim()[1] * 0.52,
         f"pY419 residue\n{res.mean():.1f} ± {res.std(ddof=1):.1f} Å", ha="left", va="center",
         fontsize=P.TS["small"], color=P.MUTED, linespacing=1.35)
axD.text(pho.mean() - 0.7, axD.get_ylim()[1] * 0.78,
         f"phosphate group\n{pho.mean():.1f} ± {pho.std(ddof=1):.1f} Å\n"
         f"(closest of {len(pho)}: {pho.min():.1f} Å)", ha="right", va="center",
         fontsize=P.TS["small"], color=P.PURPLE, linespacing=1.35)

for ax in (axA, axB, axC, axD):
    ax.tick_params(length=3, width=0.8)

out = f"{OUTDIR}/figSI_interface"
fig.savefig(out + ".png", bbox_inches="tight")
fig.savefig(out + ".pdf", bbox_inches="tight")
print("wrote", out + ".{png,pdf}")
