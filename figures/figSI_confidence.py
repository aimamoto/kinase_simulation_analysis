#!/usr/bin/env python
"""Figure S5 — MAC differences are not an AlphaFold3 prediction-confidence artifact.

Visual companion to Supplementary Note 11 / Supplementary Tables 2-3.
A: per-condition MAC, raw -> confidence-partialled (20 conditions, 3 series).
B: the confound is real — AF3 confidence differs strongly between the same conditions.
C: contrast concordance — FDR before vs after partialling, all 57 pairwise contrasts.
D: outcome tally per series.

Source: addons/mac_confidence_control.py outputs (--covar-mode common) in
260718_rerun_260606_v7_csk-src/addons/mac_confidence_out_v7r3/, plus per-model confidence
(CSK-SRC: 260718 rerun master; CDK1: its own Phase7 metadata).
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import palette as P
P.apply_rc(matplotlib)
TS = P.TS

SRCB = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
CDKBASE = os.environ.get("ALLOQUANT_CDK1", "CDK1-CCNB1_output")  # AlloQuant CDK1-CCNB1 output dir
MC = os.environ.get("ALLOQUANT_MAC_CONFIDENCE_OUT",
                    f"{SRCB}/addons/mac_confidence_out_v7r3")
OUT = os.environ.get("ALLOQUANT_FIGDIR", ".")

# ---- labels ---------------------------------------------------------------
def design(cond):
    """'csk-wtcat-holo\\nsrc-wtcat-apo\\n[CSK = Holo]' -> 'csk-holo/src-apo'."""
    return (str(cond).replace("wtcat-", "").replace("\n", "/")
            .split("[")[0].strip().rstrip("/"))

# the four thermodynamically favoured designs carried in main-text Fig. 4
MAIN2 = {"csk-apo/src-apo": "apo / apo", "csk-holo/src-apo": "CSK-ATP",
         "csk-holo/src-holo": "both-ATP", "csk-holo/src-py159-holo": "primed (pY419)"}
# CDK1 conditions, labelled exactly as in main-text Fig. 2A
MAIN1 = {"cdk1 [CDK1 = Apo]": "apo",
         "cdk1 ccnb1 [CDK1 = Apo]": "+CCNB1",
         "cdk1 ccnb1 [CDK1 = Holo]": "+CCNB1+ATP",
         "cdk1-pt14-py15-pt161 ccnb1 [CDK1 = Holo]": "triple-P",
         "cdk1-pt161 ccnb1 [CDK1 = Holo]": "pT161 (active)",
         "cdk1 [CDK1 = Holo]": "+ATP (no cyclin)"}
flat = lambda s: " ".join(str(s).split())

def lab2(cond):
    """CSK/SRC row label; main-text conditions get their Fig. 4 alias."""
    d = design(cond)
    state = str(cond).split("[")[-1].rstrip("]").split("=")[-1].strip()
    return (MAIN2[d], True) if d in MAIN2 else (f"{d} [{state}]", False)

def lab1(cond):
    return (MAIN1[flat(cond)], True) if flat(cond) in MAIN1 else (flat(cond), False)

SERIES = [("CSK", "CSK (receiver)", P.BLUE_D, "CSK_mac_raw_vs_partial.csv",
           "CSK_mac_pairwise_tests.csv", lab2),
          ("SRC", "SRC (loader)", P.PURPLE, "SRC_mac_raw_vs_partial.csv",
           "SRC_mac_pairwise_tests.csv", lab2),
          ("CDK1", "CDK1", P.AMBER, "CDK1_mac_raw_vs_partial.csv",
           "CDK1_mac_pairwise_tests.csv", lab1)]

# ---- per-model confidence -------------------------------------------------
conf_master = (pd.read_csv(f"{SRCB}/master_kinase_analysis_results_v7r3.csv")
               [["Simulation_ID", "ipTM", "PAE_Mean_AB"]]
               .drop_duplicates(subset=["Simulation_ID"]))

def confidence_by_condition(key):
    """(DataFrame[Condition, mean, sd], covariate name, Kruskal-Wallis p) per series."""
    if key == "CDK1":
        df = pd.read_csv(f"{CDKBASE}/plots_and_stats_CDK1_GMM/"
                         "Phase7_Complete_Structural_Metadata.csv")
        col = "pTM"
    else:
        df = pd.read_csv(f"{SRCB}/plots_and_stats_{key}_GMM/Phase7_Complete_Structural_Metadata.csv")
        df = df.drop(columns=["ipTM", "PAE_Mean_AB"], errors="ignore").merge(
            conf_master, on="Simulation_ID", how="left")
        col = "ipTM"
    g = df.groupby("Condition")[col]
    samples = [s.dropna().values for _, s in g if s.notna().sum() > 1]
    kw = stats.kruskal(*samples)[1] if len(samples) > 1 else np.nan
    return g.agg(["mean", "std"]).reset_index(), col, kw

# ---- assemble rows (descending raw MAC within series, as in Supp. Table 2) --
rows, blocks = [], []
for key, series_lab, colr, macf, testf, labfn in SERIES:
    mac = pd.read_csv(f"{MC}/{macf}").sort_values("MAC_raw", ascending=False)
    cbc, ccol, kw = confidence_by_condition(key)
    cbc["_k"] = cbc["Condition"].map(flat)
    cmean = dict(zip(cbc["_k"], cbc["mean"]))
    csd = dict(zip(cbc["_k"], cbc["std"]))
    start = len(rows)
    for r in mac.itertuples(index=False):
        cond = getattr(r, "Condition")
        text, is_main = labfn(cond)
        rows.append(dict(series=key, colour=colr, label=text, main=is_main,
                         raw=r.MAC_raw, part=r.MAC_partial, delta=r.delta, n=r.N,
                         conf=cmean.get(flat(cond), np.nan),
                         conf_sd=csd.get(flat(cond), np.nan)))
    blocks.append(dict(key=key, label=series_lab, colour=colr, ccol=ccol, kw=kw,
                       i0=start, i1=len(rows)))
rows = pd.DataFrame(rows)

# y positions: one blank slot between series blocks
ypos = np.zeros(len(rows))
slot = 0
for b in blocks:
    for i in range(b["i0"], b["i1"]):
        ypos[i] = slot
        slot += 1
    b["yc"] = (ypos[b["i0"]] + ypos[b["i1"] - 1]) / 2
    slot += 1
ytop = slot - 1
ypos = ytop - ypos                      # first row at the top
for b in blocks:
    b["yc"] = ytop - b["yc"]
rows["y"] = ypos

# ---- contrasts ------------------------------------------------------------
tests = []
for key, series_lab, colr, macf, testf, labfn in SERIES:
    t = pd.read_csv(f"{MC}/{testf}")
    t["series"], t["colour"] = key, colr
    tests.append(t)
tests = pd.concat(tests, ignore_index=True)
tests["sig_raw"] = tests.fdr_raw < 0.05
tests["surv"] = (tests.fdr_partial < 0.05) & tests.direction_preserved
tests["lx"] = -np.log10(tests.fdr_raw.clip(lower=1e-18))
tests["ly"] = -np.log10(tests.fdr_partial.clip(lower=1e-18))

# ---- figure ---------------------------------------------------------------
fig = plt.figure(figsize=(7.4, 6.1))
outer = fig.add_gridspec(1, 2, width_ratios=[1.72, 1.10], wspace=0.30,
                         left=0.238, right=0.975, top=0.925, bottom=0.075)
gsL = outer[0].subgridspec(1, 2, width_ratios=[1.40, 0.42], wspace=0.09)
gsR = outer[1].subgridspec(2, 1, height_ratios=[1.16, 1.0], hspace=0.72)
axA = fig.add_subplot(gsL[0, 0])
axB = fig.add_subplot(gsL[0, 1], sharey=axA)
axC = fig.add_subplot(gsR[0, 0])
axD = fig.add_subplot(gsR[1, 0])

def panel(ax, letter, title=None, dx=-0.30):
    ax.text(dx, 1.0, letter, transform=ax.transAxes, fontsize=TS["panel"],
            fontweight="bold", va="bottom", ha="left")
    if title:
        ax.text(0.47, 1.005, title, transform=ax.transAxes, fontsize=TS["small"],
                fontweight="bold", va="bottom", ha="center")

# --- A: raw -> partialled MAC ---
for r in rows.itertuples(index=False):
    axA.plot([r.raw, r.part], [r.y, r.y], color=r.colour, lw=1.5, alpha=0.55,
             solid_capstyle="round", zorder=2)
    axA.plot(r.raw, r.y, "o", ms=5.2, mfc=P.WHITE, mec=r.colour, mew=1.3, zorder=3)
    axA.plot(r.part, r.y, "o", ms=5.2, mfc=r.colour, mec=r.colour, zorder=4)
    if abs(r.delta) >= 0.02:            # only the two confidence-associated ensembles
        axA.annotate(f"Δ {r.delta:+.3f}", (r.part, r.y), textcoords="offset points",
                     xytext=(-7, 0), ha="right", va="center",
                     fontsize=TS["small"] - 0.5, color=P.RED_D, fontweight="bold")
axA.set_yticks(rows["y"])
axA.set_yticklabels(rows["label"], fontsize=TS["small"] - 0.5)
for tick, is_main in zip(axA.get_yticklabels(), rows["main"]):
    tick.set_color(P.INK if is_main else P.MUTED)
    if is_main:
        tick.set_fontweight("bold")
axA.set_ylim(-0.8, ytop + 0.8)
axA.set_xlim(0.085, 0.415)
axA.set_xticks(np.arange(0.10, 0.41, 0.05))
axA.set_xticklabels([f"{v:.2f}" for v in np.arange(0.10, 0.41, 0.05)])
axA.set_xlabel("Global network rigidity (MAC)", fontsize=TS["small"])
axA.grid(axis="y", visible=False)
axA.tick_params(axis="x", labelsize=TS["small"])
for b in blocks:                        # series header in the blank slot above each block
    axA.text(0.092, ypos[b["i0"]] + 0.78, b["label"], va="center", ha="left",
             fontsize=TS["small"], fontweight="bold", color=b["colour"])
axA.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5.2, mfc=P.WHITE,
                           mec=P.INK, mew=1.3, label="raw"),
                    Line2D([], [], marker="o", ls="none", ms=5.2, color=P.INK,
                           label="confidence-partialled")],
           loc="lower right", frameon=False, fontsize=TS["small"] - 0.5,
           handletextpad=0.3, labelspacing=0.25)
panel(axA, "A", "MAC before and after partialling", dx=-0.365)

# --- B: the confound is real ---
for r in rows.itertuples(index=False):
    axB.errorbar(r.conf, r.y, xerr=r.conf_sd, fmt="s", ms=3.4, color=r.colour,
                 ecolor=r.colour, elinewidth=1.0, capsize=1.6, zorder=3)
axB.set_xlim(0.795, 0.965)
axB.set_xticks([0.82, 0.92])
axB.tick_params(axis="x", labelsize=TS["small"] - 0.5)
axB.tick_params(axis="y", left=False, labelleft=False)
axB.grid(axis="y", visible=False)
axB.set_xlabel("AF3 confidence", fontsize=TS["small"])
for b in blocks:                        # covariate + omnibus test, in the gap above
    exp = int(np.floor(np.log10(b["kw"])))
    axB.text(0.800, ypos[b["i0"]] + 0.55,
             f"{b['ccol']}\nKW $p$<10$^{{{exp + 1}}}$", ha="left", va="center",
             fontsize=TS["small"] - 1.5, color=b["colour"], linespacing=1.1)
panel(axB, "B", dx=-0.62)

# --- C: contrast concordance ---
thr = -np.log10(0.05)
hi = 17.6
axC.axhspan(thr, hi, xmin=0, xmax=1, color=P.FILL_LT, alpha=0.55, zorder=0)
axC.plot([0, hi], [0, hi], ls=":", lw=0.8, color=P.GRID, zorder=1)
axC.axhline(thr, ls="--", lw=0.8, color=P.INK, zorder=2)
axC.axvline(thr, ls="--", lw=0.8, color=P.INK, zorder=2)
for key, _, colr, *_ in SERIES:
    t = tests[tests.series == key]
    keep, rev = t[t.direction_preserved], t[~t.direction_preserved]
    axC.scatter(keep.lx, keep.ly, s=26, color=colr, edgecolor=P.INK, lw=0.4,
                alpha=0.9, zorder=4, label=key)
    axC.scatter(rev.lx, rev.ly, s=30, marker="X", color=P.WHITE, edgecolor=colr,
                lw=1.1, zorder=4)
axC.set_xlim(-0.6, hi)
axC.set_ylim(-0.6, hi)
axC.set_xlabel("$-$log$_{10}$ $p$.adj  (raw)", fontsize=TS["small"])
axC.set_ylabel("$-$log$_{10}$ $p$.adj  (partialled)", fontsize=TS["small"])
axC.tick_params(labelsize=TS["small"] - 0.5)
axC.text(3.8, thr - 0.55, "lost", ha="left", va="top",
         fontsize=TS["small"] - 1, color=P.RED_D, style="italic")
axC.text(0.62, hi - 0.4, "gained", ha="center", va="top", rotation=90,
         fontsize=TS["small"] - 1, color=P.GREEN_D, style="italic")
axC.legend(handles=[Line2D([], [], marker="o", ls="none", ms=5.0, mec=P.INK, mew=0.4,
                           color=c, label=k) for k, _, c, *_ in SERIES]
                   + [Line2D([], [], marker="X", ls="none", ms=5.4, mfc=P.WHITE,
                             mec=P.INK, mew=1.1, label="dir. reversed")],
           loc="lower right", fontsize=TS["small"] - 1, handletextpad=0.1,
           borderpad=0.25, labelspacing=0.18, frameon=True, framealpha=1.0,
           facecolor=P.WHITE, edgecolor="none")
panel(axC, "C", f"Contrast concordance ($n$={len(tests)})", dx=-0.28)

# --- D: outcome tally ---
CATS = [("kept", P.INK2), ("newly sig.", P.GREEN_D), ("lost", P.RED_D),
        ("n.s. both", P.GRID)]
tally = []
for key, *_ in [(s[0],) for s in SERIES]:
    t = tests[tests.series == key]
    tally.append([int((t.sig_raw & t.surv).sum()), int((~t.sig_raw & t.surv).sum()),
                  int((t.sig_raw & ~t.surv).sum()), int((~t.sig_raw & ~t.surv).sum())])
tally = np.array(tally)
yb = np.arange(len(SERIES))[::-1]
left = np.zeros(len(SERIES))
for j, (cat, colr) in enumerate(CATS):
    axD.barh(yb, tally[:, j], left=left, height=0.6, color=colr,
             edgecolor=P.INK, lw=0.5, label=cat, zorder=3)
    for i, v in enumerate(tally[:, j]):
        if v:
            axD.text(left[i] + v / 2, yb[i], str(v), ha="center", va="center",
                     fontsize=TS["small"] - 0.5, color=P.WHITE if j < 3 else P.INK,
                     fontweight="bold", zorder=4)
    left += tally[:, j]
axD.set_yticks(yb)
axD.set_yticklabels([s[0] for s in SERIES], fontsize=TS["small"])
for tick, b in zip(axD.get_yticklabels(), blocks):
    tick.set_color(b["colour"])
    tick.set_fontweight("bold")
axD.set_xlabel("pairwise contrasts", fontsize=TS["small"])
axD.set_xlim(0, 22)
axD.set_ylim(-0.55, len(SERIES) - 1 + 1.75)   # blank band at the top for the key
axD.tick_params(axis="x", labelsize=TS["small"] - 0.5)
axD.grid(axis="y", visible=False)
axD.legend(loc="upper center", bbox_to_anchor=(0.5, 1.02), frameon=False,
           fontsize=TS["small"] - 1.5, ncol=2, columnspacing=0.6, handlelength=0.9,
           handletextpad=0.3, labelspacing=0.15, borderpad=0.0)
panel(axD, "D", "Verdict tally", dx=-0.28)

for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/FigureS5_confidence_control.{ext}",
                dpi=300 if ext == "png" else None)
print("saved FigureS5_confidence_control.{png,pdf}")

# ---- numbers for the legend ----------------------------------------------
print("\n--- panel A/B rows ---")
print(rows[["series", "label", "n", "raw", "part", "delta", "conf", "conf_sd"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
print("\n--- Kruskal-Wallis on confidence ---")
for b in blocks:
    print(f"  {b['key']}: {b['ccol']} p = {b['kw']:.3g}")
print("\n--- verdict tally (kept / newly sig / lost / n.s. both) ---")
for (key, *_), t in zip(SERIES, tally):
    print(f"  {key}: {list(t)}  (n={t.sum()}, sig_raw={t[0] + t[2]})")
print(f"  TOTAL contrasts = {tally.sum()}")
print("\n--- direction reversals ---")
for r in tests[~tests.direction_preserved].itertuples(index=False):
    print(f"  {r.series}: {flat(r.group_1)} vs {flat(r.group_2)}  "
          f"raw={r.fdr_raw:.3g} part={r.fdr_partial:.3g}")
print("\n--- headline contrasts ---")
for r in tests.itertuples(index=False):
    g1, g2 = flat(r.group_1), flat(r.group_2)
    if (r.series == "CSK" and {design(r.group_1), design(r.group_2)}
            == {"csk-holo/src-apo", "csk-holo/src-py159-holo"}) or \
       (r.series == "CDK1" and {g1, g2} == {"cdk1 [CDK1 = Apo]", "cdk1 ccnb1 [CDK1 = Apo]"}):
        print(f"  {r.series}: {g1} vs {g2}\n     MAC {r.MAC_raw_1:.3f}/{r.MAC_raw_2:.3f} -> "
              f"{r.MAC_partial_1:.3f}/{r.MAC_partial_2:.3f}  "
              f"p.adj {r.fdr_raw:.2g} -> {r.fdr_partial:.2g}")
