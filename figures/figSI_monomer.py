#!/usr/bin/env python
"""Figure S6 — the active CSK conformer does not pre-exist in the unbound kinase.

Control for the conformational-selection framing of Model 2. Companion to the new
Supplementary Note (monomer control); register row A8.

A: alphaC-In occupancy, nucleotide-matched (every ensemble is CSK + ATP): unbound CSK,
   docked/unprimed, docked/primed. Exact binomial 95% CIs.
B: why the unbound n = 200 is real -- pairwise Calpha RMSD shows the holo arm is a genuine
   ensemble while the apo arm collapses to a single converged structure.
The across-build comparison (unbound apo answers 100% / 0.5% / 99.5% alphaC-In in three
AlphaFold3 runs of identical sequence, while the holo arm gives 0/25, 0/200, 1/200) was
built as a panel C and dropped 2026-08-06; it survives as one sentence of the Note. Set
SHOW_BUILDS = True to restore it.

Sources
  unbound monomer : 260805_CSK_40seeds_borealis (weights 2524e951..., the bundle used for
                    every production run in this study), plus 260324_CSK and
                    260805_CSK_40seeds_cell for the build comparison in C only.
  docked complex  : 260718_rerun_260606_v7_csk-src Phase7_Complete_Structural_Metadata.csv,
                    grouped by Condition_reviewed (never by directory name).
"""
import os
import sys
import glob

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import palette as P
P.apply_rc(matplotlib)
TS = P.TS

SRCBASE = os.environ.get("ALLOQUANT_AF3_ROOT", "af3_output")  # parent dir holding the AF3 run directories
BOREALIS = f"{SRCBASE}/260805_CSK_40seeds_borealis"
CELL = f"{SRCBASE}/260805_CSK_40seeds_cell"
MARCH = f"{SRCBASE}/260324_CSK/af_analysis_ready"
DOCKED = (f"{SRCBASE}/260718_rerun_260606_v7_csk-src/plots_and_stats_CSK_GMM/"
          "Phase7_Complete_Structural_Metadata.csv")
CACHE = os.path.join(HERE, ".figSI_monomer_rmsd.npz")
OUT = os.environ.get("ALLOQUANT_FIGDIR", ".")

MASTER = "master_kinase_analysis_results_v7r3.csv"
SHOW_BUILDS = False   # panel C, dropped 2026-08-06


# ---- alphaC-In counts -----------------------------------------------------
def monomer_counts(run_dir):
    """{'apo': (n_In, n_total), 'holo': (...)} for a monomer run directory."""
    d = pd.read_csv(os.path.join(run_dir, MASTER))
    d["cond"] = d["Directory"].str.split("/").str[0]
    out = {}
    for cond, key in (("csk-wtcat_apo", "apo"), ("csk-wtcat_holo", "holo")):
        s = d[d["cond"] == cond]
        out[key] = (int((s["C_Helix"] == "In").sum()), len(s))
    return out


def docked_counts():
    d = pd.read_csv(DOCKED)
    d["design"] = (d["Condition_reviewed"].astype(str)
                   .str.replace("wtcat-", "", regex=False)
                   .str.replace("\n", "/", regex=False)
                   .str.split("[").str[0].str.strip().str.rstrip("/"))
    out = {}
    for design in ("csk-holo/src-holo", "csk-holo/src-py159-holo"):
        s = d[d["design"] == design]
        out[design] = (int((s["C_Helix"] == "In").sum()), len(s))
    return out


def ci(k, n):
    lo, hi = stats.binomtest(k, n).proportion_ci(method="exact")
    return 100 * k / n, 100 * lo, 100 * hi


# ---- pairwise Calpha RMSD -------------------------------------------------
def read_ca(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    cols = [l.strip()[len("_atom_site."):] for l in lines if l.startswith("_atom_site.")]
    ia, ix = cols.index("label_atom_id"), cols.index("Cartn_x")
    xyz = [(float(p[ix]), float(p[ix + 1]), float(p[ix + 2]))
           for p in (l.split() for l in lines if l.startswith("ATOM "))
           if p[ia] == "CA"]
    return np.asarray(xyz, dtype=np.float64)


def pairwise_rmsd(run_dir, cond, chunk=2000):
    """All pairwise best-fit Calpha RMSDs (frame-independent), closed-form Kabsch."""
    files = sorted(glob.glob(os.path.join(run_dir, cond, "*", "model.cif")))
    X = np.stack([read_ca(f) for f in files])           # (N, L, 3)
    X -= X.mean(axis=1, keepdims=True)
    n, L, _ = X.shape
    sq = (X ** 2).sum(axis=(1, 2))                       # per-structure |X|^2
    i, j = np.triu_indices(n, k=1)
    out = np.empty(len(i))
    for s in range(0, len(i), chunk):
        a, b = i[s:s + chunk], j[s:s + chunk]
        H = np.einsum("pki,pkj->pij", X[a], X[b])        # (p, 3, 3)
        sv = np.linalg.svd(H, compute_uv=False)
        d = np.sign(np.linalg.det(H))
        tr = sv[:, 0] + sv[:, 1] + d * sv[:, 2]          # sign-corrected trace
        out[s:s + chunk] = np.sqrt(np.maximum(sq[a] + sq[b] - 2 * tr, 0) / L)
    return out


def get_rmsd():
    if os.path.exists(CACHE):
        z = np.load(CACHE)
        return {k: z[k] for k in z.files}
    print("[*] computing pairwise RMSD (cached after first run)...")
    d = {"apo": pairwise_rmsd(BOREALIS, "csk-wtcat_apo"),
         "holo": pairwise_rmsd(BOREALIS, "csk-wtcat_holo")}
    np.savez_compressed(CACHE, **d)
    return d


# ---- assemble -------------------------------------------------------------
bor, cel, mar = monomer_counts(BOREALIS), monomer_counts(CELL), monomer_counts(MARCH)
dock = docked_counts()
rmsd = get_rmsd()

UNPRIMED, PRIMED = "csk-holo/src-holo", "csk-holo/src-py159-holo"
bars = [("unbound\nCSK", *bor["holo"]),
        ("docked\nunprimed", *dock[UNPRIMED]),
        ("docked\nprimed", *dock[PRIMED])]

ft = stats.fisher_exact([[bor["holo"][0], bor["holo"][1] - bor["holo"][0]],
                         [dock[UNPRIMED][0], dock[UNPRIMED][1] - dock[UNPRIMED][0]]])

print("\n=== numbers used in Fig. S6 ===")
for lab, k, n in bars:
    m, lo, hi = ci(k, n)
    print(f"  {lab.replace(chr(10), ' '):38s} {k:3d}/{n:3d} = {m:5.1f}%  95% CI [{lo:.2f}, {hi:.2f}]")
print(f"  Fisher unbound vs docked-unprimed: p = {ft.pvalue:.4f}")
for nm, c in (("March 260324", mar), ("borealis", bor), ("cell", cel)):
    print(f"  {nm:14s} apo {c['apo'][0]:3d}/{c['apo'][1]:3d}   holo {c['holo'][0]:3d}/{c['holo'][1]:3d}")
print(f"  pairwise RMSD  apo  mean {rmsd['apo'].mean():.2f} A (max {rmsd['apo'].max():.2f})")
print(f"  pairwise RMSD  holo mean {rmsd['holo'].mean():.2f} A (max {rmsd['holo'].max():.2f})")

ncol = 3 if SHOW_BUILDS else 2
fig = plt.figure(figsize=(7.4 if SHOW_BUILDS else 5.1, 3.0))
gs = fig.add_gridspec(1, ncol,
                      width_ratios=[1.02, 1.00, 1.06][:ncol],
                      wspace=0.46 if SHOW_BUILDS else 0.40,
                      left=0.085 if SHOW_BUILDS else 0.115,
                      right=0.985, top=0.84, bottom=0.28)

# --- A: occupancy ----------------------------------------------------------
axA = fig.add_subplot(gs[0, 0])
xs = np.arange(len(bars))
cols = [P.MUTED, P.BLUE, P.GREEN_D]
for x, (lab, k, n), c in zip(xs, bars, cols):
    m, lo, hi = ci(k, n)
    axA.bar(x, m, width=0.62, color=c, edgecolor=P.INK, linewidth=0.7, zorder=3)
    axA.errorbar(x, m, yerr=[[m - lo], [hi - m]], fmt="none", ecolor=P.INK,
                 elinewidth=0.9, capsize=3, zorder=4)
    axA.text(x, hi + 1.4, f"{k}/{n}", ha="center", va="bottom",
             fontsize=TS["value"], color=P.INK, zorder=5)
axA.set_xticks(xs)
axA.set_xticklabels([b[0] for b in bars], fontsize=TS["small"])
axA.set_ylabel("CSK αC-In (% of models)")
axA.set_ylim(0, 42)
axA.set_title("Nucleotide-matched (all CSK + ATP)", fontsize=TS["tick"],
              color=P.MUTED, pad=3)
# significance bracket, unbound vs docked-unprimed
y = 12.5
axA.plot([0, 0, 1, 1], [y, y + 1.1, y + 1.1, y], lw=0.8, color=P.INK, zorder=5)
axA.text(0.5, y + 1.7, f"$p$ = {ft.pvalue:.3f}", ha="center", va="bottom",
         fontsize=TS["anno"], color=P.INK)
axA.text(-0.02, 1.14, "A", transform=axA.transAxes, fontsize=TS["panel"],
         fontweight="bold", va="top", ha="right")

# --- B: is it an ensemble? -------------------------------------------------
axB = fig.add_subplot(gs[0, 1])
data = [rmsd["holo"], rmsd["apo"]]
labs = ["+ ATP\n(holo)", "apo"]
bp = axB.boxplot(data, vert=True, widths=0.55, showfliers=False, patch_artist=True,
                 medianprops=dict(color=P.INK, lw=1.1),
                 boxprops=dict(lw=0.8, edgecolor=P.INK),
                 whiskerprops=dict(lw=0.8, color=P.INK),
                 capprops=dict(lw=0.8, color=P.INK))
for patch, c in zip(bp["boxes"], [P.BLUE, P.FILL_LT]):
    patch.set_facecolor(c)
axB.set_xticklabels(labs, fontsize=TS["small"])
axB.set_ylabel("pairwise Cα RMSD (Å)")
top = max(rmsd["holo"].max(), rmsd["apo"].max())
axB.set_ylim(0, top * 1.40)
axB.set_title("Unbound CSK, n = 200 each", fontsize=TS["tick"], color=P.MUTED, pad=3)
# direct labels above each box rather than leader lines
for xpos, arr, note in ((1, rmsd["holo"], "genuine\nensemble"),
                        (2, rmsd["apo"], "one structure")):
    axB.text(xpos, arr.max() + top * 0.06, f"{arr.mean():.2f} Å\n{note}", ha="center",
             va="bottom", fontsize=TS["anno"], color=P.INK, linespacing=1.25)
axB.text(-0.02, 1.14, "B", transform=axB.transAxes, fontsize=TS["panel"],
         fontweight="bold", va="top", ha="right")

# --- C: build dependence (retired; SHOW_BUILDS) -----------------------------
if SHOW_BUILDS:
  axC = fig.add_subplot(gs[0, 2])
  builds = [("2026-03\nparams A", mar), ("2026-08\nparams A", bor), ("2026-08\nparams B", cel)]
  w = 0.36
  xs = np.arange(len(builds))
  for off, key, c, lab in ((-w / 2, "apo", P.FILL_LT, "apo"),
                           (+w / 2, "holo", P.BLUE, "+ ATP (holo)")):
      vals = [100 * b[key][0] / b[key][1] for _, b in builds]
      axC.bar(xs + off, vals, width=w, color=c, edgecolor=P.INK, linewidth=0.7,
              label=lab, zorder=3)
      for x, v, (_, b) in zip(xs + off, vals, builds):
          axC.text(x, v + 3, f"{b[key][0]}/{b[key][1]}", ha="center", va="bottom",
                   fontsize=TS["small"], color=P.INK, rotation=90, zorder=5)
  axC.set_xticks(xs)
  axC.set_xticklabels([b[0] for b in builds], fontsize=TS["small"])
  axC.set_ylabel("CSK αC-In (% of models)")
  axC.set_ylim(0, 158)
  axC.set_yticks([0, 25, 50, 75, 100])
  axC.set_title("Unbound CSK across AlphaFold3 builds", fontsize=TS["tick"],
                color=P.MUTED, pad=3)
  axC.legend(frameon=False, fontsize=TS["small"], loc="upper center",
             bbox_to_anchor=(0.5, 1.005), ncol=2, handlelength=1.0, columnspacing=0.8,
             handletextpad=0.4)
  # this study's production runs all use params A
  axC.annotate("this study", xy=(1, 0), xytext=(1, -0.205), textcoords=("data", "axes fraction"),
               ha="center", va="top", fontsize=TS["anno"], color=P.MUTED)
  axC.text(-0.02, 1.14, "C", transform=axC.transAxes, fontsize=TS["panel"],
           fontweight="bold", va="top", ha="right")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"figS6_monomer.{ext}"), dpi=300,
                bbox_inches="tight", facecolor="white")
print(f"\n[*] wrote {OUT}/figS6_monomer.png / .pdf")
