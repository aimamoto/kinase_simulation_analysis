#!/usr/bin/env python3
"""
mac_confidence_control.py -- is MAC a model-confidence artefact?

Recomputes the Global Network Rigidity (MAC) of the AlloQuant pipeline while
partialling per-model AlphaFold3 confidence (ipTM / pTM / interface PAE) out of every
edge of the Spearman distance-correlation matrix, and re-tests the between-condition
contrasts on the partialled edges.

Rationale
---------
MAC is the mean absolute Spearman correlation across the structural-metric network of an
ensemble. A reviewer can reasonably object that condition-to-condition MAC differences
reflect differences in AlphaFold3 confidence (better-predicted ensembles being more
internally self-consistent) rather than genuine mechanical coupling. This script tests
that directly by recomputing every edge as a partial Spearman correlation conditioned on
the confidence covariates, then re-running the same Mann-Whitney contrast the pipeline
uses on the resulting edge distributions.

Method
------
Edges use the closed-form partial-correlation identity, applied sequentially, one
covariate at a time -- the same approach already used in addons/paired_chain_coupling.py
for the cross-chain coupling analysis (Supplementary Note 10):

    rho(x,y | z) = (r_xy - r_xz * r_yz) / sqrt((1 - r_xz^2) * (1 - r_yz^2))

To stay comparable with the published MAC values the script mirrors the R engine
(modules/multimer_core_engine.R, Phase 5 / Phase 6b):
  * only "*_Dist" columns enter the network (angles break on periodicity),
  * within a group a column is used if it has >= MIN_NONNA non-NA values,
  * correlations that cannot be computed are set to 0,
  * MAC = mean of |rho| over the upper triangle,
  * conditions are contrasted by Mann-Whitney on the edge distributions.

IMPORTANT interpretation note
-----------------------------
Partialling removes variance, so |rho| -- and therefore MAC -- falls for EVERY group.
That is arithmetic, not evidence of weakening. The claim under test is that the
*between-condition contrast* survives: direction preserved and still significant.

Monomeric conditions have no interface, hence no ipTM and no PAE. For those the script
falls back to whichever covariates are populated (typically pTM alone) and records what
was actually used per group in the `covars_used` column. Read that column before
interpreting any monomer-vs-dimer contrast (e.g. CDK1 apo monomer vs CDK1-CCNB1).

Usage
-----
  # CSK-SRC (confidence lives in the 260718 rerun directory)
  python3 addons/mac_confidence_control.py \
      --metadata plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv \
      --confidence /media/akira/argentee/structure/SRC/af3_output/\
260718_rerun_260606_v7_csk-src/master_kinase_analysis_results_v7r2.csv \
      --label CSK --out-dir addons/mac_confidence_out

  # intrinsic (per-state) MAC instead of per-condition
  ... --group-col Macro_State --label CSK_states

  # CDK1-CCNB1 (confidence is already in its own master table)
  python3 addons/mac_confidence_control.py \
      --metadata <CDK1 Phase7 metadata csv> \
      --confidence /media/akira/argentee/structure/CDK1/af3_output/\
260517_CDK1-CCNB1/master_kinase_analysis_results_v7r2.csv \
      --label CDK1 --out-dir addons/mac_confidence_out
"""

import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

MIN_NONNA = 5      # matches the R engine's viability rule for a distance column
EPS = 1e-12


# ----------------------------------------------------------------------------- helpers
def bh_fdr(p):
    """Benjamini-Hochberg. NaNs are carried through as NaN."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order] * len(q) / (np.arange(len(q)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty_like(ranked)
    adj[order] = np.clip(ranked, 0, 1)
    out[ok] = adj
    return out


def spearman(a, b):
    """Spearman rho on complete pairs; NaN if under-determined or degenerate."""
    m = ~(pd.isna(a) | pd.isna(b))
    if m.sum() < 3:
        return np.nan
    x, y = a[m], b[m]
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    r, _ = stats.spearmanr(x, y)
    return r if np.isfinite(r) else np.nan


def partial_spearman(a, b, covars):
    """Sequential closed-form partial Spearman of a,b given a list of covariate arrays.

    All series are first restricted to rows complete in a, b and every covariate, so each
    correlation in the identity is computed on the same sample -- required for the
    closed form to be valid.
    """
    mask = ~(pd.isna(a) | pd.isna(b))
    for z in covars:
        mask &= ~pd.isna(z)
    if mask.sum() < 3 + len(covars):
        return np.nan
    x, y = a[mask], b[mask]
    r = spearman(x, y)
    if np.isnan(r):
        return np.nan
    for z in covars:
        zz = z[mask]
        rxz, ryz = spearman(x, zz), spearman(y, zz)
        if np.isnan(rxz) or np.isnan(ryz):
            continue                      # covariate is constant here: nothing to remove
        denom = np.sqrt(max((1 - rxz**2) * (1 - ryz**2), 0.0))
        if denom < EPS:
            return np.nan
        r = (r - rxz * ryz) / denom
        r = float(np.clip(r, -1.0, 1.0))
    return r


def edge_vectors(sub, dist_cols, covars_present):
    """Return (raw_edges, partial_edges, cols_used) for one group.

    Mirrors the R engine: unusable correlations become 0 so that MAC stays comparable
    with the published Phase 5 / Phase 6b values.
    """
    cols = [c for c in dist_cols if sub[c].notna().sum() >= MIN_NONNA]
    if len(cols) < 2:
        return None, None, cols
    cov_arrays = [sub[c].to_numpy(dtype=float) for c in covars_present]
    raw, par = [], []
    for i, j in itertools.combinations(range(len(cols)), 2):
        a = sub[cols[i]].to_numpy(dtype=float)
        b = sub[cols[j]].to_numpy(dtype=float)
        r = spearman(a, b)
        p = partial_spearman(a, b, cov_arrays) if cov_arrays else r
        raw.append(0.0 if np.isnan(r) else abs(r))
        par.append(0.0 if np.isnan(p) else abs(p))
    return np.array(raw), np.array(par), cols


# -------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True,
                    help="Phase7_Complete_Structural_Metadata.csv (per-model rows, has Condition)")
    ap.add_argument("--confidence", required=True,
                    help="master CSV carrying populated ipTM / pTM / PAE columns")
    ap.add_argument("--join-key", default="Simulation_ID")
    ap.add_argument("--group-col", default="Condition",
                    help="Condition (Phase 5 global MAC) or Macro_State (Phase 6b intrinsic MAC)")
    ap.add_argument("--covars", nargs="*", default=["ipTM", "PAE_Mean_AB"],
                    help="confidence covariates to partial out")
    ap.add_argument("--covar-mode", choices=["common", "per-group"], default="common",
                    help="'common' (default, CORRECT for between-group contrasts): partial out only "
                         "the covariates populated in EVERY group, so all groups lose the same "
                         "variance. 'per-group': use whatever is available in each group -- this "
                         "makes MAC_partial NON-COMPARABLE across groups and must not be used for "
                         "contrasts. Matters whenever a monomeric condition (no interface, hence no "
                         "ipTM/PAE) is compared with a complex.")
    ap.add_argument("--label", default="run")
    ap.add_argument("--out-dir", default="addons/mac_confidence_out")
    ap.add_argument("--phase5", default=None,
                    help="optional Phase5_Global_Network_Density.csv, to verify MAC_raw reproduces the pipeline")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    meta = pd.read_csv(args.metadata)
    conf = pd.read_csv(args.confidence)

    if args.group_col not in meta.columns:
        sys.exit(f"[error] group column '{args.group_col}' not in metadata")

    # ---- attach confidence -------------------------------------------------------
    cov_cols = [c for c in args.covars if c in conf.columns]
    missing = [c for c in args.covars if c not in conf.columns]
    if missing:
        print(f"[warn] covariate(s) absent from the confidence table: {missing}")
    if not cov_cols:
        sys.exit("[error] none of the requested covariates exist in the confidence table.")

    conf_small = conf[[args.join_key] + cov_cols].drop_duplicates(subset=[args.join_key])
    n_before = len(meta)
    df = meta.merge(conf_small, on=args.join_key, how="left", suffixes=("", "_conf"))
    # if metadata already carried (empty) covariate columns, prefer the joined values
    for c in cov_cols:
        if f"{c}_conf" in df.columns:
            df[c] = df[f"{c}_conf"]
    matched = df[cov_cols].notna().any(axis=1).sum()
    print(f"[*] {n_before} metadata rows; confidence attached to {matched} "
          f"({matched / n_before * 100:.1f}%)")
    if matched == 0:
        sys.exit("[error] join produced no confidence values -- check --join-key and --confidence.")

    dist_cols = [c for c in df.columns if c.endswith("_Dist")]
    print(f"[*] {len(dist_cols)} distance metrics; covariates: {cov_cols}")

    # ---- A. precondition: does confidence differ between groups at all? ----------
    rows = []
    for g, sub in df.groupby(args.group_col):
        rec = {args.group_col: g, "N": len(sub)}
        for c in cov_cols:
            rec[f"{c}_mean"] = sub[c].mean()
            rec[f"{c}_sd"] = sub[c].std()
        rows.append(rec)
    conf_by_group = pd.DataFrame(rows)
    for c in cov_cols:
        samples = [s[c].dropna().values for _, s in df.groupby(args.group_col)]
        samples = [s for s in samples if len(s) > 1]
        if len(samples) > 1:
            try:
                h, p = stats.kruskal(*samples)
                print(f"[*] Kruskal-Wallis across groups on {c}: H={h:.1f}, p={p:.3g}")
            except ValueError:
                pass
    conf_by_group.to_csv(
        os.path.join(args.out_dir, f"{args.label}_confidence_by_group.csv"), index=False)

    # ---- B. raw vs partial MAC ---------------------------------------------------
    # Establish the covariate set BEFORE looping. Under the default 'common' mode every
    # group is partialled on the identical covariates, which is what makes MAC_partial
    # comparable between groups. Falling back per-group would have e.g. the CDK1 apo
    # monomer lose only pTM while CDK1-CCNB1 loses ipTM+pTM+PAE -- different amounts of
    # removed variance, so the resulting contrast would be meaningless.
    per_group_avail = {
        g: [c for c in cov_cols if sub[c].notna().sum() >= MIN_NONNA]
        for g, sub in df.groupby(args.group_col)
    }
    if args.covar_mode == "common":
        common = [c for c in cov_cols if all(c in v for v in per_group_avail.values())]
        dropped = [c for c in cov_cols if c not in common]
        if dropped:
            lacking = {g: sorted(set(dropped) - set(v)) for g, v in per_group_avail.items()
                       if set(dropped) - set(v)}
            print(f"[*] covar-mode=common: dropping {dropped} -- not populated in every group "
                  f"(e.g. {list(lacking)[:2]}). Partialling on {common or 'NOTHING'}.")
        if not common:
            sys.exit("[error] no covariate is populated in every group. Restrict --group-col, "
                     "or pass --covars pTM, or accept --covar-mode per-group (contrasts invalid).")
    else:
        common = None
        print("[warn] covar-mode=per-group: MAC_partial is NOT comparable between groups; "
              "pairwise contrasts below are unreliable.")

    store, mac_rows = {}, []
    for g, sub in df.groupby(args.group_col):
        present = common if common is not None else per_group_avail[g]
        raw, par, cols = edge_vectors(sub, dist_cols, present)
        if raw is None:
            print(f"[warn] group '{g}': fewer than 2 viable metrics, skipped")
            continue
        store[g] = (raw, par)
        mac_rows.append({
            args.group_col: g, "N": len(sub), "n_metrics": len(cols), "n_edges": len(raw),
            "MAC_raw": raw.mean(), "MAC_partial": par.mean(),
            "delta": par.mean() - raw.mean(),
            "covars_used": ";".join(present) if present else "NONE",
        })
    mac = pd.DataFrame(mac_rows).sort_values("MAC_raw", ascending=False)
    mac.to_csv(os.path.join(args.out_dir, f"{args.label}_mac_raw_vs_partial.csv"), index=False)

    print("\n=== MAC: raw vs confidence-partialled ===")
    print(mac.to_string(index=False,
                        float_format=lambda v: f"{v:.4f}" if isinstance(v, float) else v))

    # optional check that MAC_raw reproduces the published pipeline value
    if args.phase5 and os.path.exists(args.phase5):
        p5 = pd.read_csv(args.phase5)
        key = "Condition" if "Condition" in p5.columns else p5.columns[0]
        chk = mac.merge(p5[[key, "Global_Coupling_Score"]],
                        left_on=args.group_col, right_on=key, how="inner")
        if len(chk):
            chk["abs_diff"] = (chk.MAC_raw - chk.Global_Coupling_Score).abs()
            print(f"\n[*] reproduction check vs Phase5: max |MAC_raw - published| = "
                  f"{chk.abs_diff.max():.5f} over {len(chk)} groups")

    # ---- C. pairwise contrasts on raw and partial edges --------------------------
    tests = []
    for g1, g2 in itertools.combinations(store.keys(), 2):
        r1, p1 = store[g1]
        r2, p2 = store[g2]
        try:
            _, praw = stats.mannwhitneyu(r1, r2, alternative="two-sided")
        except ValueError:
            praw = np.nan
        try:
            _, ppar = stats.mannwhitneyu(p1, p2, alternative="two-sided")
        except ValueError:
            ppar = np.nan
        tests.append({
            "group_1": g1, "group_2": g2,
            "MAC_raw_1": r1.mean(), "MAC_raw_2": r2.mean(),
            "MAC_partial_1": p1.mean(), "MAC_partial_2": p2.mean(),
            "direction_preserved": bool(
                np.sign(r1.mean() - r2.mean()) == np.sign(p1.mean() - p2.mean())),
            "p_raw": praw, "p_partial": ppar,
        })
    tst = pd.DataFrame(tests)
    if len(tst):
        tst["fdr_raw"] = bh_fdr(tst.p_raw)
        tst["fdr_partial"] = bh_fdr(tst.p_partial)
        tst["survives"] = (tst.direction_preserved) & (tst.fdr_partial < 0.05)
        tst = tst.sort_values("fdr_partial")
        tst.to_csv(os.path.join(args.out_dir, f"{args.label}_mac_pairwise_tests.csv"),
                   index=False)
        n_sig_raw = int((tst.fdr_raw < 0.05).sum())
        n_surv = int(tst.survives.sum())
        print(f"\n=== contrasts: {n_sig_raw} significant raw (FDR<0.05); "
              f"{n_surv} survive partialling with direction preserved ===")
        cols = ["group_1", "group_2", "MAC_raw_1", "MAC_raw_2",
                "MAC_partial_1", "MAC_partial_2", "fdr_raw", "fdr_partial", "survives"]
        print(tst[cols].head(15).to_string(index=False,
              float_format=lambda v: f"{v:.4g}" if isinstance(v, float) else v))

    # ---- figure ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        m = mac.copy()
        lbl = [str(s).replace("\\n", " / ").replace("\n", " / ")[:44] for s in m[args.group_col]]
        y = np.arange(len(m))[::-1]
        fig, ax = plt.subplots(figsize=(9, 0.45 * len(m) + 2.2))
        ax.barh(y + 0.19, m.MAC_raw, height=0.36, label="MAC (raw)", color="#3498db")
        ax.barh(y - 0.19, m.MAC_partial, height=0.36,
                label=f"MAC | {', '.join(cov_cols)}", color="#e67e22")
        ax.set_yticks(y)
        ax.set_yticklabels(lbl, fontsize=8)
        ax.set_xlabel("Mean Absolute Correlation (MAC)")
        ax.set_title(f"{args.label}: network rigidity before and after partialling "
                     f"out AlphaFold3 confidence", fontsize=11)
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(axis="x", alpha=0.3, linestyle=":")
        fig.text(0.01, 0.01,
                 "Partialling removes variance, so MAC falls for every group; the test is "
                 "whether the between-condition contrast survives.",
                 fontsize=7, style="italic", color="#555555")
        fig.tight_layout(rect=[0, 0.035, 1, 1])
        out_png = os.path.join(args.out_dir, f"{args.label}_mac_raw_vs_partial.png")
        fig.savefig(out_png, dpi=300)
        print(f"\n[*] figure written: {out_png}")
    except ImportError:
        print("[warn] matplotlib unavailable; figure skipped")

    print(f"[*] outputs in {args.out_dir}/")


if __name__ == "__main__":
    main()
