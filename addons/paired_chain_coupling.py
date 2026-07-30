#!/usr/bin/env python
"""
paired_chain_coupling.py  —  COMPLEMENTARY ADD-ON to the kinase/allostery pipeline.

NOT part of the core pipeline (does not touch modules/ or the wrappers). This is a
post-hoc analysis that runs on the Module 2 outputs (Phase7 per-chain metric tables).

Concept
-------
Every multimeric AlphaFold3 model contains BOTH kinase domains in ONE structure, so the
per-model geometric measurements of the two chains are 1:1 paired. This tool tests whether
the two domains' conformations are coupled at the single-structure level — a question the
condition-grouped and GMM-state analyses cannot answer.

Method (matches the manuscript Methods)
---------------------------------------
1. Join the two per-chain Phase7 tables on the shared model id (Simulation_ID).
2. For each (chainA metric x chainB metric) pair, compute Spearman rho *within each
   experimental condition*, then combine across conditions by inverse-variance-weighted
   Fisher-z meta-analysis. The within-condition design removes correlation induced purely
   by the shared experimental grouping (the dominant confound).
3. Optionally recompute as PARTIAL Spearman controlling for per-model AF3 confidence
   (ipTM, mean interface PAE) via the closed-form partial-correlation identity — rules out
   a model-quality artifact.
4. FDR-correct (Benjamini-Hochberg) across the full metric grid.

Outputs (to --out-dir)
----------------------
  paired_coupling_full_grid.csv     all chainA x chainB pairs, raw + partial rho, p, FDR
  paired_coupling_homologous.csv    same-named metric pairs (the "conformational mirroring" set)

Usage
-----
  conda run -n main python addons/paired_chain_coupling.py \
      --chainA plots_and_stats_CSK_GMM/Phase7_Complete_Structural_Metadata.csv \
      --chainB plots_and_stats_SRC_GMM/Phase7_Complete_Structural_Metadata.csv \
      --out-dir addons/coupling_out

Requires the confidence columns (ipTM/PAE_*) to be populated for the --confidence control;
if they are all-NaN (sims predating the AF3-assessment upgrade) pass --confidence "" to skip.
"""
import argparse, os, sys
import numpy as np, pandas as pd
from scipy import stats


def spearman_n(a, b):
    d = pd.concat([a, b], axis=1).dropna()
    if len(d) < 7:
        return np.nan, len(d)
    return stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0], len(d)


def meta_raw(groups, x, y, min_n):
    zs, ws = [], []
    for g in groups:
        r, n = spearman_n(g[x], g[y])
        if np.isnan(r) or n < min_n:
            continue
        zs.append(np.arctanh(np.clip(r, -0.999, 0.999)))
        ws.append(n - 3)
    if not zs:
        return np.nan, np.nan, 0
    zs, ws = np.array(zs), np.array(ws)
    zb = (zs * ws).sum() / ws.sum()
    se = 1.0 / np.sqrt(ws.sum())
    return np.tanh(zb), 2 * (1 - stats.norm.cdf(abs(zb / se))), len(zs)


def meta_partial(groups, x, y, covars, min_n):
    """Meta-analyzed partial Spearman(x, y | covars) via closed form, one covar at a time
    applied sequentially (adequate for near-independent low-variance confidence covariates)."""
    zs, ws = [], []
    for g in groups:
        rxy, n = spearman_n(g[x], g[y])
        if np.isnan(rxy) or n < min_n:
            continue
        rp = rxy
        ok = True
        for z in covars:
            rxz, _ = spearman_n(g[x], g[z])
            ryz, _ = spearman_n(g[y], g[z])
            if np.isnan(rxz) or np.isnan(ryz):
                ok = False
                break
            den = np.sqrt((1 - rxz ** 2) * (1 - ryz ** 2))
            if den < 1e-9:
                ok = False
                break
            rp = (rp - rxz * ryz) / den
        if not ok:
            continue
        zs.append(np.arctanh(np.clip(rp, -0.999, 0.999)))
        ws.append(n - 3 - len(covars))
    if not zs:
        return np.nan, np.nan, 0
    zs, ws = np.array(zs), np.array(ws)
    zb = (zs * ws).sum() / ws.sum()
    se = 1.0 / np.sqrt(ws.sum())
    return np.tanh(zb), 2 * (1 - stats.norm.cdf(abs(zb / se))), len(zs)


def bh_fdr(pvals):
    p = np.asarray(pvals, float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    return out


def numeric_shared(a, b, key, exclude):
    cols = []
    for c in a.columns:
        if c in exclude or c == key or c not in b.columns:
            continue
        if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c]):
            if a[c].notna().sum() > 50 and a[c].nunique() > 8 and b[c].nunique() > 8:
                cols.append(c)
    return cols


def main():
    ap = argparse.ArgumentParser(description="Paired cross-chain conformational coupling (pipeline add-on).")
    ap.add_argument("--chainA", required=True, help="Phase7 metadata CSV for chain A (e.g. CSK/enzyme)")
    ap.add_argument("--chainB", required=True, help="Phase7 metadata CSV for chain B (e.g. SRC/substrate)")
    ap.add_argument("--key", default="Simulation_ID", help="shared per-model join key")
    ap.add_argument("--condition-col", default="Condition_reviewed", help="within-condition grouping column")
    ap.add_argument("--confidence", nargs="*", default=["ipTM", "PAE_Mean_AB"],
                    help="per-model confidence covariates to partial out (pass nothing / '' to skip)")
    ap.add_argument("--min-n", type=int, default=25, help="min models per condition to include")
    ap.add_argument("--out-dir", default="addons/coupling_out")
    args = ap.parse_args()

    A = pd.read_csv(args.chainA)
    B = pd.read_csv(args.chainB)
    if args.condition_col not in A.columns:
        sys.exit(f"[error] condition column '{args.condition_col}' not in chainA table.")
    covars = [c for c in (args.confidence or []) if c and c in A.columns and A[c].notna().sum() > 0]
    dropped = [c for c in (args.confidence or []) if c and c not in covars]
    if dropped:
        print(f"[warn] confidence covariate(s) unavailable/all-NaN, skipping: {dropped}")

    # model-level (per-structure, identical for both chains) columns must NOT enter the per-chain
    # metric grid — they trivially self-correlate at rho=1.0 and pollute FDR.
    MODEL_LEVEL = {"ipTM", "pTM", "PAE_A_to_B", "PAE_B_to_A", "PAE_Mean_AB", "PAE_Mean",
                   "Pct_High_Conf_Interface"}
    exclude = {args.condition_col, "Chain", "Type", "Role", "Directory", "File", "Partner",
               "CoFactor_Name", "Macro_State"} | set(covars) | MODEL_LEVEL
    mets = numeric_shared(A, B, args.key, exclude)
    if not mets:
        sys.exit("[error] no shared numeric metric columns found between the two tables.")
    print(f"[info] {len(mets)} shared metrics; joining on '{args.key}'.")

    aa = A[[args.key, args.condition_col] + mets + covars].rename(columns={m: "A_" + m for m in mets}) \
          .rename(columns={args.condition_col: "cond"}).drop_duplicates(args.key)
    bb = B[[args.key] + mets].rename(columns={m: "B_" + m for m in mets}).drop_duplicates(args.key)
    p = aa.merge(bb, on=args.key)
    groups = [g for _, g in p.groupby("cond") if len(g) >= args.min_n]
    print(f"[info] {len(p)} paired models; {len(groups)} conditions with n>={args.min_n}.")

    rows = []
    for ma in mets:
        for mb in mets:
            r_raw, p_raw, k = meta_raw(groups, "A_" + ma, "B_" + mb, args.min_n)
            if np.isnan(p_raw):
                continue
            if covars:
                r_par, p_par, _ = meta_partial(groups, "A_" + ma, "B_" + mb, covars, args.min_n)
            else:
                r_par, p_par = np.nan, np.nan
            rows.append(dict(chainA_metric=ma, chainB_metric=mb, homologous=(ma == mb),
                             rho_raw=r_raw, p_raw=p_raw, rho_partial=r_par, p_partial=p_par,
                             n_conditions=k))
    df = pd.DataFrame(rows)
    df["fdr_raw"] = bh_fdr(df["p_raw"])
    if covars:
        df["fdr_partial"] = bh_fdr(df["p_partial"].fillna(1.0))
    df = df.sort_values("p_raw").reset_index(drop=True)

    os.makedirs(args.out_dir, exist_ok=True)
    full = os.path.join(args.out_dir, "paired_coupling_full_grid.csv")
    homo = os.path.join(args.out_dir, "paired_coupling_homologous.csv")
    df.to_csv(full, index=False)
    df[df["homologous"]].sort_values("rho_raw", ascending=False).to_csv(homo, index=False)
    print(f"[done] wrote {full}\n[done] wrote {homo}")
    print("\nTop homologous (same-metric) couplings:")
    show = df[df["homologous"]].sort_values("rho_raw", ascending=False).head(10)
    cols = ["chainA_metric", "rho_raw", "p_raw"] + (["rho_partial", "p_partial"] if covars else [])
    print(show[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
