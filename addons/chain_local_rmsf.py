#!/usr/bin/env python3
"""Chain-local per-residue Ca RMSF across an AlphaFold3 ensemble.

Each chain is fitted on its own conformation in the FIRST model of the condition, using that
chain's Ca atoms only, so relative motion of one subunit against the other is excluded and the
value reports internal deformation alone. Per-residue RMSF is then taken about the chain's mean
position across the models of the condition, over residues present in every model.

This is deliberately NOT the shared-frame superposition used elsewhere for a whole complex: the
two kinds of value are not comparable, because a shared frame retains inter-chain motion.

Input:  an AlphaFold3 output directory holding one subdirectory per condition, each containing
        seed-*_sample-*/model.cif.
Output: one row per (condition, chain) with the median, 90th percentile, the fraction of
        residues above a mobility threshold, and the mean pLDDT; optionally the full
        per-residue profile.

Usage:
    python3 chain_local_rmsf.py --root "$ALLOQUANT_CDK1" \
        --conditions a-cdk1_b-ccnb1-166_0atp=+CCNB1 \
                     a-cdk1_b-ccnb1-166_1atp=+CCNB1+ATP \
                     a-cdk1-pt161_b-ccnb1-166_1atp="pT161 (active)" \
                     a-cdk1-pt14-py15-pt161_b-ccnb1-166_1atp=triple-P \
        --chains A=CDK1 B=CCNB1 --out chain_local_rmsf.csv
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

CA = "CA"


def read_ca(path):
    """Return {chain: (seq_ids, coords Nx3, plddt N)} for Ca atoms of polymer residues.

    Parses the _atom_site loop positionally after reading its column order from the header, so
    it does not assume AlphaFold3's current column layout.
    """
    cols, rows, in_loop = [], [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".", 1)[1])
                in_loop = True
                continue
            if in_loop:
                if line.startswith(("ATOM", "HETATM")):
                    rows.append(line.split())
                elif rows:
                    break
    if not cols or not rows:
        raise ValueError(f"no _atom_site records in {path}")
    ix = {name: i for i, name in enumerate(cols)}
    need = ("label_atom_id", "auth_asym_id", "label_seq_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv")
    missing = [n for n in need if n not in ix]
    if missing:
        raise ValueError(f"{path}: _atom_site lacks {missing}")

    out = {}
    for f in rows:
        if f[ix["label_atom_id"]] != CA:
            continue
        seq = f[ix["label_seq_id"]]
        if seq in (".", "?"):
            continue
        # Polymer membership is decided by label_seq_id, not by group_PDB: AlphaFold3 writes
        # modified residues (TPO, PTR, SEP) as HETATM although they are part of the chain,
        # while true ligands (ATP, MG) carry no label_seq_id. Filtering on group_PDB would
        # silently drop the phosphosites and shorten the phosphorylated chains.
        ch = f[ix["auth_asym_id"]]
        rec = out.setdefault(ch, ([], [], []))
        rec[0].append(int(seq))
        rec[1].append((float(f[ix["Cartn_x"]]), float(f[ix["Cartn_y"]]), float(f[ix["Cartn_z"]])))
        rec[2].append(float(f[ix["B_iso_or_equiv"]]))
    return {c: (np.array(s), np.array(xyz, float), np.array(b, float))
            for c, (s, xyz, b) in out.items()}


def kabsch_fit(mobile, target):
    """Least-squares superposition of mobile onto target; returns the transformed mobile."""
    mc, tc = mobile.mean(0), target.mean(0)
    P, Q = mobile - mc, target - tc
    V, _, Wt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    return (P @ (V @ D @ Wt)) + tc


def condition_rmsf(cond_dir, chains, threshold):
    """Chain-local RMSF for every requested chain of one condition."""
    models = sorted(glob.glob(os.path.join(cond_dir, "seed-*_sample-*", "model.cif")))
    if not models:
        return [], {}
    parsed = [read_ca(m) for m in models]
    ref = parsed[0]                                   # the first model IS the reference frame

    rows, profiles = [], {}
    for ch, label in chains.items():
        present = [p for p in parsed if ch in p]
        if len(present) < 2:
            continue
        # residues resolved in every model of the condition
        common = set(present[0][ch][0])
        for p in present[1:]:
            common &= set(p[ch][0])
        common = np.array(sorted(common))

        def pick(p):
            seq, xyz, b = p[ch]
            keep = np.isin(seq, common)
            order = np.argsort(seq[keep])
            return xyz[keep][order], b[keep][order]

        ref_xyz, _ = pick(ref)
        stack, plddt = [], []
        for p in present:
            xyz, b = pick(p)
            stack.append(kabsch_fit(xyz, ref_xyz))
            plddt.append(b)
        stack = np.stack(stack)                       # models x residues x 3
        mean_pos = stack.mean(0)
        rmsf = np.sqrt(((stack - mean_pos) ** 2).sum(-1).mean(0))

        rows.append(dict(
            Chain=label, Chain_ID=ch, N_Models=len(present), Residues=len(common),
            Median_RMSF_A=round(float(np.median(rmsf)), 3),
            P90_RMSF_A=round(float(np.percentile(rmsf, 90)), 3),
            Pct_Residues_Above_Threshold=round(float(100 * (rmsf > threshold).mean()), 1),
            Mean_pLDDT=round(float(np.mean(plddt)), 1),
        ))
        profiles[label] = pd.DataFrame(dict(Chain=label, Residue=common,
                                            RMSF_A=np.round(rmsf, 4)))
    return rows, profiles


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True,
                    help="AlphaFold3 output directory holding the condition subdirectories")
    ap.add_argument("--conditions", nargs="+", required=True,
                    help="DIRNAME=Label pairs, in the order they should be reported")
    ap.add_argument("--chains", nargs="+", required=True, help="CHAINID=Label pairs")
    ap.add_argument("--threshold", type=float, default=0.5,
                    help="mobility threshold in Angstrom for the percentage column (default 0.5)")
    ap.add_argument("--out", default="chain_local_rmsf.csv")
    ap.add_argument("--profiles-out", default=None,
                    help="optional CSV of the full per-residue RMSF profiles")
    args = ap.parse_args()

    def pairs(items):
        out = {}
        for it in items:
            k, _, v = it.partition("=")
            out[k] = v or k
        return out

    conditions, chains = pairs(args.conditions), pairs(args.chains)

    summary, profiles = [], []
    for dirname, label in conditions.items():
        cond_dir = os.path.join(args.root, dirname)
        if not os.path.isdir(cond_dir):
            print(f"  ! missing, skipped: {dirname}")
            continue
        rows, prof = condition_rmsf(cond_dir, chains, args.threshold)
        for r in rows:
            summary.append(dict(Condition=label, **r))
            print(f"  {label:16s} {r['Chain']:6s} n={r['N_Models']:3d} "
                  f"res={r['Residues']:4d} median={r['Median_RMSF_A']:.3f} "
                  f"p90={r['P90_RMSF_A']:.3f} "
                  f">{args.threshold}A={r['Pct_Residues_Above_Threshold']:.1f}% "
                  f"pLDDT={r['Mean_pLDDT']:.1f}")
        for lab, df in prof.items():
            profiles.append(df.assign(Condition=label))

    if not summary:
        raise SystemExit("no conditions produced output")
    pd.DataFrame(summary).to_csv(args.out, index=False)
    print(f"wrote {args.out}")
    if args.profiles_out and profiles:
        pd.concat(profiles)[["Condition", "Chain", "Residue", "RMSF_A"]].to_csv(
            args.profiles_out, index=False)
        print(f"wrote {args.profiles_out}")


if __name__ == "__main__":
    main()
