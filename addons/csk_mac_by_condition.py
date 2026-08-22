"""Per-condition MAC decomposition for all four biological CSK conditions (Fig. 4A rebuild, register F7).

Extends csk_stiffen_v2.py from the two-condition contrast to all four conditions plotted in Fig. 4A,
and emits a CSV the figure script reads.

WHY NOT the register's original F7 plan. F7 (written 2026-08-03) proposed overlaying the pipeline's
per-state INTRINSIC MAC on the global MAC bars. That is the comparison the same memo later forbids
(§5): intrinsic MAC is anchored on `universal_dist_cols` = 16 metrics / 120 edges, global MAC on
`any_dist_cols` filtered per condition = 16-21 metrics / 120-210 edges. Plotting one against the other
puts an invalid comparison into the figure.

WHAT THIS DOES INSTEAD. For each condition, fix that condition's own global-MAC panel, then recompute
the same statistic separately within each occupied macro-state ON THAT FIXED PANEL. Those values are
panel-matched to the condition's own bar, so the within-condition comparison in the panel is valid.
Cross-condition comparison of the bars remains panel-dependent (register C8/S10) and is a legend matter.

Reproduces the pipeline's Phase5_Global_Network_Density.csv values as a check.
"""
import os
import numpy as np, pandas as pd

BASE = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
OUT = os.environ.get("ALLOQUANT_CSK_MAC_CSV", "csk_mac_by_condition.csv")

m = pd.read_csv(f'{BASE}/master_kinase_analysis_results_v7r3.csv')
m['sid'] = m.Simulation_ID.str.replace('_model$', '', regex=True)
csk = m[m.Type == 'CSK-WT'].copy()
cs = pd.read_csv(f'{BASE}/plots_and_stats_CSK_GMM/Phase6_State_Assignments.csv')
cs['sid'] = cs.Simulation_ID.str.replace('_model$', '', regex=True)
csk = csk.merge(cs[['sid', 'Macro_State', 'Condition_reviewed']], on='sid', how='left')
csk['cond'] = (csk.Condition_reviewed.str.replace(r'\s+', ' ', regex=True).str.strip()
               .str.replace('wtcat-', '', regex=False).str.replace(' ', '/', regex=False))

DIST = [c for c in csk.columns if c.endswith('_Dist')]

def panel(df):
    """The columns the pipeline's global MAC would use for this condition."""
    X = df[DIST].apply(pd.to_numeric, errors='coerce')
    return [c for c in X.columns if X[c].notna().sum() >= 3 and X[c].nunique(dropna=True) > 1]

def mac_on(df, cols):
    X = df[cols].apply(pd.to_numeric, errors='coerce')
    X = X.loc[:, X.nunique(dropna=True) > 1]
    if X.shape[1] < 2 or len(X) < 3:
        return np.nan, X.shape[1]
    C = np.nan_to_num(X.corr(method='spearman').values)
    iu = np.triu_indices_from(C, 1)
    return float(np.abs(C[iu]).mean()), X.shape[1]

CONDS = [('csk-apo/src-apo', 'apo/apo'), ('csk-holo/src-apo', 'CSK-ATP'),
         ('csk-holo/src-holo', 'both-ATP'), ('csk-holo/src-py159-holo', 'primed')]
MIN_N = 8   # smallest cell for which the statistic is reported at all

rows = []
for key, label in CONDS:
    sub = csk[csk['cond'] == key]
    cols = panel(sub)
    pooled, k = mac_on(sub, cols)
    occ = sub.Macro_State.value_counts()
    print(f'\n{label:9s} ({key})  n={len(sub)}  panel={k} metrics, {k*(k-1)//2} edges')
    print(f'  global MAC (pooled)            {pooled:.4f}')
    per = []
    for st, n in occ.items():
        if n < MIN_N:
            continue
        v, kk = mac_on(sub[sub.Macro_State == st], cols)
        per.append((st, v, n))
        print(f'    {st:9s} n={n:3d}  within-state MAC {v:.4f}  ({100*n/len(sub):.0f}% of condition)')
        rows.append(dict(cond=key, label=label, kind='state', state=st, mac=v, n=n,
                         pct=100*n/len(sub), panel=k))
    if per:
        w = np.average([v for _, v, _ in per], weights=[n for _, _, n in per])
        print(f'  occupancy-weighted within-state {w:.4f}   (difference {pooled-w:+.4f})')
    else:
        w = np.nan
    dropped = sorted(set(occ.index) - {s for s, _, _ in per})
    if dropped:
        tot = int(occ[dropped].sum())
        print(f'  NOT plotted (n < {MIN_N}): {", ".join(f"{s} n={occ[s]}" for s in dropped)}'
              f'  = {tot} models, {100*tot/len(sub):.0f}% of the condition')
    rows.append(dict(cond=key, label=label, kind='global', state='', mac=pooled, n=len(sub),
                     pct=100.0, panel=k))
    rows.append(dict(cond=key, label=label, kind='within_wt', state='', mac=w, n=len(sub),
                     pct=100.0, panel=k))

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'\nwrote {OUT}')
