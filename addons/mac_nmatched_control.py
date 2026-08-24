"""n-matched permutation control for the pooled-vs-per-state MAC comparison in Fig. 4A.

WHY THIS EXISTS. MAC (mean |Spearman rho|) is strongly upward-biased at small n, so an open
circle in Fig. 4A cannot be read against its bar by height. Supplementary Note 10 already states
the analytic floor for INDEPENDENT variables, sqrt(2/[pi(n-1)]). That floor answers "does this
cell carry information", but it is the wrong reference for "is this state less coupled than the
ensemble it belongs to": a random subset of a MIXED condition carries the mixture's own induced
correlation, which is far above the independent-variable floor. In the primed condition the
analytic floor at n=42 is 0.125 while a random 42-model slice of that same condition reads 0.257.

WHAT THIS DOES. For each condition, fix that condition's own global-MAC panel (as in
csk_mac_by_condition.py), then test each occupied macro-state against NREP random subsets of
the same condition of the same size, state labels ignored. Reports the per-state z and one-sided
p, and the same test for the occupancy-weighted statistic using disjoint random subsets of the
observed state sizes.

SCOPE. This compares a state with the mixture it sits in. It does NOT compare one state between
two conditions, which remains impossible (register S13): every CSK state well populated in one
condition is nearly empty in the other. The Note 10 failure-to-detect caveat is unaffected.

Emits csk_mac_nmatched.csv.
"""
import os
import numpy as np, pandas as pd

BASE = os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
OUT = os.environ.get("ALLOQUANT_CSK_NMATCHED_CSV", "csk_mac_nmatched.csv")
NREP = 2000
SEED = 20260824
MIN_N = 8   # same threshold as csk_mac_by_condition.py / Fig. 4A

m = pd.read_csv(f'{BASE}/master_kinase_analysis_results_v7r3.csv')
m['sid'] = m.Simulation_ID.str.replace('_model$', '', regex=True)
csk = m[m.Type == 'CSK-WT'].copy()
cs = pd.read_csv(f'{BASE}/plots_and_stats_CSK_GMM/Phase6_State_Assignments.csv')
cs['sid'] = cs.Simulation_ID.str.replace('_model$', '', regex=True)
csk = csk.merge(cs[['sid', 'Macro_State', 'Condition_reviewed']], on='sid', how='left')
csk['cond'] = (csk.Condition_reviewed.str.replace(r'\s+', ' ', regex=True).str.strip()
               .str.replace('wtcat-', '', regex=False).str.replace(' ', '/', regex=False))
DIST = [c for c in csk.columns if c.endswith('_Dist')]
rng = np.random.default_rng(SEED)


def mac(X):
    X = X.loc[:, X.nunique(dropna=True) > 1]
    if X.shape[1] < 2 or len(X) < 3:
        return np.nan
    C = np.nan_to_num(X.corr(method='spearman').values)
    iu = np.triu_indices_from(C, 1)
    return float(np.abs(C[iu]).mean())


def floor(n):
    """Analytic mean |rho| between independent variables (Supplementary Note 10)."""
    return float(np.sqrt(2.0 / (np.pi * (n - 1))))


CONDS = [('csk-apo/src-apo', 'apo/apo'), ('csk-holo/src-apo', 'CSK-ATP'),
         ('csk-holo/src-holo', 'both-ATP'), ('csk-holo/src-py159-holo', 'primed')]

rows = []
for key, label in CONDS:
    sub = csk[csk['cond'] == key]
    Xf = sub[DIST].apply(pd.to_numeric, errors='coerce')
    cols = [c for c in Xf.columns if Xf[c].notna().sum() >= 3 and Xf[c].nunique(dropna=True) > 1]
    X = sub[cols].apply(pd.to_numeric, errors='coerce').reset_index(drop=True)
    st = sub.Macro_State.reset_index(drop=True)
    pooled = mac(X)
    occ = [(s, n) for s, n in st.value_counts().items() if n >= MIN_N]
    print('=' * 104)
    print(f'{label:9s} n={len(X)}  panel={len(cols)} metrics / {len(cols)*(len(cols)-1)//2} edges'
          f'  pooled MAC {pooled:.4f}')
    rows.append(dict(cond=key, label=label, kind='pooled', state='', n=len(X), mac=pooled,
                     null_mu=np.nan, null_lo=np.nan, null_hi=np.nan, z=np.nan, p_lower=np.nan,
                     indep_floor=floor(len(X)), panel=len(cols)))

    for s, n in occ:
        obs = mac(X[st == s])
        null = np.array([mac(X.iloc[rng.choice(len(X), n, replace=False)]) for _ in range(NREP)])
        z = (obs - null.mean()) / null.std()
        p = (np.sum(null <= obs) + 1) / (NREP + 1)
        lo, hi = np.percentile(null, [2.5, 97.5])
        print(f'   {s:9s} n={n:3d}  MAC {obs:.4f}   size-matched null {null.mean():.4f} '
              f'({lo:.4f}-{hi:.4f})   indep floor {floor(n):.4f}   z={z:+.2f}  p={p:.4f}')
        rows.append(dict(cond=key, label=label, kind='state', state=s, n=n, mac=obs,
                         null_mu=null.mean(), null_lo=lo, null_hi=hi, z=z, p_lower=p,
                         indep_floor=floor(n), panel=len(cols)))

    ns = [n for _, n in occ]
    obs_w = np.average([mac(X[st == s]) for s, _ in occ], weights=ns)
    nullw = []
    for _ in range(NREP):
        perm = rng.permutation(len(X))
        vals, i = [], 0
        for n in ns:
            vals.append(mac(X.iloc[perm[i:i + n]])); i += n
        nullw.append(np.average(vals, weights=ns))
    nullw = np.array(nullw)
    zw = (obs_w - nullw.mean()) / nullw.std()
    pw = (np.sum(nullw <= obs_w) + 1) / (NREP + 1)
    lo, hi = np.percentile(nullw, [2.5, 97.5])
    print(f'   weighted  n={sum(ns):3d}  MAC {obs_w:.4f}   size-matched null {nullw.mean():.4f} '
          f'({lo:.4f}-{hi:.4f})   z={zw:+.2f}  p={pw:.4f}')
    rows.append(dict(cond=key, label=label, kind='weighted', state='', n=sum(ns), mac=obs_w,
                     null_mu=nullw.mean(), null_lo=lo, null_hi=hi, z=zw, p_lower=pw,
                     indep_floor=np.nan, panel=len(cols)))

pd.DataFrame(rows).to_csv(OUT, index=False)
print(f'\nwrote {OUT}   (NREP={NREP}, seed={SEED})')
