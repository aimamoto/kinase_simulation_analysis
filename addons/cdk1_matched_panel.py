"""Is the CDK1 'thaw' an artifact of comparing different metric panels?

Global MAC's panel varies by condition (apo mono = 91 edges, +CCNB1 = 120, holo = 210), so the
headline 0.365 -> 0.217 compares a 14-metric network with a 16-metric one. Recompute on a matched
panel = the metrics viable in BOTH conditions of each contrast.

R semantics replicated exactly (multimer_core_engine.R ~L497-508):
  keep cols with >=3 non-NA; spearman with pairwise.complete.obs; NA correlations -> 0;
  MAC = mean(abs(upper triangle)). Note R RETAINS zero-variance columns (their correlations
  become NA -> 0), so they drag MAC down; do not drop them.
"""
import os
import numpy as np, pandas as pd, itertools
B=os.environ.get("ALLOQUANT_CDK1", "CDK1-CCNB1_output")  # AlloQuant CDK1-CCNB1 output dir
m=pd.read_csv(f'{B}/master_kinase_analysis_results_v7r3.csv')
m['sid']=m.Simulation_ID.str.replace('_model$','',regex=True)
st=pd.read_csv(f'{B}/plots_and_stats_CDK1_GMM/Phase6_State_Assignments.csv')
st['sid']=st.Simulation_ID.str.replace('_model$','',regex=True)
st['rev']=st.Condition_reviewed.str.replace('\n','/',regex=False)
d=m[m.Type=='CDK1'].merge(st[['sid','rev','Macro_State']],on='sid',how='inner')
DIST=[c for c in d.columns if c.endswith('_Dist')]
for c in DIST: d[c]=pd.to_numeric(d[c],errors='coerce')

def panel(df):                      # R: select(where(~sum(!is.na(.)) >= 3))
    return [c for c in DIST if df[c].notna().sum()>=3]
def mac(df, cols):
    X=df[cols]
    C=X.corr(method='spearman').values.astype(float)
    C=np.nan_to_num(C)              # R: c_mat[is.na(c_mat)] <- 0
    iu=np.triu_indices(len(cols),1)
    return float(np.abs(C[iu]).mean()), len(iu[0])

gl=pd.read_csv(f'{B}/plots_and_stats_CDK1_GMM/Phase5_Global_Network_Density.csv')
gl['rev']=gl.Condition.str.split('\n').str[:-1].str.join('/')
ref=dict(zip(gl.rev,gl.Global_Coupling_Score))
print('STEP 1 - reproduce the published per-condition global MAC (validates the replication)')
print(f'  {"condition":34s}{"n":>5s}{"metrics":>9s}{"edges":>7s}{"mine":>9s}{"published":>11s}')
P={}
for c,g in d.groupby('rev'):
    p=panel(g); v,e=mac(g,p); P[c]=(p,g)
    pub=ref.get(c,np.nan)
    flag='' if (np.isnan(pub) or abs(v-pub)<5e-4) else '   <-- MISMATCH'
    print(f'  {c:34s}{len(g):5d}{len(p):9d}{e:7d}{v:9.4f}{pub:11.4f}{flag}')

print('\nSTEP 2 - the thaw and the other Fig. 2A transitions, raw vs matched panel')
PAIRS=[('cdk1','cdk1/ccnb1','apo -> +CCNB1  (THE THAW)'),
       ('cdk1/ccnb1','cdk1/ccnb1','')]
conds=list(P)
def show(a,b,label):
    pa,ga=P[a]; pb,gb=P[b]
    matched=[c for c in DIST if c in pa and c in pb]
    va,_=mac(ga,pa); vb,_=mac(gb,pb)
    ma,_=mac(ga,matched); mb,_=mac(gb,matched)
    print(f'  {label}')
    print(f'    raw     : {va:.4f} ({len(pa)} metrics) -> {vb:.4f} ({len(pb)} metrics)   D = {vb-va:+.4f}')
    print(f'    matched : {ma:.4f} -> {mb:.4f}   on {len(matched)} shared metrics   D = {mb-ma:+.4f}')
    print(f'    dropped from the larger panel: {sorted(set(pa)^set(pb)) if set(pa)^set(pb) else "none"}')
    return vb-va, mb-ma
# locate the two condition keys for the thaw
apo=[c for c in conds if c.strip()=='cdk1']
ccnb=[c for c in conds if c.strip()=='cdk1/ccnb1']
if apo and ccnb: show(apo[0],ccnb[0],'apo monomer -> +CCNB1   (THE THAW, 0.365 -> 0.217)')
print()
for a,b in itertools.combinations(conds,2):
    pa,pb=P[a][0],P[b][0]
    if set(pa)!=set(pb):
        matched=[c for c in DIST if c in pa and c in pb]
        va,_=mac(P[a][1],pa); vb,_=mac(P[b][1],pb)
        ma,_=mac(P[a][1],matched); mb,_=mac(P[b][1],matched)
        print(f'  {a} vs {b}: raw D={vb-va:+.4f} -> matched D={mb-ma:+.4f} ({len(matched)} shared)')

# ---------------------------------------------------------------------------
print('\n'+'='*74)
print('SAME CHECK FOR THE CSK SERIES (Fig. 4A)')
print('='*74)
BS=os.environ.get("ALLOQUANT_SRC", "CSK-SRC_output")  # AlloQuant CSK-SRC output dir
ms=pd.read_csv(f'{BS}/master_kinase_analysis_results_v7r3.csv')
ms['sid']=ms.Simulation_ID.str.replace('_model$','',regex=True)
sts=pd.read_csv(f'{BS}/plots_and_stats_CSK_GMM/Phase6_State_Assignments.csv')
sts['sid']=sts.Simulation_ID.str.replace('_model$','',regex=True)
sts['rev']=sts.Condition_reviewed.str.replace('\n','/',regex=False)
dc=ms[ms.Type=='CSK-WT'].merge(sts[['sid','rev']],on='sid',how='inner')
D2=[c for c in dc.columns if c.endswith('_Dist')]
for c in D2: dc[c]=pd.to_numeric(dc[c],errors='coerce')
def panel2(g): return [c for c in D2 if g[c].notna().sum()>=3]
def mac2(g,cols):
    C=np.nan_to_num(g[cols].corr(method='spearman').values.astype(float))
    iu=np.triu_indices(len(cols),1)
    return float(np.abs(C[iu]).mean())
FIG4=[('apo/apo','csk-wtcat-apo/src-wtcat-apo'),('CSK-ATP','csk-wtcat-holo/src-wtcat-apo'),
      ('both-ATP','csk-wtcat-holo/src-wtcat-holo'),('primed','csk-wtcat-holo/src-wtcat-py159-holo')]
Q={}
for nm,c in FIG4:
    g=dc[dc.rev==c]; Q[nm]=(panel2(g),g)
    print(f'  {nm:10s} n={len(g):3d}  {len(Q[nm][0]):2d} metrics  global MAC {mac2(g,Q[nm][0]):.4f}')
print('\n  pairwise, raw vs matched panel:')
for a,b in itertools.combinations([n for n,_ in FIG4],2):
    pa,ga=Q[a]; pb,gb=Q[b]
    matched=[c for c in D2 if c in pa and c in pb]
    va,vb=mac2(ga,pa),mac2(gb,pb)
    ma,mb=mac2(ga,matched),mac2(gb,matched)
    tag='  (already matched)' if set(pa)==set(pb) else ''
    print(f'    {a:9s} -> {b:9s}  raw D={vb-va:+.4f}   matched D={mb-ma:+.4f}  on {len(matched)} shared{tag}')
