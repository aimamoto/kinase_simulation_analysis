"""Does the CDK1 'thaw' survive the same mixture decomposition that broke CSK's 'stiffening'?

If pooled MAC ~= occupancy-weighted within-state MAC for the CDK1 conditions, the thaw is a
genuine within-state coupling change and Model 1 is unaffected. Panel = `_Dist` only, matching
multimer_core_engine.R L308-309.
"""
import os
import numpy as np, pandas as pd
B=os.environ.get("ALLOQUANT_CDK1", "CDK1-CCNB1_output")  # AlloQuant CDK1-CCNB1 output dir
m=pd.read_csv(f'{B}/master_kinase_analysis_results_v7r3.csv')
m['sid']=m.Simulation_ID.str.replace('_model$','',regex=True)
st=pd.read_csv(f'{B}/plots_and_stats_CDK1_GMM/Phase6_State_Assignments.csv')
st['sid']=st.Simulation_ID.str.replace('_model$','',regex=True)
gl=pd.read_csv(f'{B}/plots_and_stats_CDK1_GMM/Phase5_Global_Network_Density.csv')
intr=pd.read_csv(f'{B}/plots_and_stats_CDK1_GMM/Phase6b_State_Intrinsic_Rigidity_MAC.csv')
IM=dict(zip(intr.Macro_State,intr.Intrinsic_MAC))
cdk=m[m.Type.str.contains('CDK1',case=False,na=False)].merge(
    st[['sid','Macro_State','Condition_reviewed']],on='sid',how='left')
DIST=[c for c in cdk.columns if c.endswith('_Dist')]
def mac(df):
    X=df[DIST].apply(pd.to_numeric,errors='coerce')
    X=X.loc[:,(X.notna().sum()>=3)&(X.nunique(dropna=True)>1)]
    if X.shape[1]<2 or len(X)<3: return np.nan
    C=np.nan_to_num(X.corr(method='spearman').values)
    return float(np.abs(C[np.triu_indices_from(C,1)]).mean())

print('pipeline Global_Coupling_Score for reference:')
g=gl.copy(); g['Condition']=g.Condition.str.replace('\n',' | ',regex=False)
print(g.to_string(index=False)); print()
print(f'{"condition":52s}{"pooled":>8s}{"within(occ-wt)":>16s}{"intrinsic(occ-wt)":>19s}{"n":>5s}')
for cond,sub in cdk.groupby('Condition_reviewed'):
    if len(sub)<10: continue
    occ=sub.Macro_State.value_counts()
    ws=[(s,mac(sub[sub.Macro_State==s]),c) for s,c in occ.items() if c>=8]
    if not ws: continue
    wm=np.average([v for _,v,_ in ws],weights=[c for _,_,c in ws])
    oi=sum(c/len(sub)*IM[s] for s,c in occ.items() if s in IM)
    lbl=cond.replace('\n',' / ')[:50]
    print(f'  {lbl:50s}{mac(sub):8.4f}{wm:16.4f}{oi:19.4f}{len(sub):5d}')
    print(f'      states: '+', '.join(f'{s} n={c} MAC={v:.3f}' for s,v,c in ws))
