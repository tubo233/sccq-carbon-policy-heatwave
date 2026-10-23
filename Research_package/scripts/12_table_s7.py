"""从区域汇总重建甲烷敏感性表 S7。"""
from pathlib import Path
import pandas as pd
P=Path(__file__).resolve().parents[1]
d=pd.read_csv(P/'data/fair_supplement/crossover_emergence_summary.csv').set_index(['case','region'])
old=pd.read_csv(P/'data/tableS7_methane_sensitivity.csv')
scales={'SC–CQ share':'own_share','China':'national','World':'universal'}
rows=[]
for _,r in old.iterrows():
    v=r['variant']
    kind='all' if v.startswith('All DPEC') else ('CH4only_V1' if v.startswith('Methane only') else 'allCH4_'+v[:2])
    case=scales[r['scale']]+('__' if v.startswith('Methane only') else '__OTNZ__')+kind
    sc=d.loc[(case,'sichuan')];cq=d.loc[(case,'chongqing')]
    rows.append(dict(scale=r['scale'],variant=v,SC_dTmax_2035_median=round(sc.d2035,5),SC_p05=round(sc.d2035_p05,5),SC_p95=round(sc.d2035_p95,5),SC_share_warming_2035=round(sc.share_warming_2035,3),crossover_SC=sc.crossover,crossover_CQ=cq.crossover,toe20_SC=sc.toe20_median,toe20_CQ=cq.toe20_median))
pd.DataFrame(rows).to_csv(P/'data/tableS7_methane_sensitivity.csv',index=False)
print('Rebuilt Table S7:',len(rows),'rows')
