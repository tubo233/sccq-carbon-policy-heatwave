"""相对于发布包的路径。"""
from pathlib import Path
import os, sys, json
from datetime import datetime, timezone
EXP = Path(__file__).resolve().parents[1]
ROOT = EXP.parent
IN = ROOT / 'inputs/fair'
PROV = IN
EMIS = ROOT / 'inputs/emissions'
CEADS = EMIS / 'ceads_history.csv'
GADM = ROOT / 'inputs/boundaries/gadm41_CHN_2.json'
FAIR_OUT = EXP / 'scales'
HW_OUT = EXP / 'heatwave'
LEGACY = ROOT / 'weather'
SSP245DIR = ROOT / 'inputs/climate'
sys.dont_write_bytecode = True
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(key,'2')
os.environ.setdefault('MPLBACKEND','Agg')
def save(p,v):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,ensure_ascii=False,indent=2,default=str),encoding='utf-8')
def now():
    return datetime.now(timezone.utc).isoformat()
def kappa_table():
    import pandas as pd
    d=pd.read_csv(IN/'regional_response_coefficients.csv')
    d=d[d.season.eq('annual') & d.window.eq('years_61_80_primary') & d.variable.isin(['tasmax','tasmin'])]
    return d.pivot_table(index=['region','variable'],columns='model',values='kappa')
