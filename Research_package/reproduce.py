"""从随包输入或预计算集合复现论文输出。"""
import argparse, os, subprocess, sys
from pathlib import Path
P=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--stage',choices=['figures','metrics','scenarios','fair','weather','air-quality','scale-fair'],default='figures')
a=parser.parse_args()
jobs={
 'figures':['scripts/03_paper_figures.py','scripts/04_tables.py','scripts/10_si_data_figures.py','scripts/11_plot_supplement_scale_series.py'],
 'metrics':['models/scales/analyze_scales.py','models/main/04_metrics.py','scripts/01_paper_numbers.py','scripts/09_analyze_supplement.py','scripts/12_table_s7.py'],
 'scenarios':['models/main/02_scenarios.py'],
 'fair':['models/main/03_run_fair.py','scripts/08_fair_methane_benchmarks.py'],
 'weather':['models/heatwave/observed_noise.py','models/heatwave/run_response_curve.py'],
 'scale-fair':['models/scales/run_fair_scales.py'],
 'air-quality':['models/air_quality/03_pollution_channel.py']}
env=dict(os.environ,PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1',MPLBACKEND='Agg')
for file in jobs[a.stage]:
 print('Running',file,flush=True)
 subprocess.run([sys.executable,str(P/file)],cwd=P,env=env,check=True)
print('Completed',a.stage)
