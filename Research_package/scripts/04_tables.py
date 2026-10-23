"""表 2、表 S9 和表 S5B 的格式化数值（打印并保存为 JSON）。

读 models/main 的输出和 data/。
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
D = HERE / 'models/main'
DATA = HERE / 'data'
CN = {'sichuan': '四川', 'chongqing': '重庆'}
m = pd.read_csv(D / 'results/metrics_FaIR_lambda_sensitivity.csv', encoding='utf-8-sig')
em = pd.read_csv(D / 'results/emission_metrics.csv', encoding='utf-8-sig')
pdays = pd.read_csv(DATA / 'person_days_2035.csv')
st = pd.read_csv(DATA / 'sensitivity_table_2035.csv')
hist = pd.read_csv(D / 'history_gdp_energy_co2.csv')


def v(r, y, k, f='{:.2f}', ci=True):
    row = m[(m.region == CN[r]) & (m.year == y) & (m.metric == k)].iloc[0]
    s = f.format(row['median'])
    return s + (f" [{f.format(row.p05)}, {f.format(row.p95)}]" if ci else '')


T2 = {}
for r in ['sichuan', 'chongqing']:
    for y in [2030, 2035]:
        e = em[(em.region == CN[r]) & (em.year == y)].iloc[0]
        c = {}
        c['Energy use NPI / NPR / NPLS (Mtce)'] = ' / '.join(f"{e[f'energy_{s}_10k_tce'] / 100:.0f}" for s in ['NPI', 'NPR', 'NPLS'])
        c['CO2 NPI / NPR / NPLS (Mt)'] = ' / '.join(f"{e[f'co2_{s}_mt']:.0f}" for s in ['NPI', 'NPR', 'NPLS'])
        c['Cumulative CO2 since 2006, reduction vs NPI: NPR / NPLS (%)'] = f"{e.avoided_cum_NPR_pct_of_NPI:.1f} / {e.avoided_cum_NPLS_pct_of_NPI:.1f}"
        c['W_NPI (°C)'] = v(r, y, 'W_NPI')
        c['A, NPI − NPR (°C)'] = v(r, y, 'A_NPI_NPR', '{:.4f}')
        c['A, NPI − NPLS (°C)'] = v(r, y, 'A_NPI_NPLS', '{:.4f}')
        c['A, NPR − NPLS (°C)'] = v(r, y, 'A_NPR_NPLS', '{:.5f}')
        c['P, NPR / NPLS (%)'] = f"{v(r, y, 'P_NPR', '{:.2f}', False)} / {v(r, y, 'P_NPLS', '{:.2f}')}"
        c['C, NPR / NPLS (%)'] = f"{v(r, y, 'C_NPR', '{:.1f}', False)} / {v(r, y, 'C_NPLS', '{:.1f}')}"
        c['G (%)'] = v(r, y, 'G_local_gain', '{:.1f}')
        c['ΔHWF, NPI − NPLS (days/yr)'] = v(r, y, 'dHWF_NPI_NPLS', '{:.2f}')
        c['Own-province view: A NPI − NPLS (°C) / C (%) / G (%)'] = f"{v(r, y, 'own_A_NPI_NPLS', '{:.4f}', False)} / {v(r, y, 'own_C_NPLS', '{:.1f}', False)} / {v(r, y, 'own_G_local_gain', '{:.1f}', False)}"
        T2[f'{r} {y}'] = c
T2 = pd.DataFrame(T2)
print(T2.to_string())
p = pdays.set_index(['region', 'metric'])
print(p.round(3).to_string())

T3 = {}
for r in ['sichuan', 'chongqing']:
    s = st[st.region == r]
    for metric, f in [('A_NPI_NPLS', '{:.4f}'), ('A_NPR_NPLS', '{:.5f}'), ('P_NPLS', '{:.2f}'), ('C_NPLS', '{:.1f}')]:
        x = s[s.metric == metric].set_index('source')
        col = {}
        for src in x.index:
            lo, hi = x.loc[src, 'lo'], x.loc[src, 'hi']
            col[src] = '—' if np.isnan(lo) else (f.format(lo) if np.isclose(lo, hi) else f'{f.format(lo)}–{f.format(hi)}')
        T3[f'{r} {metric}'] = col
T3 = pd.DataFrame(T3)
print(T3.to_string())
h = hist.set_index(['region', 'year'])
for r in ['sichuan', 'chongqing']:
    print(r, 'energy 2005/2022', h.loc[(r, 2005), 'energy_10k_tce_equal'] / 100, h.loc[(r, 2022), 'energy_10k_tce_equal'] / 100,
          'co2 2005/2022', h.loc[(r, 2005), 'co2_mt_ceads'], h.loc[(r, 2022), 'co2_mt_ceads'])
obs = pd.read_csv(HERE / 'models/heatwave/observed_annual_tmax_hwf_1990_2024.csv')
o = obs[obs.weighting.eq('population2020')].pivot_table(index='year', columns='region', values='hwf')
print('top-3 HWF years:', {r: o[r].nlargest(3).round(1).to_dict() for r in o.columns})
json.dump(dict(T2=T2.to_dict(), T3=T3.to_dict()), open(DATA / 'table2_tableS5B_S9_values.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
