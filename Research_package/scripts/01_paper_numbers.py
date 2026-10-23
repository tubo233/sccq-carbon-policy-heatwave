"""补充材料里的逐年序列、历史段累计值、人日和分省指标都从这里出。

读主情景的 global_T.npz、情景路径和指标表，指标定义与 models/main/04_metrics.py 一致。
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
EXP = HERE / 'models'
D = EXP / 'main'
sys.path.insert(0, str(EXP / 'common'))
from paths import IN, HW_OUT  # noqa: E402

OUT = HERE / 'data'; OUT.mkdir(parents=True, exist_ok=True)
REG = {'sichuan': 'SC', 'chongqing': 'CQ'}
GISTEMP_LAMBDA = {'sichuan': 2.078933, 'chongqing': 1.296840}
POP2035 = {'sichuan': 81.2097e6, 'chongqing': 30.2372e6}
MAIN = dict(NPI='NPI_aeei1.0', NPR='NPR', NPLS='NPLS_gcam_mid', ZERO='zero_SCCQ')


def lambdas():
    d = pd.read_csv(IN / 'regional_response_coefficients.csv')
    d = d[d.variable.eq('tasmax') & d.season.eq('annual')]
    out = {}
    for r in REG:
        k = d[d.region.eq(r)]
        out[r] = pd.DataFrame(dict(model=list(k.model) + ['Observational regression on GISTEMP'],
                                   window=list(k.window) + ['—'],
                                   lam=np.append(k.kappa.to_numpy(), GISTEMP_LAMBDA[r])))
    return out


def hw_curve(r):
    c = pd.read_csv(HW_OUT / 'response_curve.csv')
    c = c[c.background.eq('recent') & c.region.eq(r) & c.definition.eq('annual90') & c.kind.eq('day') & c.weighting.eq('population2020')].sort_values('offset_tx')
    x, y = c.offset_tx.to_numpy(), c.hwf_mean.to_numpy()
    hi = (y[-1] - y[-2]) / (x[-1] - x[-2]); lo = (y[1] - y[0]) / (x[1] - x[0])
    return lambda w: np.where(w > x[-1], y[-1] + hi * (w - x[-1]), np.where(w < x[0], y[0] + lo * (w - x[0]), np.interp(w, x, y)))


def q(a, p, axis=None):
    return np.nanquantile(a, p, axis=axis)


def main():
    z = np.load(D / 'fair/global_T.npz'); T = z['T']; yrs = z['years']; cases = list(z['cases'])
    base = T[(yrs >= 2015) & (yrs <= 2024), cases.index('NPR'), :].mean(axis=0)
    g = {c: T[:, cases.index(c), :] - base[None, :] for c in cases}          # (年份, 成员)
    L = lambdas()
    pd.concat([v.assign(region=r) for r, v in L.items()]).to_csv(OUT / 'lambda_values.csv', index=False)
    rows = []
    for r, tag in REG.items():
        lam = L[r].lam.to_numpy()[:, None, None]                               # (λ, 1, 1)
        views = {'combined': (MAIN['NPI'], MAIN['NPR'], MAIN['NPLS'], MAIN['ZERO']),
                 'own': (f"{MAIN['NPI']}__{tag}only", 'NPR', f"{MAIN['NPLS']}__{tag}only", f'zero__{tag}only')}
        hwf = hw_curve(r)
        for view, (npi, npr, npls, zero) in views.items():
            W = {k: lam * g[c][None] for k, c in dict(NPI=npi, NPR=npr, NPLS=npls, ZERO=zero).items()}   # (λ, 年份, 成员)
            A_npr, A_npls = W['NPI'] - W['NPR'], W['NPI'] - W['NPLS']
            own_warm = W['NPI'] - W['ZERO']
            with np.errstate(divide='ignore', invalid='ignore'):
                met = dict(W_NPI=W['NPI'], W_NPR=W['NPR'], W_NPLS=W['NPLS'], A_NPI_NPR=A_npr, A_NPI_NPLS=A_npls,
                           A_NPR_NPLS=A_npls - A_npr, P_NPR=100 * A_npr / W['NPI'], P_NPLS=100 * A_npls / W['NPI'],
                           C_NPR=100 * A_npr / own_warm, C_NPLS=100 * A_npls / own_warm,
                           G=100 * (A_npls - A_npr) / A_npr,
                           dHWF_NPI_NPLS=hwf(W['NPI']) - hwf(W['NPLS']), dHWF_NPI_NPR=hwf(W['NPI']) - hwf(W['NPR']),
                           dHWF_NPR_NPLS=hwf(W['NPR']) - hwf(W['NPLS']))
            for k, v in met.items():
                flat = np.moveaxis(v, 1, 0).reshape(len(yrs), -1)               # (年份, λ×成员)
                med, lo, hi = q(flat, .5, 1), q(flat, .05, 1), q(flat, .95, 1)
                for i, y in enumerate(yrs):
                    rows.append(dict(region=r, view=view, metric=k, year=int(y), median=med[i], p05=lo[i], p95=hi[i]))
    ts = pd.DataFrame(rows)
    ts.to_csv(OUT / 'timeseries_metrics_1990_2060.csv', index=False)

    # 分时段累计排放差
    sp = pd.read_csv(D / 'scenario_paths_energy_co2.csv')
    em = []
    for r in REG:
        e = sp[sp.region.eq(r)].pivot_table(index='year', columns='scenario', values='co2_mt')
        en = sp[sp.region.eq(r)].pivot_table(index='year', columns='scenario', values='energy_10k_tce')
        for a, b in [(2006, 2022), (2023, 2035), (2006, 2035)]:
            s = e.loc[a:b]
            em.append(dict(region=r, period=f'{a}-{b}', cum_NPI_Gt=s[MAIN['NPI']].sum() / 1e3, cum_NPR_Gt=s['NPR'].sum() / 1e3,
                           cum_NPLS_Gt=s[MAIN['NPLS']].sum() / 1e3,
                           NPI_minus_NPR_Gt=(s[MAIN['NPI']] - s['NPR']).sum() / 1e3,
                           NPR_minus_NPLS_Mt=(s['NPR'] - s[MAIN['NPLS']]).sum(),
                           energy_NPI_end=en.loc[b, MAIN['NPI']], energy_NPR_end=en.loc[b, 'NPR'], energy_NPLS_end=en.loc[b, MAIN['NPLS']]))
    em = pd.DataFrame(em); em.to_csv(OUT / 'cumulative_emissions_by_period.csv', index=False)

    # 人日：2035 年人口 × 主情景热浪天数差的中位数
    m = pd.read_csv(D / 'results/metrics_FaIR_lambda_sensitivity.csv', encoding='utf-8-sig')
    rn = {'sichuan': '四川', 'chongqing': '重庆'}
    pdays = []
    for r in REG:
        for k in ['dHWF_NPI_NPR', 'dHWF_NPI_NPLS', 'dHWF_NPR_NPLS']:
            row = m[(m.region == rn[r]) & (m.year == 2035) & (m.metric == k)].iloc[0]
            pdays.append(dict(region=r, metric=k, days_median=row['median'], days_p05=row.p05, days_p95=row.p95,
                              person_days_million_median=row['median'] * POP2035[r] / 1e6,
                              person_days_million_p05=row.p05 * POP2035[r] / 1e6, person_days_million_p95=row.p95 * POP2035[r] / 1e6))
    pd.DataFrame(pdays).to_csv(OUT / 'person_days_2035.csv', index=False)

    pd.set_option('display.width', 220)
    pick = ts[(ts.view == 'combined') & ts.year.isin([2022, 2030, 2035, 2050, 2060]) & ts.metric.isin(['W_NPI', 'A_NPI_NPR', 'A_NPI_NPLS', 'P_NPLS', 'C_NPR', 'C_NPLS', 'dHWF_NPI_NPLS'])]
    print(pick.pivot_table(index=['region', 'metric'], columns='year', values='median').round(5).to_string())
    print(ts[(ts.view == 'own') & (ts.year == 2035) & ts.metric.isin(['A_NPI_NPLS', 'A_NPR_NPLS', 'P_NPLS', 'C_NPLS', 'G'])].round(5).to_string())
    print(em.round(3).to_string())
    print(pd.DataFrame(pdays).round(3).to_string())


if __name__ == '__main__':
    main()
