"""尺度分析：841 FaIR 成员 × 3 个 CMIP6 模式 κ 的联合分布 → 区域人口加权 Tmax 变化 → 热浪天数（同模式背景响应曲线）。

符号约定：ΔT = T(政策/情景) − T(SSP2-4.5 参照)；负值 = 更凉。历史情景：ΔT = T(实际) − T(反事实)。
噪声尺度（观测 1990-2024 去趋势，人口加权）：Tmax σ 与 HWF σ 见 models/heatwave/observed_noise_scales.csv。
显现阈值：单年 = 1σ；20 年均值 = 2σ/√20（年际独立近似；观测 lag-1 自相关接近 0）。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *
import numpy as np
import pandas as pd

RUN = FAIR_OUT / 'run'
OUTA = EXP / 'exposure' / 'results'
OUTB = EXP / 'scales' / 'results'
MODELS = ['GFDL-ESM4', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
REGIONS = ['sichuan', 'chongqing']
RNAME = {'sichuan': '四川', 'chongqing': '重庆'}
POP2035 = {'sichuan': 8120.9693e4, 'chongqing': 3023.7169e4}   # 2035 年人口（人）


def load_fair():
    z = np.load(RUN / 'global_temperature.npz')
    return z['T'], z['years'], list(z['cases']), z['aerosol_erf']


def kappas():
    k = kappa_table()
    return {r: np.array([k.loc[(r, 'tasmax'), m] for m in MODELS]) for r in REGIONS}


def curves():
    c = pd.read_csv(HW_OUT / 'response_curve.csv')
    c = c[c.definition.eq('annual90') & c.weighting.eq('population2020')]
    out = {}
    for (bg, r, kind), g in c.groupby(['background', 'region', 'kind']):
        g = g.sort_values('offset_tx')
        out[(bg, r, kind)] = (g.offset_tx.to_numpy(), g.hwf_mean.to_numpy(), g.hwf_sd_members.to_numpy())
    return out


def dhwf(curve, dt):
    x, y, _ = curve
    y0 = np.interp(0.0, x, y)
    lo = (y[1] - y[0]) / (x[1] - x[0]); hi = (y[-1] - y[-2]) / (x[-1] - x[-2])
    v = np.interp(dt, x, y)
    v = np.where(dt < x[0], y[0] + lo * (dt - x[0]), v)
    v = np.where(dt > x[-1], y[-1] + hi * (dt - x[-1]), v)
    return v - y0


def noise():
    n = pd.read_csv(HW_OUT / 'observed_noise_scales.csv')
    n = n[n.weighting.eq('population2020')].set_index(['region', 'variable']).detrended_sd
    return {r: dict(t=float(n.loc[(r, 'tmax_mean')]), h=float(n.loc[(r, 'hwf')])) for r in REGIONS}


def q(x, p):
    return float(np.quantile(x, p))


def regional(dg, r, K):
    """dg: (成员,) 全球 ΔT → (成员×3,) 区域 ΔT；同时返回对应模式标签索引。"""
    return (dg[:, None] * K[r][None, :]).ravel(), np.tile(np.arange(3), len(dg))


LEVERS_2035 = [
    ('川渝自身·仅CO2', 'ssp245__own_share__OTNZ__co2only', 2035),
    ('川渝自身·全物种', 'ssp245__own_share__OTNZ__all', 2035),
    ('川渝自身·全物种+CH4', 'ssp245__own_share__OTNZ__allCH4', 2035),
    ('川渝·仅清洁空气', 'ssp245__own_share__CA__all', 2035),
    ('全国·仅CO2', 'ssp245__national__OTNZ__co2only', 2035),
    ('全国·全物种', 'ssp245__national__OTNZ__all', 2035),
    ('全国·全物种+CH4', 'ssp245__national__OTNZ__allCH4', 2035),
    ('全国·仅清洁空气', 'ssp245__national__CA__all', 2035),
    ('全世界同等力度·仅CO2', 'ssp245__universal__OTNZ__co2only', 2035),
    ('全世界同等力度·全物种', 'ssp245__universal__OTNZ__all', 2035),
    ('全世界同等力度·全物种+CH4', 'ssp245__universal__OTNZ__allCH4', 2035),
    ('全世界·仅清洁空气', 'ssp245__universal__CA__all', 2035),
    ('全球SSP1-1.9 − SSP2-4.5', 'ssp119__ref', 2035),
    ('全球SSP1-2.6 − SSP2-4.5', 'ssp126__ref', 2035),
    ('全球SSP3-7.0 − SSP2-4.5', 'ssp370__ref', 2035),
    ('历史上界·川渝2005-12指数增长延续(2022)', 'hist__SCCQ_exp2005_2012', 2022),
    ('历史上界·川渝2005-12线性增长延续(2022)', 'hist__SCCQ_lin2005_2012', 2022),
    ('历史上界·全国多排20%(2022)', 'hist__China_plus20pct', 2022),
]


def delta_global(T, yrs, cases, case, year):
    ref = cases.index('ssp245__ref')
    d = T[yrs == year, cases.index(case), :][0] - T[yrs == year, ref, :][0]
    return -d if case.startswith('hist') else d


def table_a(T, yrs, cases, K, C, N):
    rows = []
    for label, case, year in LEVERS_2035:
        dg = delta_global(T, yrs, cases, case, year)
        for r in REGIONS:
            dr, mi = regional(dg, r, K)
            dh = np.concatenate([dhwf(C[(m, r, 'day')], dr[mi == i]) for i, m in enumerate(MODELS)])
            base = np.mean([np.interp(0, *C[(m, r, 'day')][:2]) for m in MODELS])
            rows.append(dict(lever=label, case=case, year=year, region=RNAME[r],
                             dT_global_median=q(dg, .5), dT_global_p05=q(dg, .05), dT_global_p95=q(dg, .95),
                             dTmax_regional_median=q(dr, .5), dTmax_p05=q(dr, .05), dTmax_p95=q(dr, .95),
                             share_members_cooling=float((dr < 0).mean()),
                             dHWF_days_median=q(dh, .5), dHWF_p05=q(dh, .05), dHWF_p95=q(dh, .95),
                             dPersonDays_million_median=q(dh, .5) * POP2035[r] / 1e6,
                             abs_dTmax_over_Tnoise=abs(q(dr, .5)) / N[r]['t'],
                             abs_dHWF_over_HWFnoise=abs(q(dh, .5)) / N[r]['h'],
                             hwf_2035_background_mean=base))
    return pd.DataFrame(rows)


def committed(C):
    rows = []
    d = pd.read_csv(SSP245DIR / 'ssp245_period_monthly_deltas.csv')
    d = d[d.target_year.eq(2035) & d.variable.eq('tasmax')].groupby(['model', 'region']).delta_c.mean()
    for r in REGIONS:
        for kind in ['day', 'night', 'compound']:
            rec = C[('recent', r, kind)]; h0 = np.interp(0, rec[0], rec[1]); s0 = np.interp(0, rec[0], rec[2])
            for m in MODELS:
                cm = C[(m, r, kind)]; h1 = np.interp(0, cm[0], cm[1])
                rows.append(dict(region=RNAME[r], kind=kind, model=m, dTmax_2035_vs_2015_2024=float(d.loc[(m, r)]),
                                 hwf_recent=h0, hwf_recent_member_sd=s0, hwf_2035=h1, dHWF=h1 - h0, dHWF_pct=100 * (h1 - h0) / h0,
                                 person_days_2035_billion=h1 * POP2035[r] / 1e9, dperson_days_billion=(h1 - h0) * POP2035[r] / 1e9))
    return pd.DataFrame(rows)


def series_b(T, yrs, cases, K, C, N):
    years = np.arange(2025, 2101)
    rows = []; toe = []
    targets = [c for c in cases if c.startswith('ssp245__') and ('__OTNZ__' in c or '__EPNZ__' in c or '__CA__' in c)] + ['ssp119__ref', 'ssp126__ref', 'ssp370__ref']
    for case in targets:
        for r in REGIONS:
            cool_med = []; cool_p17 = []; hw_med = []; hw_p17 = []
            for y in years:
                dg = delta_global(T, yrs, cases, case, y)
                dr, mi = regional(dg, r, K)
                dh = np.concatenate([dhwf(C[(m, r, 'day')], dr[mi == i]) for i, m in enumerate(MODELS)])
                rows.append(dict(case=case, region=RNAME[r], year=int(y), dTmax_median=q(dr, .5), dTmax_p05=q(dr, .05), dTmax_p17=q(dr, .17),
                                 dTmax_p83=q(dr, .83), dTmax_p95=q(dr, .95), dHWF_median=q(dh, .5), dHWF_p05=q(dh, .05), dHWF_p95=q(dh, .95),
                                 share_cooling=float((dr < 0).mean())))
                cool_med.append(-q(dr, .5)); cool_p17.append(-q(dr, .83)); hw_med.append(-q(dh, .5)); hw_p17.append(-q(dh, .83))
            cool_med = np.array(cool_med); cool_p17 = np.array(cool_p17); hw_med = np.array(hw_med); hw_p17 = np.array(hw_p17)

            def first_persistent(x, thr):
                ok = x >= thr
                for i in range(len(ok)):
                    if ok[i:].all():
                        return int(years[i])
                return None
            thr = dict(T_single=N[r]['t'], T_20yr=2 * N[r]['t'] / np.sqrt(20), H_single=N[r]['h'], H_20yr=2 * N[r]['h'] / np.sqrt(20))
            toe.append(dict(case=case, region=RNAME[r],
                            crossover_median_cooling=first_persistent(cool_med, 1e-12),
                            ToE_T20yr_median=first_persistent(cool_med, thr['T_20yr']), ToE_T20yr_likely=first_persistent(cool_p17, thr['T_20yr']),
                            ToE_Tsingle_median=first_persistent(cool_med, thr['T_single']),
                            ToE_H20yr_median=first_persistent(hw_med, thr['H_20yr']), ToE_H20yr_likely=first_persistent(hw_p17, thr['H_20yr']),
                            ToE_Hsingle_median=first_persistent(hw_med, thr['H_single']),
                            cooling_2035=cool_med[years == 2035][0], cooling_2050=cool_med[years == 2050][0], cooling_2060=cool_med[years == 2060][0],
                            cooling_2100=cool_med[years == 2100][0], **{f'thr_{k}': v for k, v in thr.items()}))
    return pd.DataFrame(rows), pd.DataFrame(toe)


def self_capture(T, yrs, cases):
    rows = []
    for st in ['OTNZ', 'EPNZ']:
        for ss in ['co2only', 'all', 'allCH4']:
            for y in [2035, 2050, 2060, 2100]:
                own = delta_global(T, yrs, cases, f'ssp245__own_share__{st}__{ss}', y)
                nat = delta_global(T, yrs, cases, f'ssp245__national__{st}__{ss}', y)
                uni = delta_global(T, yrs, cases, f'ssp245__universal__{st}__{ss}', y)
                rows.append(dict(stringency=st, species=ss, year=y, own_median=q(own, .5), national_median=q(nat, .5), universal_median=q(uni, .5),
                                 own_over_universal_pct=100 * q(own, .5) / q(uni, .5), national_over_universal_pct=100 * q(nat, .5) / q(uni, .5)))
    return pd.DataFrame(rows)


def main():
    OUTA.mkdir(parents=True, exist_ok=True); OUTB.mkdir(parents=True, exist_ok=True)
    T, yrs, cases, A = load_fair(); K = kappas(); C = curves(); N = noise()
    ta = table_a(T, yrs, cases, K, C, N); ta.to_csv(OUTA / 'A1_2035_scale_decomposition.csv', index=False, encoding='utf-8-sig')
    cm = committed(C); cm.to_csv(OUTA / 'A2_committed_exposure_2035.csv', index=False, encoding='utf-8-sig')
    sb, toe = series_b(T, yrs, cases, K, C, N)
    sb.to_csv(OUTB / 'B1_regional_timeseries_2025_2100.csv', index=False, encoding='utf-8-sig')
    toe.to_csv(OUTB / 'B2_time_of_emergence.csv', index=False, encoding='utf-8-sig')
    sc = self_capture(T, yrs, cases); sc.to_csv(OUTB / 'B3_self_capture_ratio.csv', index=False, encoding='utf-8-sig')
    save(EXP / 'scales' / 'analysis_inputs.json', dict(kappa_tasmax={r: dict(zip(MODELS, K[r].tolist())) for r in REGIONS}, noise=N, pop2035=POP2035))
    pd.set_option('display.width', 250)
    print(ta[['lever', 'region', 'dTmax_regional_median', 'dTmax_p05', 'dTmax_p95', 'dHWF_days_median', 'abs_dHWF_over_HWFnoise']].round(5).to_string())
    print(cm.round(2).to_string())
    print(toe.to_string())
    print(sc.round(4).to_string())


if __name__ == '__main__':
    main()
