"""按 Kaya 恒等式生成川渝 NPI / NPR / NPLS 的能源与 CO2 路径（2006-2060）。

GDP 增速按年鉴和已有预测接续；NPI 用 AEEI 1%（敏感性 0.5–5%）并冻结 CF，NPR 用国家强度目标，
NPLS 用 GCAM 中情景或省级目标。具体构造和敏感性见 README / Table S1C。
"""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
YEARS = np.arange(2005, 2061)
LAB = ROOT / 'inputs/gcam/gcam_china_projections_2025_2060.csv'


def gdp_path(h, reg, shift=0.0):
    g = h.set_index('year').gdp_index_1978_100.dropna()
    idx = {y: g[y] for y in g.index if y >= 2005}
    real = {'sichuan': {2023: .060, 2024: .057}, 'chongqing': {}}[reg]
    last = max(idx)
    rates = {}
    for y in range(last + 1, 2061):
        if y in real:
            r = real[y]
        elif y == 2025:
            r = .050
        elif y <= 2030:
            r = {'sichuan': .0647, 'chongqing': .0510}[reg]
        elif y <= 2035:
            r = {'sichuan': .0570, 'chongqing': .0485}[reg]
        elif y <= 2045:
            r = .04
        else:
            r = .03
        rates[y] = r
    out = {}
    # 低增长敏感性只从 2026 年起下调；2025 年及以前保留年鉴指数与实际/官方增速
    for y in YEARS:
        if y in idx and (shift == 0 or y <= 2025):
            out[y] = idx[y]
        else:
            prev = out[y - 1]
            r = rates.get(y, (idx[y] / idx[y - 1] - 1) if y in idx else .05)
            out[y] = prev * (1 + r - (shift if y >= 2026 else 0))
    return pd.Series(out)


def lab_paths():
    # GCAM-China 2025-2060 预测（四部门），已衔接到 CEADs 2020 年排放
    x = pd.read_csv(LAB)
    return {(reg, lab): g.groupby('year').co2_mt.sum() for (reg, lab), g in x.groupby(['region', 'scenario'])}


_dpec = pd.read_csv(EMIS / 'DPEC7383239_CO2.csv').set_index('year')['on-time peak-net zero-clean air']


def peak_and_ndc(c, cut2035=.085):
    """国家约束：2030 年后不超过 2030 年水平；2035 年比峰值低 cut2035（国家 2035 NDC：比峰值降 7-10%）；
    2036 年后按 DPEC 全国 on-time peak-net zero 路径的相对变化外延（仅作延伸）。"""
    peak = max(c[y] for y in range(2023, 2031))
    y_peak = max((y for y in range(2023, 2031)), key=lambda y: c[y])
    target35 = peak * (1 - cut2035)
    out = dict(c)
    for y in range(y_peak + 1, 2036):
        lin = peak + (target35 - peak) * (y - y_peak) / (2035 - y_peak)
        out[y] = min(c[y], lin)
    for y in range(2036, 2061):
        out[y] = out[2035] * _dpec[y] / _dpec[2035]
    return out


def interp_to(series, years):
    return np.interp(years, series.index.to_numpy(float), series.to_numpy(float))


def build(shift=0.0):
    h = pd.read_csv(HERE / 'history_gdp_energy_co2.csv')
    labs = lab_paths()
    rows = []; meta = {}
    for reg in ['sichuan', 'chongqing']:
        hr = h[h.region.eq(reg)].set_index('year')
        G = gdp_path(hr.reset_index(), reg, shift)
        E_obs = hr.energy_10k_tce_equal; C_obs = hr.co2_mt_ceads
        EI05 = E_obs[2005] / G[2005]; CF05 = C_obs[2005] / E_obs[2005]
        CI = (C_obs / hr.gdp_index_1978_100).dropna()
        paths = {}
        # --- NPR（国家政策响应）---
        EI20, CI20 = E_obs[2020] / G[2020], C_obs[2020] / G[2020]
        EI22, CI22 = E_obs[2022] / G[2022], C_obs[2022] / G[2022]

        def target_path(x22, x25, rate5, years):
            out = {}
            for y in years:
                if y <= 2022:
                    continue
                if y <= 2025:
                    x25b = min(x25, x22 * (1 - .01) ** 3)          # 若已提前达标，至少再降1%/年
                    out[y] = x22 * (x25b / x22) ** ((y - 2022) / 3)
                else:
                    out[y] = out[2025] * (1 - rate5) ** ((y - 2025) / 5)
            return out
        ei_npr = target_path(EI22, EI20 * (1 - .135), .135, YEARS)
        ci_npr = target_path(CI22, CI20 * (1 - .18), .173, YEARS)
        e_npr = {y: (E_obs[y] if y <= 2022 else G[y] * ei_npr[y]) for y in YEARS}
        c_npr = peak_and_ndc({y: (C_obs[y] if y <= 2022 else G[y] * ci_npr[y]) for y in YEARS})
        paths['NPR'] = (e_npr, c_npr, 'observed to 2022; national average intensity targets; national 2030 peak and 2035 NDC (-8.5% from peak); DPEC on-time net-zero after 2035')
        # --- NPI：AEEI（自主技术进步）---
        # 5%/年 的上限覆盖了 2005-2022 年能耗强度实际下降幅度（4.2-4.3%/年）
        for a in [.005, .01, .015, .02, .03, .04, .05]:
            e = {y: G[y] * EI05 * (1 - a) ** (y - 2005) for y in YEARS}
            c = {y: e[y] * CF05 for y in YEARS}
            paths[f'NPI_aeei{a * 100:.1f}'] = (e, c, f'no intensity policy after 2005; AEEI {a * 100:.1f}%/yr; CF frozen 2005')
        # 碳强度因子按观测值走到 2022 年后冻结，这样 2005 年前已规划的水电等历史减碳
        # 都算到非政策因素
        CF_obs = (C_obs / E_obs).dropna()
        for a in [.01, .04]:
            e = {y: G[y] * EI05 * (1 - a) ** (y - 2005) for y in YEARS}
            c = {y: e[y] * (CF_obs[y] if y <= 2022 else CF_obs[2022]) for y in YEARS}
            paths[f'NPI_aeei{a * 100:.1f}_cfobs'] = (e, c, f'no intensity policy after 2005; AEEI {a * 100:.1f}%/yr; CF observed to 2022, held after')
        # 重庆 2005 年碳强度因子敏感性（能源核算口径）；四川不变
        for f_cf in [0.9, 1.1]:
            k = f_cf if reg == 'chongqing' else 1.0
            e = {y: G[y] * EI05 * .99 ** (y - 2005) for y in YEARS}
            c = {y: e[y] * CF05 * k for y in YEARS}
            paths[f'NPI_aeei1.0_cqcf{f_cf:.1f}'] = (e, c, f'as NPI_aeei1.0; Chongqing 2005 carbon factor x{f_cf:.1f}')
        # NPI-2023：2022 前观测，2023 起 AEEI 1%，CF 冻结在 2022
        EI22_, CF22 = E_obs[2022] / G[2022], C_obs[2022] / E_obs[2022]
        e = {y: (E_obs[y] if y <= 2022 else G[y] * EI22_ * .99 ** (y - 2022)) for y in YEARS}
        c = {y: (C_obs[y] if y <= 2022 else e[y] * CF22) for y in YEARS}
        paths['NPI_from2023'] = (e, c, 'observed to 2022; no intensity policy from 2023; AEEI 1%/yr')
        # NPI 统计趋势窗口
        for w0 in [1997, 2000, 2002]:
            b = np.polyfit(np.arange(w0, 2006), np.log(CI.loc[w0:2005].to_numpy()), 1)[0]
            c = {y: G[y] * CI[2005] * np.exp(b * (y - 2005)) for y in YEARS}
            e = {y: c[y] / CF05 for y in YEARS}
            paths[f'NPI_trend{w0}'] = (e, c, f'carbon-intensity log-linear trend {w0}-2005 = {100 * (np.exp(b) - 1):+.2f}%/yr')
        # --- NPLS：GCAM 实验室情景 ---
        ei_pl = target_path(EI22, EI20 * (1 - .14), .14, YEARS)
        e_pl = {y: (E_obs[y] if y <= 2022 else G[y] * ei_pl[y]) for y in YEARS}
        for lab in ['mid', 'low', 'high']:
            s = labs[(reg, lab)]
            sy = pd.Series(interp_to(s, YEARS), index=YEARS)
            c = {}
            for y in YEARS:
                if y <= 2022:
                    c[y] = C_obs[y]
                elif y < 2025:
                    c[y] = C_obs[2022] + (sy[2025] - C_obs[2022]) * (y - 2022) / 3
                else:
                    c[y] = sy[y]
            paths[f'NPLS_gcam_{lab}'] = (e_pl, c, f'lab GCAM-China {lab} scenario (4-sector CO2), connected to CEADs 2020')
        # --- NPLS：省级目标 ---
        ci_pl = target_path(CI22, CI20 * (1 - .195), .173, YEARS)
        if reg == 'sichuan':
            cap30 = CI[2005] * .30
            for y in range(2026, 2061):
                lim = ci_pl[2025] * (cap30 / ci_pl[2025]) ** ((y - 2025) / 5) if y <= 2030 else None
                if lim is not None:
                    ci_pl[y] = min(ci_pl[y], lim)
            for y in range(2031, 2061):
                ci_pl[y] = ci_pl[2030] * (1 - .173) ** ((y - 2030) / 5)
        c = peak_and_ndc({y: (C_obs[y] if y <= 2022 else G[y] * ci_pl[y]) for y in YEARS})
        paths['NPLS_provtarget'] =(e_pl, c, 'provincial documents: EI -14% (2021-25), CI -19.5% (assumed allocation), SC 2030 CI -70% vs 2005; stricter of local vs national')
        for name, (e, c, desc) in paths.items():
            for y in YEARS:
                rows.append(dict(region=reg, scenario=name, year=int(y), gdp_index=G[y], energy_10k_tce=e[y], co2_mt=c[y],
                                 observed=(y <= 2022 and name in ('NPR',) or (y <= 2022 and name.startswith('NPLS')) or (y <= 2022 and name == 'NPI_from2023')),
                                 definition=desc))
        meta[reg] = dict(EI2005=EI05, CF2005_tco2_per_tce=CF05 * 100, CI2005=CI[2005], CI2020=CI20, CI2022=CI22,
                         CI_2022_vs_2020_pct=100 * (CI22 / CI20 - 1), EI_2022_vs_2020_pct=100 * (EI22 / EI20 - 1))
    return pd.DataFrame(rows), meta


def main():
    df, meta = build(0.0)
    df.to_csv(HERE / 'scenario_paths_energy_co2.csv', index=False)
    low, _ = build(0.01)
    low.to_csv(HERE / 'scenario_paths_energy_co2_lowGDP.csv', index=False)
    save(HERE / 'scenario_meta.json', meta)
    p = df[df.year.isin([2005, 2020, 2025, 2030, 2035, 2060])].pivot_table(index=['region', 'scenario'], columns='year', values='co2_mt').round(0)
    print(p.to_string())
    cum = df[df.year.between(2006, 2035)].groupby(['region', 'scenario']).co2_mt.sum().unstack(0).round(0)
    print('cumulative 2006-2035 Mt'); print(cum.to_string())
    print(json.dumps(meta, indent=1, ensure_ascii=False)[:1200])


if __name__ == '__main__':
    main()
