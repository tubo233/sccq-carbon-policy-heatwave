"""补充表 S2–S6（CSV）和图 S1–S3（PNG）。

图 S5 只在本地放有 CNEMC 派生汇总时生成，不随包分发；表 S1、S7、S8 和图 S4 在别处生成。
"""
import importlib.util
import json
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
EXP = HERE / 'models'
D = EXP / 'main'
HW = EXP / 'heatwave'
C = EXP / 'air_quality'
SI = HERE / 'data' / 'si'; SI.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location('pf', HERE / 'scripts/03_paper_figures.py')
pf = importlib.util.module_from_spec(spec); spec.loader.exec_module(pf)
REG = {'sichuan': 'Sichuan', 'chongqing': 'Chongqing'}
CN = {'sichuan': '四川', 'chongqing': '重庆'}
MODELS = ['GFDL-ESM4', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
# Okabe-Ito 配色，与主图一致（四川蓝、重庆朱红）
BLUE, ORANGE, GREEN, VERMILLION, PURPLE = '#0072B2', '#E69F00', '#009E73', '#D55E00', '#CC79A7'
GEN_COL, OBS_COL = BLUE, ORANGE


# ---------------------------------------------------------------- 表格
def table_s2():
    h = pd.read_csv(D / 'history_gdp_energy_co2.csv')
    rows = []
    for y in range(2005, 2025):
        r = {'Year': y}
        for reg, t in [('sichuan', 'SC'), ('chongqing', 'CQ')]:
            x = h[(h.region == reg) & (h.year == y)].iloc[0]
            r[f'{t} GDP index (1978 = 100)'] = '' if pd.isna(x.gdp_index_1978_100) else f'{x.gdp_index_1978_100:.0f}'
            e = '' if pd.isna(x.energy_10k_tce_equal) else f'{x.energy_10k_tce_equal / 100:.1f}'
            if e and reg == 'chongqing' and x.energy_equal_status != 'official':
                e += ' d'
            r[f'{t} energy (Mtce)'] = e
            r[f'{t} CO₂ (Mt)'] = '' if pd.isna(x.co2_mt_ceads) else f'{x.co2_mt_ceads:.0f}'
        rows.append(r)
    pd.DataFrame(rows).to_csv(SI / 'tableS2_energy_gdp_co2.csv', index=False, encoding='utf-8-sig')


def table_s3():
    lam = pd.read_csv(HERE / 'data/lambda_values.csv')
    w = {'years_41_60': 'years 41–60', 'years_61_80_primary': 'years 61–80', 'years_81_100': 'years 81–100'}
    lam['window'] = lam.window.map(lambda s: w.get(s, s))
    t = lam.pivot_table(index=['model', 'window'], columns='region', values='lam', sort=False).reset_index()
    t = t[['model', 'window', 'sichuan', 'chongqing']]
    t.columns = ['Source', 'Window', 'λ Sichuan', 'λ Chongqing']
    for c in ['λ Sichuan', 'λ Chongqing']:
        t[c] = t[c].map(lambda v: f'{v:.3f}')
    t.to_csv(SI / 'tableS3_lambda.csv', index=False, encoding='utf-8-sig')


NPI_LABEL = {'NPI_aeei0.5': 'AEEI 0.5% yr⁻¹', 'NPI_aeei1.0': 'AEEI 1% yr⁻¹ (main)', 'NPI_aeei1.5': 'AEEI 1.5% yr⁻¹',
             'NPI_aeei2.0': 'AEEI 2% yr⁻¹', 'NPI_aeei3.0': 'AEEI 3% yr⁻¹', 'NPI_aeei4.0': 'AEEI 4% yr⁻¹', 'NPI_aeei5.0': 'AEEI 5% yr⁻¹',
             'NPI_aeei1.0_cfobs': 'AEEI 1% yr⁻¹, CF observed to 2022', 'NPI_aeei4.0_cfobs': 'AEEI 4% yr⁻¹, CF observed to 2022',
             'NPI_aeei1.0_cqcf0.9': 'AEEI 1% yr⁻¹, CQ 2005 CF −10%', 'NPI_aeei1.0_cqcf1.1': 'AEEI 1% yr⁻¹, CQ 2005 CF +10%',
             'NPI_trend1997': 'CI trend 1997–2005', 'NPI_trend2000': 'CI trend 2000–2005',
             'NPI_trend2002': 'CI trend 2002–2005', 'NPI_from2023': 'Constraints removed from 2023',
             'NPI_aeei1.0_lowGDP': 'AEEI 1% yr⁻¹, GDP −1 pp from 2026'}
NPI_ORDER = ['NPI_aeei1.0', 'NPI_aeei0.5', 'NPI_aeei1.5', 'NPI_aeei2.0', 'NPI_aeei3.0', 'NPI_aeei4.0', 'NPI_aeei5.0',
             'NPI_aeei1.0_cfobs', 'NPI_aeei4.0_cfobs', 'NPI_aeei1.0_cqcf0.9', 'NPI_aeei1.0_cqcf1.1',
             'NPI_trend1997', 'NPI_trend2000', 'NPI_trend2002', 'NPI_from2023']
NPLS_LABEL = {'NPLS_gcam_mid': 'GCAM mid (main)', 'NPLS_gcam_low': 'GCAM low', 'NPLS_gcam_high': 'GCAM high',
              'NPLS_provtarget': 'Provincial targets', 'NPLS_provtarget_lowGDP(NPR_lowGDP)': 'Provincial targets, GDP −1 pp'}


def table_s4():
    sp = pd.read_csv(D / 'scenario_paths_energy_co2.csv')
    rows = []
    for s in NPI_ORDER:
        r = {'NPI construction': NPI_LABEL[s]}
        for reg, t in [('sichuan', 'SC'), ('chongqing', 'CQ')]:
            d = sp[(sp.region == reg)]
            a = d[d.scenario == s].set_index('year').co2_mt; o = d[d.scenario == 'NPR'].set_index('year').co2_mt
            m = re.search(r'= ([+-−]?\d+\.\d+)%/yr', d[d.scenario == s].definition.iloc[0])
            r[f'{t} CI trend (% yr⁻¹)'] = m.group(1).replace('-', '−') if m else '—'
            r[f'{t} excess 2006–2022 (Gt CO₂)'] = f'{(a.loc[2006:2022] - o.loc[2006:2022]).sum() / 1e3:.2f}'
            r[f'{t} CO₂ 2035 (Mt)'] = f'{a[2035]:.0f}'
        rows.append(r)
    pd.DataFrame(rows).to_csv(SI / 'tableS4_npi_constructions.csv', index=False, encoding='utf-8-sig')


def table_s5():
    """单因素表：每个 NPI 构造配主 NPLS，每个 NPLS 变体配主 NPI，再加低 GDP 情形。"""
    m = pd.read_csv(D / 'results/metrics_construction_sensitivity.csv', encoding='utf-8-sig')
    combos = [(n, 'NPLS_gcam_mid') for n in NPI_ORDER] + [('NPI_aeei1.0', s) for s in ['NPLS_gcam_low', 'NPLS_gcam_high', 'NPLS_provtarget']]
    combos.append(('NPI_aeei1.0_lowGDP', 'NPLS_provtarget_lowGDP(NPR_lowGDP)'))
    rows = []
    for npi, npls in combos:
        g = m[(m.NPI == npi) & (m.NPLS == npls)]
        s = g[g.region == '四川'].set_index('year'); c = g[g.region == '重庆'].set_index('year')
        rows.append({'NPI': NPI_LABEL.get(npi, npi), 'NPLS': NPLS_LABEL.get(npls, npls),
                     'SC A 2030 (°C)': f'{s.A_NPI_NPLS[2030]:.4f}', 'SC A 2035 (°C)': f'{s.A_NPI_NPLS[2035]:.4f}',
                     'SC A 2060 (°C)': f'{s.A_NPI_NPLS[2060]:.4f}', 'CQ A 2035 (°C)': f'{c.A_NPI_NPLS[2035]:.4f}',
                     'SC ΔHWF 2035 (days yr⁻¹)': f'{s.dHWF_NPI_NPLS[2035]:.3f}', 'CQ ΔHWF 2035 (days yr⁻¹)': f'{c.dHWF_NPI_NPLS[2035]:.3f}',
                     'NPR − NPLS 2035 SC (°C)': f'{s.A_NPR_NPLS[2035]:.5f}', 'P 2035 (%)': f'{s.P_NPLS[2035]:.2f}',
                     'C 2035 (%)': f'{s.C_NPLS[2035]:.1f}', 'G 2035 (%)': f'{s.G_local_gain[2035]:.1f}'})
    pd.DataFrame(rows).to_csv(SI / 'tableS5_construction_sensitivity.csv', index=False, encoding='utf-8-sig')


def table_s10():
    """历史时段 NPI − NPR 的均值差（两省同时改情景）。"""
    sp = pd.read_csv(D / 'scenario_paths_energy_co2.csv')
    ts = pd.read_csv(HERE / 'data/timeseries_metrics_1990_2060.csv')
    rows = []
    for a, b in [(2006, 2010), (2011, 2015), (2016, 2022)]:
        r = {'Period': f'{a}–{b}'}
        for reg, t in [('sichuan', 'SC'), ('chongqing', 'CQ')]:
            d = sp[sp.region == reg]
            npi = d[d.scenario == 'NPI_aeei1.0'].set_index('year').co2_mt; npr = d[d.scenario == 'NPR'].set_index('year').co2_mt
            x = ts[(ts.region == reg) & (ts.view == 'combined') & ts.year.between(a, b)]
            r[f'{t} excess CO₂ (Mt yr⁻¹)'] = f'{(npi.loc[a:b] - npr.loc[a:b]).mean():.0f}'
            r[f'{t} A, NPI − NPR (°C)'] = f"{x[x.metric == 'A_NPI_NPR']['median'].mean():.5f}"
            r[f'{t} ΔHWF (days yr⁻¹)'] = f"{x[x.metric == 'dHWF_NPI_NPR']['median'].mean():.3f}"
        rows.append(r)
    pd.DataFrame(rows).to_csv(SI / 'tableS10_historical_periods.csv', index=False, encoding='utf-8-sig')


def table_s11():
    """单省情景（只改一个省）下各 NPI 构造的结果。"""
    o = pd.read_csv(D / 'results/metrics_construction_sensitivity_by_province.csv', encoding='utf-8-sig')
    o = o[(o.year == 2035) & (o.NPLS == 'NPLS_gcam_mid')]
    rows = []
    for npi in NPI_ORDER:
        r = {'NPI construction': NPI_LABEL[npi]}
        for tag in ['SC', 'CQ']:
            x = o[(o.NPI == npi) & (o.view == f'{tag} only')].iloc[0]
            r[f'{tag} only: A, NPI − NPLS (°C)'] = f'{x.A_NPI_NPLS:.4f}'.replace('-', '−')
            r[f'{tag} only: C (%)'] = f'{x.C_NPLS:.1f}'.replace('-', '−')
        rows.append(r)
    pd.DataFrame(rows).to_csv(SI / 'tableS11_province_specific_constructions.csv', index=False, encoding='utf-8-sig')


def table_s6():
    c = pd.read_csv(HW / 'response_curve.csv')
    c = c[(c.definition == 'annual90') & (c.offset_tx == 0)]
    rows = []
    for w, wl in [('population2020', 'Population'), ('area', 'Area')]:
        for k, kl in [('day', 'Daytime'), ('night', 'Night-time'), ('compound', 'Compound')]:
            for mdl in MODELS:
                r = {'Weighting': wl, 'Heatwave type': kl, 'CMIP6 model': mdl}
                for reg, t in [('sichuan', 'SC'), ('chongqing', 'CQ')]:
                    x = c[(c.weighting == w) & (c.kind == k) & (c.region == reg)].set_index('background').hwf_mean
                    r[f'{t} 2015–2024 (days)'] = f'{x["recent"]:.1f}'
                    r[f'{t} 2035 (days)'] = f'{x[mdl]:.1f}'
                    r[f'{t} change (days; %)'] = f'+{x[mdl] - x["recent"]:.1f} ({100 * (x[mdl] / x["recent"] - 1):+.0f}%)'
                rows.append(r)
    pd.DataFrame(rows).to_csv(SI / 'tableS6_committed_exposure.csv', index=False, encoding='utf-8-sig')


# ---------------------------------------------------------------- 图
def observed_hwn():
    """逐年观测人口加权 HWF 和 HWN（格点识别，固定 1990–2005 Q90，≥3 天）。"""
    LEGACY = HERE / 'weather'
    sp = importlib.util.spec_from_file_location('legacy_common', LEGACY / 'scripts/common.py')
    old = importlib.util.module_from_spec(sp); sp.loader.exec_module(old)
    dates, meta, tx, tn = old.load_cache(); tx = np.asarray(tx)
    qx, _ = old.thresholds(tx, tx, dates, 'annual90')
    years, hwf, hwn = old.annual_metrics(tx > qx[None, :], dates)
    out = []
    for r in REG:
        k = meta.region_id.eq(r).to_numpy(); w = meta.population2020.to_numpy()[k]; w = w / w.sum()
        for i, y in enumerate(years):
            out.append(dict(region=r, year=int(y), hwf=float(hwf[i, k] @ w), hwn=float(hwn[i, k] @ w)))
    o = pd.DataFrame(out); o.to_csv(SI / 'observed_hwf_hwn_popweighted.csv', index=False)
    return o


def fig_s1():
    z = np.load(HW / 'response_members.npz', allow_pickle=True)
    cases = list(z['cases']); V = z['values']
    ci = cases.index('recent__+0.00')
    idx = lambda arr, v: list(z[arr]).index(v)
    obs_file = SI / 'observed_hwf_hwn_popweighted.csv'
    o = pd.read_csv(obs_file) if obs_file.exists() else observed_hwn()
    # 第三栏放观测年序列，顺便看 2022 在年份里的位置
    fig = plt.figure(figsize=(7.1, 5.2))
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.98, bottom=0.09, top=0.86, wspace=0.3, hspace=0.5)
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]
    stats = {}
    for j, (reg, name) in enumerate(REG.items()):
        ri = idx('regions', reg)
        sel = lambda metric: V[ci, :, ri, idx('definitions', 'annual90'), idx('kinds', 'day'), idx('weightings', 'population2020'), idx('metrics', metric)]
        mf, mn = sel('HWF'), sel('HWN')
        d = o[o.region == reg].sort_values('year')
        yv = d.year.to_numpy(); f = d.hwf.to_numpy()
        det = f - np.polyval(np.polyfit(yv, f, 1), yv) + d[(d.year >= 2015)].hwf.mean()
        dur_obs = (d.hwf / d.hwn.replace(0, np.nan)).to_numpy()
        stats[reg] = dict(member_sd=float(mf.std(ddof=1)), obs_detrended_sd=float(det.std(ddof=2)),  # 回归残差 SD，同 observed_noise_scales.csv
                          member_mean=float(mf.mean()), obs_2015_2024_mean=float(d[d.year >= 2015].hwf.mean()),
                          member_duration=float(np.nanmean(mf / np.where(mn > 0, mn, np.nan))),
                          obs_duration_1990_2024=float(np.nanmean(dur_obs)), obs_duration_2015_2024=float(np.nanmean(dur_obs[yv >= 2015])))
        # (a) 年 HWF 分布
        ax = axs[0]; pf.style(ax)
        pos = 2 * j
        for p, vals, col, lab in [(pos, mf, GEN_COL, 'Weather generator (128 members)'), (pos + 0.8, det, OBS_COL, 'Observed 1990–2024, detrended')]:
            b = ax.boxplot(vals, positions=[p], widths=0.55, patch_artist=True, showfliers=True,
                           medianprops=dict(color='white', lw=1.2), flierprops=dict(marker='o', ms=2, mfc=col, mec=col))
            for x in b['boxes']:
                x.set(facecolor=col, edgecolor=col, alpha=0.85)
            for x in b['whiskers'] + b['caps']:
                x.set(color=col, lw=0.8)
        # (b) 平均持续时间
        ax = axs[1]; pf.style(ax)
        ax.bar([pos], [stats[reg]['member_duration']], width=0.55, color=GEN_COL)
        ax.bar([pos + 0.8], [stats[reg]['obs_duration_2015_2024']], width=0.55, color=OBS_COL)
    axs[0].set_xticks([0.4, 2.4]); axs[0].set_xticklabels(list(REG.values()))
    axs[0].set_ylabel('Annual heatwave days', fontsize=7.5, color=pf.INK2)
    axs[1].set_xticks([0.4, 2.4]); axs[1].set_xticklabels(list(REG.values()))
    axs[1].set_ylabel('Mean heatwave duration (days)', fontsize=7.5, color=pf.INK2)
    fig.legend(handles=[Line2D([], [], color=GEN_COL, lw=6, label='Generated (recent climate)'),
                           Line2D([], [], color=OBS_COL, lw=6, label='Observed (detrended to 2015–2024)')],
                  frameon=False, fontsize=6.5, loc='upper left', bbox_to_anchor=(0.06, 1.0), ncol=2)
    axs[1].text(0.05, 1.04, 'observed: 2015–2024 mean', transform=axs[1].transAxes, fontsize=6.3, color=pf.INK2)
    ax = axs[2]; pf.style(ax)
    for reg, col in [('sichuan', BLUE), ('chongqing', VERMILLION)]:
        d = o[o.region == reg].sort_values('year')
        ax.plot(d.year, d.hwf, color=col, lw=1.3, marker='o', ms=2.5, label=REG[reg])
    for yr in (2006, 2022, 2024):
        v = float(o[o.year == yr].hwf.max())
        ax.annotate(str(yr), (yr, v), xytext=(0, 5), textcoords='offset points', ha='center', fontsize=6.3, color=pf.INK2)
    ax.axvspan(2014.5, 2024.5, color='#f6f6f3', lw=0, zorder=0)
    ax.set_xlim(1989.5, 2024.5); ax.set_ylabel('Observed heatwave days', fontsize=7.5, color=pf.INK2)
    ax.legend(frameon=False, fontsize=6.6, loc='upper left')
    pf.tag(axs[0], 'a', x=0.0, y=1.02); pf.tag(axs[1], 'b', x=0.0, y=1.02); pf.tag(axs[2], 'c', x=0.0, y=1.02)
    pf.save(fig, 'FigureS1_weather_generator_validation')
    json.dump(stats, open(SI / 'figureS1_stats.json', 'w'), indent=1)
    return stats


def fig_s2():
    c = pd.read_csv(HW / 'response_curve.csv')
    c = c[(c.definition == 'annual90') & (c.weighting == 'population2020')]
    fig, axs = plt.subplots(2, 3, figsize=(7.1, 4.4), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.1, top=0.81, wspace=0.28, hspace=0.25)
    cols = {'recent': pf.INK, 'GFDL-ESM4': ORANGE, 'MPI-ESM1-2-LR': GREEN, 'UKESM1-0-LL': PURPLE}
    for i, reg in enumerate(REG):
        for j, (k, kl) in enumerate([('day', 'Daytime'), ('night', 'Night-time'), ('compound', 'Compound')]):
            ax = axs[i, j]; pf.style(ax)
            for bg, col in cols.items():
                d = c[(c.region == reg) & (c.kind == k) & (c.background == bg)].sort_values('offset_tx')
                if bg == 'recent':  # 2035 各背景的范围宽度接近，只画一条带子免得线太乱
                    ax.fill_between(d.offset_tx, d.hwf_q025, d.hwf_q975, color=col, alpha=0.10, lw=0)
                ax.plot(d.offset_tx, d.hwf_mean, color=col, lw=1.5 if bg == 'recent' else 1.3,
                        label='2015–2024 climate' if bg == 'recent' else f'2035 background, {bg}')
            ax.set_title(f'{REG[reg]}, {kl.lower()}', fontsize=7.5, color=pf.INK, loc='left')
            if i == 1:
                ax.set_xlabel('Uniform Tmax shift δ (°C)', fontsize=7.2, color=pf.INK2)
            if j == 0:
                ax.set_ylabel('Heatwave days yr$^{-1}$', fontsize=7.2, color=pf.INK2)
    fig.legend(*axs[0, 0].get_legend_handles_labels(), frameon=False, fontsize=8, loc='upper left', bbox_to_anchor=(0.06, 1.0), ncol=2)
    pf.save(fig, 'FigureS2_response_curves')


def fig_s3():
    ts = pd.read_csv(HERE / 'data/timeseries_metrics_1990_2060.csv')
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.8))
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.16, top=0.86, wspace=0.28)
    for reg, col in [('sichuan', BLUE), ('chongqing', VERMILLION)]:
        for ax, metric, lab in [(axs[0], 'A_NPI_NPLS', 'Avoided warming NPI − NPLS (°C)'), (axs[1], 'C_NPLS', 'Own-emission warming avoided, C (%)')]:
            for view, ls in [('combined', (0, (4, 2))), ('own', '-')]:
                d = ts[(ts.region == reg) & (ts.view == view) & (ts.metric == metric) & (ts.year >= 2010)].sort_values('year')
                if view == 'own':
                    ax.fill_between(d.year, d.p05, d.p95, color=col, alpha=0.12, lw=0)
                ax.plot(d.year, d['median'], color=col, lw=1.4 if view == 'own' else 1.0, ls=ls)
            pf.style(ax); ax.set_ylabel(lab, fontsize=7.2, color=pf.INK2)
            ax.axvline(2022.5, color=pf.AXIS, lw=0.6); ax.axvspan(2035.5, 2060, color='#f6f6f3', lw=0, zorder=0)
            ax.set_xlim(2010, 2060)
    axs[0].legend(handles=[Line2D([], [], color=BLUE, lw=1.4, label='Sichuan'), Line2D([], [], color=VERMILLION, lw=1.4, label='Chongqing'),
                           Line2D([], [], color=pf.INK2, lw=1.4, label='only this province changes'),
                           Line2D([], [], color=pf.INK2, lw=1.0, ls=(0, (4, 2)), label='both provinces change')],
                  frameon=False, fontsize=6.3, loc='upper left')
    pf.tag(axs[0], 'a', x=0.0); pf.tag(axs[1], 'b', x=0.0)
    pf.save(fig, 'FigureS3_province_specific_runs')


def fig_s5():
    """只做描述：PM2.5、臭氧及其与热浪的同期出现。"""
    src = C / 'heatwave_day_pollution_by_year.csv'
    if not src.exists():
        print('跳过图 S5：CNEMC 派生汇总未随包提供')
        return
    h = pd.read_csv(src)
    fig = plt.figure(figsize=(7.1, 2.6))
    axs = [fig.add_axes([0.075 + i * 0.33, 0.14, 0.235, 0.72]) for i in range(3)]
    for ax, (a, b, lab) in zip(axs[:2], [('pm25_on_hw', 'pm25_all', 'PM$_{2.5}$ (μg m$^{-3}$)'), ('mda8_on_hw', 'mda8_all', 'O$_3$ MDA8 (μg m$^{-3}$)')]):
        pf.style(ax)
        ax.plot(h.year, h[a], color=VERMILLION, lw=1.4, marker='o', ms=2.5, label='heatwave days')
        ax.plot(h.year, h[b], color=pf.INK2, lw=1.2, marker='o', ms=2.5, label='all days')
        ax.set_ylabel(lab, fontsize=7.2, color=pf.INK2); ax.set_xticks([2015, 2018, 2021, 2024])
    axs[1].axhline(160, color=pf.INK2, lw=0.6, ls=(0, (3, 2)))
    axs[1].text(2015, 162, 'Grade II standard (160)', fontsize=6, color=pf.INK2, va='bottom'); axs[1].set_ylim(70, 172)
    axs[0].legend(frameon=False, fontsize=6.4, loc='upper right')
    ax = axs[2]; pf.style(ax)
    ax.bar(h.year, 100 * h.compound_share_of_hw, color=BLUE, width=0.7)
    ax.set_ylabel('Heatwave person-days with\nO$_3$ exceedance (%)', fontsize=7.2, color=pf.INK2); ax.set_xticks([2015, 2018, 2021, 2024])
    for i, t in enumerate('abc'):
        pf.tag(axs[i], t, x=0.0)
    pf.save(fig, 'FigureS5_air_quality_heatwave_days')


if __name__ == '__main__':
    todo = sys.argv[1:] or ['table_s2', 'table_s3', 'table_s4', 'table_s5', 'table_s6', 'table_s10', 'table_s11', 'fig_s1', 'fig_s2', 'fig_s3', 'fig_s5']
    for f in todo:
        print(f, globals()[f]() or '')
