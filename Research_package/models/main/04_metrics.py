"""政策情景指标与 FaIR × λ_r 敏感性。

W_s = λ_r [T_s − mean_{2015-2024} T_NPR]；A/P/G/C 分别是避免增温、占 NPI 增温比例、
NPLS 相对 NPR 的追加增益、相对川渝自身排放贡献的比例。λ_r 用 9 个 CMIP6 系数加
GISTEMP 回归值；HWF 用近期气候背景下的响应曲线换算。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'; OUT.mkdir(exist_ok=True)
REG = {'sichuan': '四川', 'chongqing': '重庆'}
GISTEMP_LAMBDA = {'sichuan': 2.078933, 'chongqing': 1.296840}
MAIN = dict(NPI='NPI_aeei1.0', NPR='NPR', NPLS='NPLS_gcam_mid')
YEARS_REPORT = [2025, 2030, 2035, 2050, 2060]


def lambdas():
    d = pd.read_csv(IN / 'regional_response_coefficients.csv')
    d = d[d.variable.eq('tasmax') & d.season.eq('annual')]
    out = {}
    for r in REG:
        k = d[d.region.eq(r)][['model', 'window', 'kappa']]
        lab = [f'{m}|{w}' for m, w in zip(k.model, k.window)] + ['GISTEMP']
        out[r] = (np.append(k.kappa.to_numpy(), GISTEMP_LAMBDA[r]), lab)
    return out


def hw_curve(r):
    c = pd.read_csv(HW_OUT / 'response_curve.csv')
    c = c[c.background.eq('recent') & c.region.eq(r) & c.definition.eq('annual90') & c.kind.eq('day') & c.weighting.eq('population2020')].sort_values('offset_tx')
    x, y = c.offset_tx.to_numpy(), c.hwf_mean.to_numpy()
    hi = (y[-1] - y[-2]) / (x[-1] - x[-2]); lo = (y[1] - y[0]) / (x[1] - x[0])
    return lambda w: np.where(w > x[-1], y[-1] + hi * (w - x[-1]), np.where(w < x[0], y[0] + lo * (w - x[0]), np.interp(w, x, y)))


def cmip_w():
    d = pd.read_csv(SSP245DIR / 'ssp245_period_monthly_deltas.csv')
    d = d[d.variable.eq('tasmax')].groupby(['region', 'target_year', 'model']).delta_c.mean()
    return {r: {int(y): d.loc[(r, y)].to_numpy() for y in [2030, 2035]} for r in REG}


CMIP_W = cmip_w()


def q(a, p):
    return float(np.nanquantile(a, p))


def main():
    z = np.load(HERE / 'fair/global_T.npz'); T = z['T']; yrs = z['years']; cases = list(z['cases'])
    base = T[(yrs >= 2015) & (yrs <= 2024), cases.index('NPR'), :].mean(axis=0)       # 每个成员自己的基期
    g = {c: T[:, cases.index(c), :] - base[None, :] for c in cases}                    # 全球相对基期增温 (年份, 成员)
    L = lambdas()
    em = pd.read_csv(HERE / 'scenario_paths_energy_co2.csv')
    rows, sens, order, table, sens_own = [], [], [], [], []
    for r, rn in REG.items():
        lam, lab = L[r]; lam_c = np.median(lam[:-1]); hwf = hw_curve(r)
        for y in YEARS_REPORT:
            yi = int(np.where(yrs == y)[0][0])
            gw = {c: g[c][yi] for c in cases}                                          # (成员,)
            def metrics(npi, npr, npls, lam_arr, zero='zero_SCCQ'):
                la = np.asarray(lam_arr)[:, None]                                        # (λ, 成员)
                W = {k: la * gw[v][None, :] for k, v in dict(NPI=npi, NPR=npr, NPLS=npls).items()}
                Wz = la * gw[zero][None, :]
                A_npr, A_npls = W['NPI'] - W['NPR'], W['NPI'] - W['NPLS']
                with np.errstate(divide='ignore', invalid='ignore'):
                    out = dict(W_NPI=W['NPI'], W_NPR=W['NPR'], W_NPLS=W['NPLS'],
                               A_NPI_NPR=A_npr, A_NPI_NPLS=A_npls, A_NPR_NPLS=A_npls - A_npr,
                               P_NPR=100 * A_npr / W['NPI'], P_NPLS=100 * A_npls / W['NPI'], P_local=100 * (A_npls - A_npr) / W['NPR'],
                               G_local_gain=np.where(np.abs(A_npr) > 1e-12, 100 * (A_npls - A_npr) / A_npr, np.nan),
                               C_NPR=100 * A_npr / (W['NPI'] - Wz), C_NPLS=100 * A_npls / (W['NPI'] - Wz),
                               HWF_NPI=hwf(W['NPI']), HWF_NPR=hwf(W['NPR']), HWF_NPLS=hwf(W['NPLS']))
                    b = hwf(np.zeros(1))[0]
                    out.update(dHWF_NPI_NPR=out['HWF_NPI'] - out['HWF_NPR'], dHWF_NPI_NPLS=out['HWF_NPI'] - out['HWF_NPLS'],
                               dHWF_NPR_NPLS=out['HWF_NPR'] - out['HWF_NPLS'],
                               PH_NPR=100 * (out['HWF_NPI'] - out['HWF_NPR']) / (out['HWF_NPI'] - b),
                               PH_NPLS=100 * (out['HWF_NPI'] - out['HWF_NPLS']) / (out['HWF_NPI'] - b))
                    out['order_ok'] = (W['NPI'] > W['NPR']) & (W['NPR'] > W['NPLS'])
                    if y in CMIP_W[r]:
                        wc = np.asarray(CMIP_W[r][y])[:, None, None]                       # (3 模式, λ, 成员)
                        out['P_NPLS_cmipdenom'] = (100 * A_npls[None] / wc).reshape(-1, A_npls.shape[-1])
                        out['P_NPR_cmipdenom'] = (100 * A_npr[None] / wc).reshape(-1, A_npr.shape[-1])
                return out
            # 主设定：交叉（841 × 10 个 λ）
            tag = {'sichuan': 'SC', 'chongqing': 'CQ'}[r]
            m_main = metrics(MAIN['NPI'], MAIN['NPR'], MAIN['NPLS'], lam)
            m_fair = metrics(MAIN['NPI'], MAIN['NPR'], MAIN['NPLS'], [lam_c])            # 只变 FaIR
            own = metrics(MAIN['NPI'] + f'__{tag}only', 'NPR', MAIN['NPLS'] + f'__{tag}only', lam, zero=f'zero__{tag}only')
            own_f = metrics(MAIN['NPI'] + f'__{tag}only', 'NPR', MAIN['NPLS'] + f'__{tag}only', [lam_c], zero=f'zero__{tag}only')
            m_main.update({'own_' + k: v for k, v in own.items() if k != 'order_ok' and not k.startswith('HWF') and 'cmip' not in k})
            m_fair.update({'own_' + k: v for k, v in own_f.items() if k != 'order_ok' and not k.startswith('HWF') and 'cmip' not in k})
            for k, v in m_main.items():
                if k == 'order_ok':
                    continue
                vf = m_fair.get(k, v)
                med_by_lam = np.nanmedian(v, axis=1) if 'cmip' not in k else np.nanmedian(v.reshape(3, len(lam), -1), axis=(0, 2))
                rows.append(dict(region=rn, year=y, metric=k, median=q(v, .5), p05=q(v, .05), p95=q(v, .95),
                                 fair_only_p05=q(vf, .05), fair_only_median=q(vf, .5), fair_only_p95=q(vf, .95),
                                 lambda_only_min=float(np.nanmin(med_by_lam)), lambda_only_max=float(np.nanmax(med_by_lam)),
                                 lambda_gistemp=float(med_by_lam[-1])))
            order.append(dict(region=rn, year=y, share_NPI_gt_NPR_gt_NPLS=float(m_main['order_ok'].mean()),
                              share_NPI_gt_NPR=float((m_main['A_NPI_NPR'] > 0).mean()), share_NPR_gt_NPLS=float((m_main['A_NPR_NPLS'] > 0).mean())))
            # 情景构造敏感性（FaIR 中位、λ 中位）
            for npi in [c for c in cases if c.startswith('NPI') and '__' not in c and 'lowGDP' not in c]:
                for npls in [c for c in cases if c.startswith('NPLS') and 'lowGDP' not in c and '__' not in c]:
                    mm = metrics(npi, 'NPR', npls, [lam_c])
                    sens.append(dict(region=rn, year=y, NPI=npi, NPLS=npls, **{k: q(mm[k], .5) for k in ['A_NPI_NPR', 'A_NPI_NPLS', 'A_NPR_NPLS', 'P_NPR', 'P_NPLS', 'P_local', 'G_local_gain', 'C_NPR', 'C_NPLS', 'dHWF_NPI_NPLS']}))
            # 每个 NPI 构造的单省情景（FaIR 中位、λ 中位）
            for npi in [c for c in cases if c.startswith('NPI') and '__' not in c and 'lowGDP' not in c]:
                for npls in ['NPLS_gcam_mid', 'NPLS_provtarget']:
                    mm = metrics(f'{npi}__{tag}only', 'NPR', f'{npls}__{tag}only', [lam_c], zero=f'zero__{tag}only')
                    sens_own.append(dict(region=rn, year=y, view=f'{tag} only', NPI=npi, NPLS=npls,
                                         **{k: q(mm[k], .5) for k in ['A_NPI_NPR', 'A_NPI_NPLS', 'A_NPR_NPLS', 'P_NPR', 'P_NPLS', 'G_local_gain', 'C_NPR', 'C_NPLS']}))
            mm = metrics('NPI_aeei1.0_lowGDP', 'NPR_lowGDP', 'NPLS_provtarget_lowGDP', [lam_c])
            sens.append(dict(region=rn, year=y, NPI='NPI_aeei1.0_lowGDP', NPLS='NPLS_provtarget_lowGDP(NPR_lowGDP)',
                             **{k: q(mm[k], .5) for k in ['A_NPI_NPR', 'A_NPI_NPLS', 'A_NPR_NPLS', 'P_NPR', 'P_NPLS', 'P_local', 'G_local_gain', 'C_NPR', 'C_NPLS', 'dHWF_NPI_NPLS']}))
    res = pd.DataFrame(rows); res.to_csv(OUT / 'metrics_FaIR_lambda_sensitivity.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(sens).to_csv(OUT / 'metrics_construction_sensitivity.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(sens_own).to_csv(OUT / 'metrics_construction_sensitivity_by_province.csv', index=False, encoding='utf-8-sig')
    pd.DataFrame(order).to_csv(OUT / 'ordering_robustness.csv', index=False, encoding='utf-8-sig')
    # 排放指标
    er = []
    for r, rn in REG.items():
        e = em[em.region.eq(r)]
        for y in [2025, 2030, 2035, 2060]:
            d = {s: e[e.scenario.eq(s)] for s in MAIN.values()}
            ann = {k: float(d[v].set_index('year').co2_mt[y]) for k, v in MAIN.items()}
            en = {k: float(d[v].set_index('year').energy_10k_tce[y]) for k, v in MAIN.items()}
            cum = {k: float(d[v][d[v].year.between(2006, y)].co2_mt.sum()) for k, v in MAIN.items()}
            er.append(dict(region=rn, year=y, **{f'energy_{k}_10k_tce': en[k] for k in MAIN}, **{f'co2_{k}_mt': ann[k] for k in MAIN},
                           **{f'cum2006_{k}_mt': cum[k] for k in MAIN},
                           avoided_cum_NPR_mt=cum['NPI'] - cum['NPR'], avoided_cum_NPLS_mt=cum['NPI'] - cum['NPLS'],
                           avoided_cum_NPR_pct_of_NPI=100 * (cum['NPI'] - cum['NPR']) / cum['NPI'],
                           avoided_cum_NPLS_pct_of_NPI=100 * (cum['NPI'] - cum['NPLS']) / cum['NPI'],
                           local_additional_mt=cum['NPR'] - cum['NPLS'], local_additional_pct_of_NPR=100 * (cum['NPR'] - cum['NPLS']) / cum['NPR'],
                           local_additional_pct_of_NPR_avoided=100 * (cum['NPR'] - cum['NPLS']) / (cum['NPI'] - cum['NPR'])))
    pd.DataFrame(er).to_csv(OUT / 'emission_metrics.csv', index=False, encoding='utf-8-sig')
    pd.set_option('display.width', 250)
    show = res[res.year.isin([2030, 2035, 2060]) & res.metric.isin(['W_NPI', 'W_NPR', 'W_NPLS', 'A_NPI_NPR', 'A_NPI_NPLS', 'A_NPR_NPLS', 'P_NPR', 'P_NPLS', 'P_local', 'G_local_gain', 'C_NPR', 'C_NPLS', 'dHWF_NPI_NPLS', 'PH_NPLS', 'P_NPLS_cmipdenom', 'own_A_NPI_NPLS', 'own_P_NPLS', 'own_C_NPLS', 'own_G_local_gain', 'own_A_NPR_NPLS'])]
    print(show.round(5).to_string())
    print(pd.DataFrame(order).round(3).to_string())
    print(pd.DataFrame(er).round(1).T.to_string())


if __name__ == '__main__':
    main()
