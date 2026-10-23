"""补充 FaIR 实验（08）的区域 Tmax 变化、交叉年份和信噪比显现。

区域变化 = 全球 ΔT × 三个 CMIP6 模式的 κ，共 841 × 3 = 2,523 组；交叉年份取中位变化
此后保持为负的年份；显现阈值用 2σ/√20（20 年均值）或 σ（单年），'likely' 取 83 分位。
σ 是观测人口加权年均 Tmax 的去趋势年际标准差（1990–2024）。
"""
import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
EXP = HERE / 'models'
sys.path.insert(0, str(EXP / 'common'))
from paths import kappa_table  # noqa: E402

RUN = HERE / 'data' / 'fair_supplement'
MODELS = ['GFDL-ESM4', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']
REGIONS = ['sichuan', 'chongqing']
YRS = np.arange(2023, 2101)


def main():
    z = np.load(RUN / 'global_temperature.npz'); T, years, cases = z['T'], z['years'], list(z['cases'])
    k = kappa_table(); K = {r: np.array([k.loc[(r, 'tasmax'), m] for m in MODELS]) for r in REGIONS}
    noise = json.load(open(EXP / 'scales/analysis_inputs.json', encoding='utf-8'))['noise']
    ref = T[:, cases.index('ssp245__ref'), :]
    sel = np.isin(years, YRS)
    ts, summ = [], []
    for c in cases:
        dg = (T[:, cases.index(c), :] - ref)[sel]                               # (年份, 成员)
        for r in REGIONS:
            dr = (dg[:, :, None] * K[r][None, None, :]).reshape(len(YRS), -1)    # (年份, 成员×模式)
            q = {p: np.quantile(dr, p / 100, axis=1) for p in [5, 17, 50, 83, 95]}
            share_cool = (dr < 0).mean(axis=1)
            for i, y in enumerate(YRS):
                ts.append(dict(case=c, region=r, year=int(y), share_cooling=share_cool[i], **{f'dTmax_p{p:02d}': q[p][i] for p in q}))
            s = noise[r]['t']; thr20 = 2 * s / np.sqrt(20)

            def persistent(x):
                for i in range(len(x)):
                    if x[i:].all():
                        return int(YRS[i])
                return None
            at = lambda y: float(q[50][YRS == y][0])
            summ.append(dict(case=c, region=r, crossover=persistent(q[50] < 0),
                             toe20_median=persistent(-q[50] >= thr20), toe20_likely=persistent(-q[83] >= thr20),
                             toe1_median=persistent(-q[50] >= s), thr20=thr20, thr1=s,
                             d2035=at(2035), d2035_p05=float(q[5][YRS == 2035][0]), d2035_p95=float(q[95][YRS == 2035][0]),
                             share_warming_2035=float(1 - share_cool[YRS == 2035][0]),
                             d2050=at(2050), d2060=at(2060), d2100=at(2100)))
    ts = pd.DataFrame(ts); summ = pd.DataFrame(summ)
    ts.to_csv(RUN / 'regional_timeseries_2023_2100.csv', index=False)
    summ.to_csv(RUN / 'crossover_emergence_summary.csv', index=False)

    # 自捕获：同力度下川渝份额 / 全球（仅 CO2，全球中位）
    g = lambda c, y: np.median(T[years == y, cases.index(c), :][0] - ref[years == y][0])
    sc = pd.DataFrame([dict(year=y, own_over_world_pct=100 * g('own_share__OTNZ__co2only', y) / g('universal__OTNZ__co2only', y),
                            china_over_world_pct=100 * g('national__OTNZ__co2only', y) / g('universal__OTNZ__co2only', y))
                       for y in [2035, 2050, 2060, 2100]])
    sc.to_csv(RUN / 'self_capture_co2only.csv', index=False)

    pd.set_option('display.width', 250)
    cols = ['case', 'crossover', 'toe20_median', 'toe20_likely', 'toe1_median', 'd2035', 'd2035_p05', 'd2035_p95', 'share_warming_2035', 'd2050', 'd2060', 'd2100']
    for r in REGIONS:
        print(r); print(summ[summ.region == r][cols].round(5).to_string())
    print(sc.round(3).to_string())


if __name__ == '__main__':
    main()
