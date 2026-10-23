"""观测噪声尺度：人口加权年均 Tmax 与人口加权（先格点识别）年热浪天数的年际标准差（1990-2024，去线性趋势）。"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *
import importlib.util
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location('legacy_weather', LEGACY / 'scripts/common.py')
old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)


def main():
    dates, meta, tx, tn = old.load_cache()
    txa = np.asarray(tx)
    train = np.asarray((dates.year >= 1990) & (dates.year <= 2005))
    q90 = np.quantile(txa[train], .9, axis=0)
    hot = txa > q90[None]
    tri = hot[:-2] & hot[1:-1] & hot[2:]
    marked = np.zeros_like(hot)
    for i in range(3):
        marked[i:i + len(tri)] |= tri
    years = np.unique(dates.year)
    rows = []
    for r in old.REGIONS:
        m = (meta.region_id == r).to_numpy()
        for wname, col in [('population2020', 'population2020'), ('area', 'area_weight')]:
            w = meta[col].to_numpy()[m]; w = w / w.sum()
            daily_t = txa[:, m] @ w
            daily_h = marked[:, m].astype(float) @ w
            for y in years:
                s = np.asarray(dates.year == y)
                rows.append(dict(region=r, weighting=wname, year=int(y), tmax_mean=float(daily_t[s].mean()), hwf=float(daily_h[s].sum())))
    df = pd.DataFrame(rows)
    df.to_csv(HW_OUT / 'observed_annual_tmax_hwf_1990_2024.csv', index=False)
    out = []
    for (r, wname), g in df.groupby(['region', 'weighting']):
        for var in ['tmax_mean', 'hwf']:
            x = g.year.to_numpy(float); y = g[var].to_numpy()
            p = np.polyfit(x, y, 1); res = y - np.polyval(p, x)
            ac1 = float(np.corrcoef(res[:-1], res[1:])[0, 1])
            out.append(dict(region=r, weighting=wname, variable=var, trend_per_decade=p[0] * 10, detrended_sd=float(res.std(ddof=2)),
                            lag1_autocorr=ac1, mean_1990_2005=float(y[x <= 2005].mean()), mean_2015_2024=float(y[x >= 2015].mean())))
    pd.DataFrame(out).to_csv(HW_OUT / 'observed_noise_scales.csv', index=False)
    print(pd.DataFrame(out).round(3).to_string())


if __name__ == '__main__':
    main()
