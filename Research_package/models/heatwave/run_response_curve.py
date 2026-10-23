"""人口加权热浪天数对区域 Tmax 均匀偏移的响应曲线。

背景：recent = 2015-2024 观测月气候；GFDL/MPI/UKESM = 同一观测基期叠加各模式 SSP2-4.5 的 2035 年（2026-2045）月变化。
偏移：对 Tmax 加 δ，对 Tmin 加 δ×(κ_tasmin/κ_tasmax)（三模式中位）。
天气：整年残差捐赠方案（128 成员，SEED=20260924+2035），所有情景共用同一残差。
热浪：1990-2005 逐格点固定阈值（annual90/annual95/warm90），严格超阈且连续 ≥3 天；先格点识别，再人口/面积加权。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *

os.environ.setdefault('NUMBA_NUM_THREADS', '8')
import importlib.util, time
import numpy as np
import pandas as pd
from numba import njit, prange
from threadpoolctl import threadpool_limits

spec = importlib.util.spec_from_file_location('legacy_weather', LEGACY / 'scripts/common.py')
old = importlib.util.module_from_spec(spec); spec.loader.exec_module(old)
OUT = HW_OUT
YEAR = 2035
MEMBERS = 128
SEED = 20260924
OFFSETS = [-0.2, -0.1, -0.05, -0.02, 0.0, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0]
BACKGROUNDS = ['recent', 'GFDL-ESM4', 'MPI-ESM1-2-LR', 'UKESM1-0-LL']


@njit(parallel=True, cache=True)
def count_cases(mid, half, delta_mid, ratio, offset, q, month, warm, inyear):
    nt, nm, ng = mid.shape; nc = delta_mid.shape[0]
    out = np.zeros((nc, nm, ng, 3, 3, 2), np.uint16)
    for g in prange(ng):
        for m in range(nm):
            runs = np.zeros((nc, 3, 3), np.int32)
            for t in range(nt):
                for c in range(nc):
                    mm = mid[t, m, g] + delta_mid[c, month[t], g]
                    hh = half[t, m, g] * ratio[c, month[t], g]
                    sx = mm + hh + offset[c, g, 0]; sn = mm - hh + offset[c, g, 1]
                    for d in range(3):
                        hd = sx > q[d, g, 0]; hn = sn > q[d, g, 1]
                        if d == 2 and not warm[t]:
                            hd = False; hn = False
                        for k in range(3):
                            hot = hd if k == 0 else (hn if k == 1 else hd and hn)
                            if hot:
                                runs[c, d, k] += 1; r = runs[c, d, k]
                                if r == 3:
                                    if inyear[t - 2]: out[c, m, g, d, k, 1] += 1
                                    for back in range(3):
                                        if inyear[t - back]: out[c, m, g, d, k, 0] += 1
                                elif r > 3 and inyear[t]:
                                    out[c, m, g, d, k, 0] += 1
                            else:
                                runs[c, d, k] = 0
    return out


def setup():
    k = kappa_table().median(axis=1)
    tn_ratio = {r: float(k.loc[(r, 'tasmin')] / k.loc[(r, 'tasmax')]) for r in old.REGIONS}
    cases = [dict(name=f'{bg}__{d:+.2f}', background=bg, offset_tx=d) for bg in BACKGROUNDS for d in OFFSETS]
    return cases, tn_ratio


def monthly_targets(tx, tn, dates, cases, grid_regions):
    periods = pd.read_csv(SSP245DIR / 'ssp245_period_monthly_deltas.csv')
    q = periods.set_index(['model', 'region', 'variable', 'target_year', 'month']).delta_c
    target_months = [12] + list(range(1, 13)) + [1]
    ng = tx.shape[1]; out = np.empty((len(cases), 14, ng, 2))
    observed = {m: np.stack([np.stack([tx[(dates.year == y) & (dates.month == m)].mean(axis=0),
                                       tn[(dates.year == y) & (dates.month == m)].mean(axis=0)], axis=-1)
                             for y in range(2015, 2025)]).mean(axis=0) for m in range(1, 13)}
    for ci, c in enumerate(cases):
        for j, m in enumerate(target_months):
            out[ci, j] = observed[m]
            if c['background'] == 'recent':
                continue
            for r in old.REGIONS:
                mask = np.asarray(grid_regions) == r
                for vi, var in enumerate(['tasmax', 'tasmin']):
                    out[ci, j, mask, vi] += q.loc[(c['background'], r, var, YEAR, m)]
    return out


def run(chunk=64):
    t0 = time.monotonic(); OUT.mkdir(exist_ok=True); chunks = OUT / 'chunks'; chunks.mkdir(exist_ok=True)
    dates, meta, tx, tn = old.load_cache(); cases, tn_ratio = setup()
    pd.DataFrame(cases).to_csv(OUT / 'case_design.csv', index=False)
    save(OUT / 'protocol.json', dict(created_utc=now(), year=YEAR, members=MEMBERS, seed=SEED + YEAR, offsets=OFFSETS,
                                    backgrounds=BACKGROUNDS, tn_over_tx_kappa_ratio=tn_ratio))
    x = old.design(dates); pinv = np.linalg.pinv(x); base = (dates.year >= 1995) & (dates.year <= 2014)
    cp = np.linalg.pinv(old.design(dates[base], False))
    target = pd.date_range(f'{YEAR - 1}-12-30', f'{YEAR + 1}-01-02')
    a, b, w = old.donor_plan(dates, target, 'whole_year', MEMBERS, SEED + YEAR)
    month = np.where(target.year < YEAR, 0, np.where(target.year > YEAR, 13, target.month)).astype(np.int64)
    inyear = np.asarray(target.year == YEAR); warm = np.asarray((target.month >= 5) & (target.month <= 9))
    nw = 2
    values = np.zeros((len(cases), MEMBERS, 2, 3, 3, nw, 2))
    for lo in range(0, len(meta), chunk):
        hi = min(lo + chunk, len(meta)); fp = chunks / f'grids_{lo:04}_{hi:04}.npz'
        gm = meta.iloc[lo:hi]
        if fp.exists():
            raw = np.load(fp)['raw']
        else:
            txg = np.array(tx[:, lo:hi]); tng = np.array(tn[:, lo:hi]); rep = old.represent(txg, tng); flat = rep.reshape(len(dates), -1)
            residual = (flat - x @ (pinv @ flat)).reshape(rep.shape)
            q = np.stack([np.stack(old.thresholds(txg, tng, dates, d), axis=-1) for d in old.DEFINITIONS])
            fc = old.design(target, False)
            mean = (fc @ (cp @ flat[base])).reshape(len(target), hi - lo, 2)
            daily = mean[:, None] + residual[a] + w[:, :, None, None] * (residual[b] - residual[a])
            mid = np.ascontiguousarray(daily[..., 0]); half = np.ascontiguousarray(np.exp(daily[..., 1]) * .5)
            zm = np.stack([mid[month == j].mean(axis=(0, 1)) for j in range(14)])
            zd = np.stack([(2 * half[month == j]).mean(axis=(0, 1)) for j in range(14)])
            targets = monthly_targets(txg, tng, dates, cases, gm.region_id.to_numpy())
            td = targets[..., 0] - targets[..., 1]
            assert np.all(td > 0)
            dm = targets.mean(axis=-1) - zm[None]; ratio = td / zd[None]
            offsets = np.zeros((len(cases), hi - lo, 2))
            for ci, c in enumerate(cases):
                for r in old.REGIONS:
                    msk = gm.region_id.eq(r).to_numpy()
                    offsets[ci, msk, 0] = c['offset_tx']; offsets[ci, msk, 1] = c['offset_tx'] * tn_ratio[r]
            raw = count_cases(mid, half, np.ascontiguousarray(dm), np.ascontiguousarray(ratio), offsets, q, month, warm, inyear)
            np.savez_compressed(fp, raw=raw)
        for ri, r in enumerate(old.REGIONS):
            for wi, col in enumerate(['population2020', 'area_weight']):
                ww = np.where(gm.region_id.eq(r), gm[col] / meta.loc[meta.region_id.eq(r), col].sum(), 0.)
                values[:, :, ri, :, :, wi, :] += np.einsum('cmgdke,g->cmdke', raw, ww, optimize=True)
        print('grids', hi, f'{time.monotonic() - t0:.0f}s', flush=True)
    np.savez_compressed(OUT / 'response_members.npz', values=values, cases=[c['name'] for c in cases], regions=old.REGIONS,
                        definitions=old.DEFINITIONS, kinds=old.KINDS, weightings=['population2020', 'area'], metrics=['HWF', 'HWN'])
    rows = []
    for ci, c in enumerate(cases):
        for ri, r in enumerate(old.REGIONS):
            for di, d in enumerate(old.DEFINITIONS):
                for ki, k in enumerate(old.KINDS):
                    for wi, wname in enumerate(['population2020', 'area']):
                        v = values[ci, :, ri, di, ki, wi, 0]
                        base_ci = [i for i, cc in enumerate(cases) if cc['background'] == c['background'] and cc['offset_tx'] == 0.0][0]
                        dv = v - values[base_ci, :, ri, di, ki, wi, 0]
                        rows.append(dict(case=c['name'], background=c['background'], offset_tx=c['offset_tx'], region=r, definition=d, kind=k,
                                         weighting=wname, hwf_mean=v.mean(), hwf_sd_members=v.std(ddof=1), hwf_q025=np.quantile(v, .025),
                                         hwf_q975=np.quantile(v, .975), delta_vs_zero_mean=dv.mean(), delta_mcse=dv.std(ddof=1) / np.sqrt(len(dv))))
    pd.DataFrame(rows).to_csv(OUT / 'response_curve.csv', index=False)


if __name__ == '__main__':
    with threadpool_limits(limits=2):
        run()
