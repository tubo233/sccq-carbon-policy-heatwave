"""三个空间尺度上的 DPEC 政策力度 FaIR 实验（841 成员，SSP2-4.5 背景），外加历史排放上界反事实。

尺度分 own_share、national、universal；扰动从 2023 年开始，2060 年后差值/比例保持。
DPEC 没有 BC/OC，就按 PM2.5 的相对变化近似；CH4 变体和历史反事实在 build_cases 里。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *
import time, gc
import numpy as np
import pandas as pd
import xarray as xr
import fair
from fair import FAIR
from fair.io import read_properties
from fair.interface import fill, initialise

YEARS = np.arange(1750, 2101)
START = 2023
OUT = FAIR_OUT
DPEC_SCEN = {'REF': 'reference', 'CA': 'clean air', 'OTNZ': 'on-time peak-net zero-clean air', 'EPNZ': 'early peak-net zero-clean air'}
DPEC_FILES = {'CO2 FFI': 'CO2', 'Sulfur': 'SO2', 'NOx': 'NOx', 'VOC': 'VOC', 'NH3': 'NH3', 'PM25': 'PM25'}
THU_VAR = {'CO2 FFI': 'Emissions|CO2', 'Sulfur': 'Emissions|Sulfur', 'NOx': 'Emissions|NOx', 'VOC': 'Emissions|VOC',
           'NH3': 'Emissions|NH3', 'BC': 'Emissions|BC', 'OC': 'Emissions|OC'}
POLLUTANTS = ['Sulfur', 'NOx', 'VOC', 'NH3', 'BC', 'OC']


def annual(series_years, values, years):
    return np.interp(years, series_years, values)


def load_dpec():
    """返回 dict[scen][species] -> 2015..2100 年数组（单位与 FaIR 一致：CO2 Gt，其余 Mt）。2060 后保持 2060 值。"""
    yrs = np.arange(2015, 2101)
    out = {k: {} for k in DPEC_SCEN}
    for sp, f in DPEC_FILES.items():
        d = pd.read_csv(EMIS / f'DPEC7383239_{f}.csv')
        for k, col in DPEC_SCEN.items():
            out[k][sp] = annual(d.year.to_numpy(), d[col].to_numpy(float), yrs)
    return yrs, out


def load_thucma():
    d = pd.read_csv(EMIS / 'THUCMA_provincial_totals_from_sector_table.csv')
    cols = [str(y) for y in range(2020, 2061, 5)]
    nat = d.groupby('Variable')[cols].sum()
    sccq = d[d.Region.isin(['Sichuan', 'Chongqing'])].groupby('Variable')[cols].sum()
    return nat, sccq, np.arange(2020, 2061, 5)


def ceads():
    d = pd.read_csv(CEADS, encoding='utf-8-sig')
    s = d.groupby('年份')['CO2总排放_百万吨'].sum(min_count=2)
    return s


def build_cases():
    species, properties = read_properties(filename=str(IN / 'species.csv'))
    ix = {s: species.index(s) for s in species}
    a = xr.open_dataarray(PROV / 'original_archive_ssps_harmonized_1750-2499.nc', engine='h5netcdf').load()
    bases = {s: a.sel(scenario=s, config='unspecified', specie=species, timepoints=YEARS + .5).values for s in ['ssp245', 'ssp126', 'ssp370', 'ssp119']}
    base = bases['ssp245']
    fy = np.arange(START, 2101)          # 扰动年份
    fi = fy - 1750                       # 行索引
    dy, dp = load_dpec()
    dsel = fy - 2015
    nat_thu, sccq_thu, thu_years = load_thucma()
    share = {sp: float(sccq_thu.loc[v, '2020'] / nat_thu.loc[v, '2020']) for sp, v in THU_VAR.items()}
    bc_oc_nat2020 = {sp: float(nat_thu.loc[THU_VAR[sp], '2020']) for sp in ['BC', 'OC']}
    bc_oc_sccq2020 = {sp: float(sccq_thu.loc[THU_VAR[sp], '2020']) for sp in ['BC', 'OC']}
    pm20 = float(dp['REF']['PM25'][2020 - 2015])

    def national_delta(s):
        """中国全国各物种差值（s - REF），FaIR 单位，长度 len(fy)。"""
        d = {}
        for sp in ['CO2 FFI', 'Sulfur', 'NOx', 'VOC', 'NH3']:
            d[sp] = dp[s][sp][dsel] - dp['REF'][sp][dsel]
        dpm = (dp[s]['PM25'][dsel] - dp['REF']['PM25'][dsel]) / pm20
        for sp in ['BC', 'OC']:
            d[sp] = bc_oc_nat2020[sp] * dpm
        return d

    def universal_delta(s):
        d = {}
        for sp in ['CO2 FFI', 'Sulfur', 'NOx', 'VOC', 'NH3']:
            ratio = dp[s][sp][dsel] / dp['REF'][sp][dsel]
            d[sp] = base[fi, ix[sp]] * (ratio - 1)
        ratio = dp[s]['PM25'][dsel] / dp['REF']['PM25'][dsel]
        for sp in ['BC', 'OC']:
            d[sp] = base[fi, ix[sp]] * (ratio - 1)
        return d

    cases, design, deltas = {}, [], []

    def add(name, arr, **meta):
        cases[name] = arr
        design.append(dict(case=name, **meta))

    for bg, b in bases.items():
        add(f'{bg}__ref', b.copy(), kind='background', scale='none', stringency='none', species_set='none')
    add('ssp245__identity', base.copy(), kind='identity', scale='none', stringency='none', species_set='none')

    # CH4 变体：DPEC 不含 CH4。净零路径通常同时削减 CH4，这里令全球 CH4 取 SSP1-2.6 路径（universal），
    # 全国按中国约占全球人为 CH4 的 16%（EDGAR 量级，近似）承担，川渝再按 CO2 份额承担（近似）。
    ch4_world = bases['ssp126'][fi, ix['CH4']] - base[fi, ix['CH4']]
    CH4_CHINA_SHARE = 0.16
    for s in ['CA', 'OTNZ', 'EPNZ']:
        nd = national_delta(s)
        ud = universal_delta(s)
        variants = {'national': nd, 'universal': ud, 'own_share': {sp: share[sp] * v for sp, v in nd.items()}}
        ch4 = {'universal': ch4_world, 'national': CH4_CHINA_SHARE * ch4_world, 'own_share': share['CO2 FFI'] * CH4_CHINA_SHARE * ch4_world}
        for scale, d in variants.items():
            for sset in (['all', 'co2only', 'allCH4'] if s != 'CA' else ['all']):
                v = base.copy()
                if sset == 'allCH4':
                    v[fi, ix['CH4']] += ch4[scale]
                for sp, dv in d.items():
                    if sset == 'co2only' and sp != 'CO2 FFI':
                        continue
                    v[fi, ix[sp]] += dv
                name = f'ssp245__{scale}__{s}__{sset}'
                add(name, v, kind='policy', scale=scale, stringency=s, species_set=sset)
                for sp, dv in d.items():
                    if sset == 'co2only' and sp != 'CO2 FFI':
                        continue
                    for y in [2030, 2035, 2050, 2060]:
                        deltas.append(dict(case=name, species=sp, year=y, delta=float(dv[y - START])))

    # 历史上界：正差值 = 反事实比实际多排（即政策“避免”的排放），取负号加到历史上得到“实际世界”相对“反事实世界”的差
    ce = ceads()
    hy = np.arange(2006, 2023)
    actual = ce.loc[hy].to_numpy(float)
    y0, y1 = 2005, 2012
    g = np.polyfit(np.arange(y0, y1 + 1), np.log(ce.loc[y0:y1].to_numpy(float)), 1)
    cf_exp = np.exp(np.polyval(g, hy))
    l = np.polyfit(np.arange(y0, y1 + 1), ce.loc[y0:y1].to_numpy(float), 1)
    cf_lin = np.polyval(l, hy)
    # 全国近似序列（仅用于 20% 上界；锚点：2006 年 6.9、2022 年 11.4 Gt，2015/2020 年取 DPEC）
    china = np.interp(hy, [2006, 2015, 2020, 2022], [6.9, 9.52, 10.38, 11.4])
    hist_meta = {}
    for name, extra in [('hist__SCCQ_exp2005_2012', (cf_exp - actual) / 1000),
                        ('hist__SCCQ_lin2005_2012', (cf_lin - actual) / 1000),
                        ('hist__China_plus20pct', 0.2 * china)]:
        v = base.copy()
        v[hy - 1750, ix['CO2 FFI']] += extra   # 反事实世界 = 实际 + 避免的排放
        add(name, v, kind='historical_counterfactual', scale=name.split('__')[1], stringency='bound', species_set='co2only')
        hist_meta[name] = dict(cumulative_extra_gtco2_2006_2022=float(extra.sum()), extra_2022_gtco2=float(extra[-1]))

    meta = dict(share_sccq_2020_thucma=share, bc_oc_national_2020=bc_oc_nat2020, bc_oc_sccq_2020=bc_oc_sccq2020,
                historical=hist_meta, start_year=START,
                ceads_sccq_2020_mt=float(ce.loc[2020]), ceads_sccq_2022_mt=float(ce.loc[2022]))
    return species, properties, cases, pd.DataFrame(design), pd.DataFrame(deltas), meta


def run(limit=None):
    t0 = time.monotonic()
    species, properties, cases, design, deltas, meta = build_cases()
    names = list(cases)
    params = pd.read_csv(IN / 'parameters.csv', index_col=0)
    original = pd.read_csv(IN / 'original_parameters.csv', index_col=0)
    ids = sorted(params.index)
    ids = ids[:limit] if limit else ids
    out = OUT / ('smoke' if limit else 'run')
    out.mkdir(parents=True, exist_ok=True)
    design.to_csv(out / 'case_design.csv', index=False)
    deltas.to_csv(out / 'emission_deltas_by_case.csv', index=False)
    save(out / 'protocol.json', dict(created_utc=now(), fair_version=fair.__version__, members=len(ids), cases=names, meta=meta,
                                    parameter_release='fair-calibrate1.4.1; no recalibration',
                                    note='扰动自2023年起；2060年后差值/比值保持2060；BC/OC按一次PM2.5相对变化近似；CH4/N2O/CO等未改变'))
    so = pd.read_csv(PROV / 'calibration1.4.1_solar_erf_timebounds.csv').set_index('year').erf.loc[1750:2101].to_numpy()
    vo = pd.read_csv(PROV / 'calibration1.4.1_volcanic_erf_timebounds.csv').set_index('timebounds').erf.loc[1750:2101].to_numpy()
    assert len(so) == len(vo) == len(YEARS) + 1
    temps, aer = [], []
    for start in range(0, len(ids), 32):
        members = ids[start:start + 32]
        f = FAIR(ch4_method='Thornhill2021')
        f.define_time(1750, 2101, 1); f.define_scenarios(names); f.define_configs(members)
        f.define_species(species, properties); f.allocate()
        f.fill_species_configs(str(IN / 'species.csv')); f.override_defaults(str(IN / 'parameters.csv'))
        fill(f.climate_configs['stochastic_run'], False)
        for n, arr in cases.items():
            f.emissions.loc[{'scenario': n}] = arr[:, None, :]
        for s, v in [('Solar', so), ('Volcanic', vo)]:
            f.forcing.loc[{'specie': s}] = v[:, None, None] * params.loc[members, f'forcing_scale[{s}]'].to_numpy()[None, None, :]
        f.forcing.loc[{'specie': 'Solar'}] += np.minimum((f.timebounds - 1750) / 270, 1)[:, None, None] * original.loc[members, 'fscale_solar_trend'].to_numpy()[None, None, :]
        initialise(f.concentration, f.species_configs['baseline_concentration'])
        for v in [f.forcing, f.temperature, f.cumulative_emissions, f.airborne_emissions, f.ocean_heat_content_change]:
            initialise(v, 0.)
        f.run(progress=False)
        raw = f.temperature.sel(layer=0).values
        assert np.isfinite(raw).all()
        t = (raw[:-1] + raw[1:]) / 2          # 年均（相邻 timebounds 平均）
        ref = names.index('ssp245__ref')
        assert float(np.abs(t[:, names.index('ssp245__identity')] - t[:, ref]).max()) == 0.0
        pre = YEARS < START
        for n in names:
            if n.startswith('ssp245__') and not n.startswith('hist'):
                assert np.array_equal(t[pre, names.index(n)], t[pre, ref]), n
        er = f.forcing.sel(specie='Aerosol-radiation interactions').values + f.forcing.sel(specie='Aerosol-cloud interactions').values
        keep = YEARS >= 1995
        temps.append(t[keep].astype(np.float32))
        aer.append(((er[:-1] + er[1:]) / 2)[keep].astype(np.float32))
        print(f'{start + len(members)}/{len(ids)} members, {len(names)} cases, {time.monotonic() - t0:.0f}s', flush=True)
        del f, raw, t; gc.collect()
    T = np.concatenate(temps, axis=2)   # (年份, 情景, 成员)
    A = np.concatenate(aer, axis=2)
    np.savez_compressed(out / 'global_temperature.npz', years=YEARS[YEARS >= 1995], cases=np.array(names), members=np.array(ids), T=T, aerosol_erf=A)
    print('done', time.monotonic() - t0)


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(); p.add_argument('--limit', type=int); a = p.parse_args()
    run(a.limit)
