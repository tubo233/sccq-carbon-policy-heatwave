"""补充 FaIR 实验（841 成员，SSP2-4.5 背景，1750–2100）：全球 SSP 基准、三尺度同力度、CH4 敏感性 V0–V4。

2026 起点的基准把各物种相对 SSP2-4.5 的差值整体后移三年；V0 用 CO2 份额代理，V1 用 EDGAR v8.0
2020 份额，V2 用 SSP1-1.9 力度，V3 为 SSP1-2.6 减半，V4 从 2030 开始。细节见 README / Table S8。
"""
import sys
import time
import gc
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parents[1]
EXP = HERE / 'models'
sys.path.insert(0, str(EXP / 'common')); sys.path.insert(0, str(EXP / 'scales'))
from paths import IN, PROV, save, now  # noqa: E402
import run_fair_scales as rs  # noqa: E402
import fair  # noqa: E402
from fair import FAIR  # noqa: E402
from fair.io import read_properties  # noqa: E402
from fair.interface import fill, initialise  # noqa: E402

YEARS = rs.YEARS
OUT = HERE / 'data' / 'fair_supplement'; OUT.mkdir(parents=True, exist_ok=True)
SH = pd.read_csv(HERE / 'data/methane_shares_EDGARv8_2020.csv', index_col=0).value
CH4_CHINA = {'V0': 0.16, 'EDGAR': float(SH['china_share_of_world'])}
CH4_SCCQ = {'V0': None, 'EDGAR': float(SH['sccq_share_of_china_grid'])}      # None 表示用 CO2 份额（THUCMA 2020）


def build():
    species, properties = read_properties(filename=str(IN / 'species.csv'))
    ix = {s: species.index(s) for s in species}
    a = xr.open_dataarray(PROV / 'original_archive_ssps_harmonized_1750-2499.nc', engine='h5netcdf').load()
    B = {s: a.sel(scenario=s, config='unspecified', specie=species, timepoints=YEARS + .5).values for s in ['ssp245', 'ssp126', 'ssp119']}
    base = B['ssp245']
    assert all(np.allclose(B[s][YEARS <= 2022], base[YEARS <= 2022], rtol=0, atol=1e-9, equal_nan=True) for s in B)
    fy = np.arange(rs.START, 2101); fi = fy - 1750; dsel = fy - 2015
    dy, dp = rs.load_dpec()
    nat_thu, sccq_thu, _ = rs.load_thucma()
    share = {sp: float(sccq_thu.loc[v, '2020'] / nat_thu.loc[v, '2020']) for sp, v in rs.THU_VAR.items()}
    bc_oc_nat = {sp: float(nat_thu.loc[rs.THU_VAR[sp], '2020']) for sp in ['BC', 'OC']}
    pm20 = float(dp['REF']['PM25'][2020 - 2015])
    cases, design = {}, []

    def add(name, arr, **meta):
        cases[name] = arr; design.append(dict(case=name, **meta))

    add('ssp245__ref', base.copy(), group='background')
    add('ssp245__identity', base.copy(), group='identity')
    for s in ['ssp119', 'ssp126']:
        add(f'{s}__start2023', B[s].copy(), group='benchmark', start=2023)
        v = base.copy(); t = np.where(YEARS >= 2026)[0]
        v[t] = base[t] + (B[s][t - 3] - base[t - 3])
        nonco2 = [i for sp, i in ix.items() if not sp.startswith('CO2')]
        w = v[np.ix_(t, nonco2)]
        v[np.ix_(t, nonco2)] = np.where(np.isfinite(w), np.maximum(w, np.minimum(base[np.ix_(t, nonco2)], 0)), w)
        add(f'{s}__start2026', v, group='benchmark', start=2026)

    def national_delta():
        d = {sp: dp['OTNZ'][sp][dsel] - dp['REF'][sp][dsel] for sp in ['CO2 FFI', 'Sulfur', 'NOx', 'VOC', 'NH3']}
        dpm = (dp['OTNZ']['PM25'][dsel] - dp['REF']['PM25'][dsel]) / pm20
        d.update({sp: bc_oc_nat[sp] * dpm for sp in ['BC', 'OC']})
        return d

    def universal_delta():
        d = {sp: base[fi, ix[sp]] * (dp['OTNZ'][sp][dsel] / dp['REF'][sp][dsel] - 1) for sp in ['CO2 FFI', 'Sulfur', 'NOx', 'VOC', 'NH3']}
        r = dp['OTNZ']['PM25'][dsel] / dp['REF']['PM25'][dsel]
        d.update({sp: base[fi, ix[sp]] * (r - 1) for sp in ['BC', 'OC']})
        return d

    nd = national_delta()
    deltas = {'universal': universal_delta(), 'national': nd, 'own_share': {sp: share[sp] * v for sp, v in nd.items()}}

    def ch4_delta(alloc, strength='ssp126', start=2023, half=False):
        lag = start - 2023
        src = B['ssp119' if strength == 'ssp119' else 'ssp126'][:, ix['CH4']]; e245 = base[:, ix['CH4']]
        m = np.ones(len(fy))
        ok = fy >= start
        m[ok] = src[fi[ok] - lag] / e245[fi[ok] - lag]
        if half:
            m = 1 - 0.5 * (1 - m)
        return alloc * e245[fi] * (m - 1)

    def alloc(scale, variant):
        if scale == 'universal':
            return 1.0
        c = CH4_CHINA['V0' if variant == 'V0' else 'EDGAR']
        if scale == 'national':
            return c
        s = CH4_SCCQ['V0' if variant == 'V0' else 'EDGAR']
        return c * (share['CO2 FFI'] if s is None else s)

    variants = {'V0': dict(strength='ssp126', start=2023), 'V1': dict(strength='ssp126', start=2023),
                'V2': dict(strength='ssp119', start=2023), 'V3': dict(strength='ssp126', start=2023, half=True),
                'V4': dict(strength='ssp126', start=2030)}
    ch4_rows = []
    for scale, d in deltas.items():
        for sset in ['co2only', 'all']:
            v = base.copy()
            for sp, dv in d.items():
                if sset == 'co2only' and sp != 'CO2 FFI':
                    continue
                v[fi, ix[sp]] += dv
            add(f'{scale}__OTNZ__{sset}', v, group='scale', scale=scale, species_set=sset)
        for vn, kw in variants.items():
            if scale == 'universal' and vn == 'V0':
                continue                                   # 全球尺度下分配方式无影响：V0 = V1
            a_ = alloc(scale, vn); dch4 = ch4_delta(a_, **kw)
            v = base.copy()
            for sp, dv in d.items():
                v[fi, ix[sp]] += dv
            v[fi, ix['CH4']] += dch4
            add(f'{scale}__OTNZ__allCH4_{vn}', v, group='methane', scale=scale, species_set='allCH4', variant=vn, alloc=a_, **kw)
            ch4_rows += [dict(case=f'{scale}__OTNZ__allCH4_{vn}', year=y, dCH4_Mt=float(dch4[y - rs.START])) for y in [2030, 2035, 2050, 2060, 2100]]
        v = base.copy(); v[fi, ix['CH4']] += ch4_delta(alloc(scale, 'V1'))
        add(f'{scale}__CH4only_V1', v, group='methane_only', scale=scale, variant='V1')
    meta = dict(share_thucma_2020=share, ch4_china_share=CH4_CHINA, ch4_sccq_share_of_china=CH4_SCCQ, variants=variants)
    return species, properties, cases, pd.DataFrame(design), pd.DataFrame(ch4_rows), meta


def run():
    t0 = time.monotonic()
    species, properties, cases, design, ch4_rows, meta = build()
    names = list(cases)
    design.to_csv(OUT / 'case_design.csv', index=False); ch4_rows.to_csv(OUT / 'methane_deltas.csv', index=False)
    params = pd.read_csv(IN / 'parameters.csv', index_col=0)
    original = pd.read_csv(IN / 'original_parameters.csv', index_col=0)
    ids = sorted(params.index)
    so = pd.read_csv(PROV / 'calibration1.4.1_solar_erf_timebounds.csv').set_index('year').erf.loc[1750:2101].to_numpy()
    vo = pd.read_csv(PROV / 'calibration1.4.1_volcanic_erf_timebounds.csv').set_index('timebounds').erf.loc[1750:2101].to_numpy()
    temps = []
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
        t = (raw[:-1] + raw[1:]) / 2
        ref = names.index('ssp245__ref')
        assert float(np.abs(t[:, names.index('ssp245__identity')] - t[:, ref]).max()) == 0.0
        for n in names:
            assert np.allclose(t[YEARS < 2023, names.index(n)], t[YEARS < 2023, ref], rtol=0, atol=1e-9), n
            if n.endswith('start2026'):
                assert np.array_equal(t[YEARS < 2026, names.index(n)], t[YEARS < 2026, ref]), n
        temps.append(t[YEARS >= 1995].astype(np.float32))
        print(f'{start + len(members)}/{len(ids)} members, {len(names)} cases, {time.monotonic() - t0:.0f}s', flush=True)
        del f, raw, t; gc.collect()
    T = np.concatenate(temps, axis=2)
    np.savez_compressed(OUT / 'global_temperature.npz', years=YEARS[YEARS >= 1995], cases=np.array(names), members=np.array(ids), T=T)
    save(OUT / 'protocol.json', dict(created_utc=now(), fair_version=fair.__version__, members=len(ids), cases=names, meta=meta,
                                     seconds=time.monotonic() - t0))
    print('done', time.monotonic() - t0)


if __name__ == '__main__':
    run()
