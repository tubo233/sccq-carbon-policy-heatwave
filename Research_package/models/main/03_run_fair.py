"""政策情景三情景的 FaIR 实验（841 组校准参数，SSP2-4.5 全球背景，1750-2060）。

参照：NPR = 实际历史 + 国家政策路径，令其等于 SSP2-4.5 背景本身；
其他情景只改变川渝 CO2 FFI：ΔE(t) = [E_s,SC + E_s,CQ − E_NPR,SC − E_NPR,CQ](t)，t ≥ 2006。
另有“川渝零排放”情景（t ≥ 2006 去掉 NPR 下川渝全部 CO2），用于计算川渝自身排放贡献的增温。
"""
import sys, time, gc
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *
import numpy as np
import pandas as pd
import xarray as xr
import fair
from fair import FAIR
from fair.io import read_properties
from fair.interface import fill, initialise

HERE = Path(__file__).resolve().parent
YEARS = np.arange(1750, 2061)
START = 2006


def cases():
    species, properties = read_properties(filename=str(IN / 'species.csv'))
    j = species.index('CO2 FFI')
    a = xr.open_dataarray(PROV / 'original_archive_ssps_harmonized_1750-2499.nc', engine='h5netcdf').load()
    base = a.sel(scenario='ssp245', config='unspecified', specie=species, timepoints=YEARS + .5).values
    main = pd.read_csv(HERE / 'scenario_paths_energy_co2.csv')
    low = pd.read_csv(HERE / 'scenario_paths_energy_co2_lowGDP.csv')
    tot = lambda d, s: d[d.scenario.eq(s)].groupby('year').co2_mt.sum()
    npr = tot(main, 'NPR')
    fy = np.arange(START, 2061); fi = fy - 1750
    out = {'NPR': base.copy(), 'identity': base.copy()}
    design = [dict(case='NPR', source='main'), dict(case='identity', source='main')]
    for s in sorted(main.scenario.unique()):
        if s == 'NPR':
            continue
        v = base.copy(); v[fi, j] += ((tot(main, s) - npr).loc[fy] / 1000).to_numpy()
        out[s] = v; design.append(dict(case=s, source='main'))
    for s in ['NPR', 'NPI_aeei1.0', 'NPLS_provtarget']:
        v = base.copy(); v[fi, j] += ((tot(low, s) - npr).loc[fy] / 1000).to_numpy()
        out[s + '_lowGDP'] = v; design.append(dict(case=s + '_lowGDP', source='lowGDP'))
    v = base.copy(); v[fi, j] -= (npr.loc[fy] / 1000).to_numpy()
    out['zero_SCCQ'] = v; design.append(dict(case='zero_SCCQ', source='main'))
    # 分省：只有一个省改变情景，另一省保持 NPR
    reg = lambda d, s, r: d[d.scenario.eq(s) & d.region.eq(r)].set_index('year').co2_mt
    # 分省：每个 NPI 构造和两个 NPLS 变体各跑一次单省情景
    single = [s for s in sorted(main.scenario.unique()) if s.startswith('NPI')] + ['NPLS_gcam_mid', 'NPLS_provtarget']
    for r, tag in [('sichuan', 'SC'), ('chongqing', 'CQ')]:
        nr = reg(main, 'NPR', r)
        for s in single:
            v = base.copy(); v[fi, j] += ((reg(main, s, r) - nr).loc[fy] / 1000).to_numpy()
            out[f'{s}__{tag}only'] = v; design.append(dict(case=f'{s}__{tag}only', source='main'))
        v = base.copy(); v[fi, j] -= (nr.loc[fy] / 1000).to_numpy()
        out[f'zero__{tag}only'] = v; design.append(dict(case=f'zero__{tag}only', source='main'))
    return species, properties, out, pd.DataFrame(design)


def run():
    t0 = time.monotonic()
    species, properties, cs, design = cases(); names = list(cs)
    params = pd.read_csv(IN / 'parameters.csv', index_col=0); original = pd.read_csv(IN / 'original_parameters.csv', index_col=0)
    ids = sorted(params.index)
    out = HERE / 'fair'; out.mkdir(exist_ok=True)
    design.to_csv(out / 'case_design.csv', index=False)
    save(out / 'protocol.json', dict(created_utc=now(), fair=fair.__version__, members=len(ids), cases=names, start=START,
                                    reference='NPR == SSP2-4.5 background'))
    so = pd.read_csv(PROV / 'calibration1.4.1_solar_erf_timebounds.csv').set_index('year').erf.loc[1750:2061].to_numpy()
    vo = pd.read_csv(PROV / 'calibration1.4.1_volcanic_erf_timebounds.csv').set_index('timebounds').erf.loc[1750:2061].to_numpy()
    temps = []
    for st in range(0, len(ids), 32):
        m = ids[st:st + 32]
        f = FAIR(ch4_method='Thornhill2021'); f.define_time(1750, 2061, 1); f.define_scenarios(names); f.define_configs(m)
        f.define_species(species, properties); f.allocate(); f.fill_species_configs(str(IN / 'species.csv')); f.override_defaults(str(IN / 'parameters.csv'))
        fill(f.climate_configs['stochastic_run'], False)
        for n, arr in cs.items():
            f.emissions.loc[{'scenario': n}] = arr[:, None, :]
        for s, v in [('Solar', so), ('Volcanic', vo)]:
            f.forcing.loc[{'specie': s}] = v[:, None, None] * params.loc[m, f'forcing_scale[{s}]'].to_numpy()[None, None, :]
        f.forcing.loc[{'specie': 'Solar'}] += np.minimum((f.timebounds - 1750) / 270, 1)[:, None, None] * original.loc[m, 'fscale_solar_trend'].to_numpy()[None, None, :]
        initialise(f.concentration, f.species_configs['baseline_concentration'])
        for v in [f.forcing, f.temperature, f.cumulative_emissions, f.airborne_emissions, f.ocean_heat_content_change]:
            initialise(v, 0.)
        f.run(progress=False)
        raw = f.temperature.sel(layer=0).values; assert np.isfinite(raw).all()
        t = (raw[:-1] + raw[1:]) / 2
        assert float(np.abs(t[:, names.index('identity')] - t[:, names.index('NPR')]).max()) == 0
        pre = YEARS < START
        for i in range(len(names)):
            assert np.array_equal(t[pre, i], t[pre, 0])
        temps.append(t[YEARS >= 1990].astype(np.float64))
        print(st + len(m), len(ids), f'{time.monotonic() - t0:.0f}s', flush=True)
        del f; gc.collect()
    T = np.concatenate(temps, axis=2)
    np.savez_compressed(out / 'global_T.npz', years=YEARS[YEARS >= 1990], cases=np.array(names), members=np.array(ids), T=T)


if __name__ == '__main__':
    run()
