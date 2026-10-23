"""EDGAR v8.0（2020）的甲烷分配份额：中国占全球人为 CH4 的比例，以及川渝占中国的比例
（0.1° 网格，格点中心落在 GADM 4.1 省界内）。models/scales 的尺度实验用 CO2 份额代理甲烷。
"""
import json
import shutil
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from shapely import contains_xy
from shapely.geometry import shape
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
ED = HERE / 'inputs/EDGAR_v8'
OUT = HERE / 'data'


def main():
    t = pd.read_excel(ED / 'EDGAR_CH4_1970_2022.xlsx', 'TOTALS BY COUNTRY', header=9)
    world = t['Y_2020'].sum() / 1e3                                   # Gg 转 Mt
    china_tab = t.loc[t.Country_code_A3.eq('CHN'), 'Y_2020'].sum() / 1e3

    # netCDF4 在 Windows 上打不开含非 ASCII 字符的路径
    tmp = Path(tempfile.mkdtemp()) / 'edgar.nc'
    shutil.copy(ED / 'v8.0_FT2022_GHG_CH4_2020_TOTALS_emi.nc', tmp)
    e = xr.open_dataset(tmp).emissions.load()                          # 每格 CH4（t）
    gj = json.load(open(HERE / 'inputs/boundaries/gadm41_CHN_2.json', encoding='utf-8'))
    prov = {}
    for f in gj['features']:
        prov.setdefault(f['properties']['NAME_1'], []).append(shape(f['geometry']))
    prov = {k: unary_union(v) for k, v in prov.items()}
    china = unary_union(list(prov.values()))
    sub = e.sel(lat=slice(15, 55), lon=slice(70, 140))
    lon, lat = np.meshgrid(sub.lon.values, sub.lat.values)
    v = sub.values
    out = dict(world_Mt=world, china_table_Mt=china_tab, china_share_of_world=china_tab / world)
    out['china_grid_Mt'] = float(v[contains_xy(china, lon, lat)].sum() / 1e6)
    for k in ['Sichuan', 'Chongqing']:
        out[f'{k}_Mt'] = float(v[contains_xy(prov[k], lon, lat)].sum() / 1e6)
    sccq = out['Sichuan_Mt'] + out['Chongqing_Mt']
    out['sccq_Mt'] = sccq
    out['sccq_share_of_china_grid'] = sccq / out['china_grid_Mt']
    out['sccq_share_of_china_table'] = sccq / china_tab
    pd.Series(out).to_csv(OUT / 'methane_shares_EDGARv8_2020.csv', header=['value'])
    print(pd.Series(out).round(4).to_string())


if __name__ == '__main__':
    main()
