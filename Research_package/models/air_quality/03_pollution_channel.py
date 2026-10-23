"""川渝 22 个地级市热浪期间空气质量（2015-2024，图 S5）。

只做描述性汇总：逐年人口加权热浪人日占比、热浪日与全日的 PM2.5 和 MDA8 臭氧、
MDA8 超二级标准（160 µg/m³）的热浪人日占比，以及臭氧超标日的人口加权占比。
输入：CNEMC 小时数据整理的城市逐日污染物序列、ERA5-Land 城市热浪序列和 2020 年城市人口。
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'common'))
from paths import *
import pandas as pd

HERE = Path(__file__).resolve().parent
O3_LIMIT = 160.0


def load(src):
    p = pd.read_csv(src / 'cnemc_sccq_daily.csv', parse_dates=['date'], encoding='utf-8-sig')
    h = pd.read_csv(src / 'city_daily_heat.csv', parse_dates=['date'], encoding='utf-8-sig')
    pop = pd.read_csv(src / 'city_population2020.csv', encoding='utf-8-sig', index_col=0).iloc[:, 0]
    d = h.merge(p, on=['date', 'city'], how='left')
    d['year'] = d.date.dt.year
    d['pop'] = d.city.map(pop)
    return d


def main(src=HERE, dst=HERE):
    d = load(src)
    d['o3x'] = (d.mda8 > O3_LIMIT).astype(float)
    d = d[d.pm25.notna() & d.mda8.notna()].copy()
    d['w'] = d['pop'] / d.groupby('date')['pop'].transform('sum')
    rows = []
    for y, g in d.groupby('year'):
        hw = g.hw_pop_share * g['pop']
        rows.append(dict(year=y,
                         hw_person_days_share=float((g.hw_pop_share * g.w).sum()),
                         pm25_on_hw=float((g.pm25 * hw).sum() / hw.sum()), pm25_all=float((g.pm25 * g['pop']).sum() / g['pop'].sum()),
                         mda8_on_hw=float((g.mda8 * hw).sum() / hw.sum()), mda8_all=float((g.mda8 * g['pop']).sum() / g['pop'].sum()),
                         compound_share_of_hw=float((g.o3x * hw).sum() / hw.sum()),
                         o3_exceed_days_popw=float((g.o3x * g.w).sum())))
    desc = pd.DataFrame(rows)
    desc.to_csv(dst / 'heatwave_day_pollution_by_year.csv', index=False, encoding='utf-8-sig')
    print(desc.round(3).to_string())


if __name__ == '__main__':
    main()
