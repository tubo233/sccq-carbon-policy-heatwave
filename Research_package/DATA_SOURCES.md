# Data sources

| Source | Version or scope | Access and use in this package |
| --- | --- | --- |
| ERA5-Land | Daily Tmax/Tmin, 1990–2024, Sichuan–Chongqing 0.1° grid | https://doi.org/10.24381/cds.e2161bac ; prepared regional daily arrays supplied in `inputs/weather` |
| WorldPop | 2020 population counts | https://www.worldpop.org ; rebinned grid population supplied with weather metadata and figure data |
| CEADs | Provincial CO₂ history | https://www.ceads.net ; historical series used for pathway anchoring |
| FaIR | v2.2.4 | https://github.com/OMS-NetZero/FAIR ; installed through requirements.txt |
| fair-calibrate | v1.4.1 | https://github.com/chrisroadmap/fair-calibrate/tree/v1.4.1 ; harmonised SSPs, calibrated parameters and natural forcing supplied |
| CMIP6 | GFDL-ESM4, MPI-ESM1-2-LR, UKESM1-0-LL | https://esgf-node.llnl.gov/projects/cmip6/ ; regional coefficients and monthly SSP2-4.5 changes supplied |
| DPEC | Cheng et al. emissions dataset | https://doi.org/10.5281/zenodo.7383239 ; national pathway CSVs supplied |
| Zhong et al. SSP2-com (file label THUCMA) | Province-level SSP2-com emission dataset; 2020 provincial totals used for the allocation shares | https://doi.org/10.5281/zenodo.20728925 ; fixed-version provincial SSP2-com totals aggregated from the sector table (`inputs/emissions/THUCMA_provincial_totals_from_sector_table.csv`), used for 2020 allocation and pathway calculations; related paper: https://doi.org/10.1038/s41467-025-62983-5 |
| EDGAR | v8.0, 2020 methane allocation | https://edgar.jrc.ec.europa.eu/dataset_ghg80 ; country table and 2020 grid supplied |
| Crippa et al. | GHG emissions of all world countries, 2023, EUR 31658 EN, JRC134504 | https://doi.org/10.2760/953322 (online edition); https://publications.jrc.ec.europa.eu/repository/handle/JRC134504 |
| GCAM-China | Laboratory scenario projections used in this article | Projections for 2025–2060 (low, mid, high), connected to CEADs emissions in 2020, supplied in inputs/gcam; model background: Calvin et al. (2019), https://doi.org/10.5194/gmd-12-677-2019 |
| Sichuan statistics | Statistical yearbooks and 2023/2024 GDP growth releases | https://tjj.sc.gov.cn ; extracted historical data supplied in models/main/history_gdp_energy_co2.csv |
| Chongqing statistics | Statistical yearbooks | https://tjj.cq.gov.cn ; extracted historical data supplied in models/main/history_gdp_energy_co2.csv |
| CNEMC | City air-quality observations, 2015–2024 | https://www.cnemc.cn/sssj/ (real-time platform: https://air.cnemc.cn:18007/); city-day series, derived annual summary, coverage table and Figure S5 are not redistributed (CNEMC legal statement: https://www.cnemc.cn/flsm/); obtain the city-day data from CNEMC and place it at models/air_quality/cnemc_sccq_daily.csv |
| GADM | v4.1 China administrative geometry | https://gadm.org/download_country.html ; `gadm41_CHN_2.json` expected in inputs/boundaries (not redistributed; GADM licence: https://gadm.org/license.html); download the China level-2 file; Figure 1 is included as an academic map |

## Sichuan GDP growth

The 2023 communiqué reports 6.0% real GDP growth; the 2024 communiqué reports 5.7%. These are the published growth rates used for the study's 2023–2024 continuation; they are not recomputed from nominal GDP levels.

- Sichuan Provincial Bureau of Statistics (2024), Statistical Communiqué of Sichuan Province on the 2023 National Economic and Social Development: https://tjj.sc.gov.cn/scstjj/c111701/2024/3/14/f403a921ad204ecfaecde2866aec3aac.shtml
- Sichuan Provincial Bureau of Statistics (2025), Statistical Communiqué of Sichuan Province on the 2024 National Economic and Social Development: https://tjj.sc.gov.cn/scstjj/c112117/2025/3/17/159a4f075f5740da8d4f75dd09d511af.shtml

## Policy documents used for the provincial targets

- Sichuan energy conservation work plan, Chuan Fu Fa [2022] No. 20: https://www.sc.gov.cn/10462/zfwjts/2022/7/21/54cc17aea9164de6a57b361367782cf5.shtml
- Chongqing energy conservation implementation plan, Yu Fu Fa [2022] No. 39: https://www.cq.gov.cn/zwgk/zfxxgkml/zfgb/2022/d18q/202210/t20221020_11207120.html
- Chongqing 14th Five-Year energy development plan, Yu Fu Ban Fa [2022] No. 48: https://www.cq.gov.cn/zwgk/zfxxgkml/szfwj/qtgw/202206/t20220615_10818266.html
- Sichuan 12th Five-Year energy development plan, Chuan Ban Fa [2011] No. 65: https://www.sc.gov.cn/10462/10883/11066/2011/10/28/10191801.shtml

The complete policy catalogue is Table S1 in the supplementary document; its two parts are also included as CSV files in data/si.
