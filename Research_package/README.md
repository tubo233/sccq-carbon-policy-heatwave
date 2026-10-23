# Code and data for the Sichuan–Chongqing carbon policy and heatwave study

This package accompanies the manuscript submitted to Frontiers in Earth Science. It contains the inputs, precomputed ensembles, analysis scripts, and data underlying the included figures and tables. The central policy comparison uses NPI, NPR and NPLS. The supplementary DPEC and SSP experiments address spatial scale, emission species and the timing of emergence.

## Contents of this deposit

This Zenodo deposit contains this README, the Python code, prepared model inputs, ensemble results and tabular outputs, and the prepared daily ERA5-Land temperature arrays (`inputs/weather`). The Supplementary Material of the article (Tables S1–S12 and Figures S1–S5) is published with the article. The scripts do not download data; the external inputs listed under “Not included in this deposit” must be supplied manually.

### Files of this record

The package consists of three zip files, each stored in the record as numbered parts of at most 10 MB (`.001`, `.002`, …):

- `Research_package_main.zip`: everything except the two daily temperature arrays.
- `Research_package_tmax.zip`: `inputs/weather/tmax_float64.npy` (daily maximum temperature).
- `Research_package_tmin.zip`: `inputs/weather/tmin_float64.npy` (daily minimum temperature).

To restore the zip files, download all parts into one folder and run `python join_parts.py`, which joins the parts and checks them against `MD5SUMS.txt`. Alternatively, join them with `copy /b Research_package_main.zip.001+Research_package_main.zip.002+… Research_package_main.zip` on Windows or `cat Research_package_main.zip.* > Research_package_main.zip` on Linux and macOS, or open the `.001` part with 7-Zip. All three zip files unpack into the same `Research_package/` folder; extract them into one directory to obtain the complete package. The two temperature arrays are read only by the `weather` stage; all other stages run from `Research_package_main.zip` alone.

## Data used by individual scripts

- **GCAM-China.** `inputs/gcam/gcam_china_projections_2025_2060.csv` contains the CO₂ emissions of the low, mid and high GCAM-China scenarios for power, industry, buildings and transport, 2025–2060, connected to CEADs emissions in 2020. `models/main/02_scenarios.py` reads these projections directly.
- **CNEMC air-quality observations.** The `air-quality` stage expects city-day series compiled from https://www.cnemc.cn as `models/air_quality/cnemc_sccq_daily.csv` (columns `date, city, pm25, pm25_hours, mda8, o3_8h_hours, no2, no2_hours`; not included, see below). The raw series, the annual summary, the coverage table and Figure S5 are not redistributed under the CNEMC terms; obtaining the data and regenerating these outputs is described below.
- **GADM 4.1 boundaries.** `scripts/03_paper_figures.py` (study-area maps of Figure 1) and `scripts/07_methane_shares.py` (methane allocation from the EDGAR grid) read the level-2 geometry for China, `gadm41_CHN_2.json`, available from https://gadm.org, in `inputs/boundaries/` (not included, see below). Figure 1 is included as an academic map; the boundary file itself is not redistributed. The resulting methane shares are supplied in `data`.

## Not included in this deposit

These inputs are not redistributed because their licences do not allow it. The code that uses them is included so they can be regenerated locally.

### CNEMC air-quality data

- Not included: `models/air_quality/cnemc_sccq_daily.csv` (city-day series), `models/air_quality/heatwave_day_pollution_by_year.csv` (annual summary), `models/air_quality/cnemc_coverage_by_city_year.csv` (coverage table), and `figures/FigureS5_air_quality_heatwave_days.*` (figure derived from these data).
- Why: the CNEMC legal statement limits use to personal use and requires permission for copying or derivative products (https://www.cnemc.cn/flsm/).
- Where to get it: CNEMC national urban air-quality real-time publishing platform (https://air.cnemc.cn:18007/) and CNEMC data pages (https://www.cnemc.cn/sssj/); historical city-day data for 2015–2024 may require a request to CNEMC. The code expects columns `date, city, pm25, pm25_hours, mda8, o3_8h_hours, no2, no2_hours`; place the file at `models/air_quality/cnemc_sccq_daily.csv`.
- How to regenerate locally: `python models/air_quality/03_pollution_channel.py` writes the annual summary, then `python scripts/10_si_data_figures.py fig_s5` draws Figure S5. These outputs are for local use; do not redistribute them without CNEMC permission.

### GADM boundaries

- Not included: `inputs/boundaries/gadm41_CHN_2.json`.
- Why: the GADM licence permits academic and other non-commercial use, including maps for academic publishing, but not redistribution or commercial use without prior permission (https://gadm.org/license.html).
- Where to get it: download the China level-2 file from https://gadm.org/download_country.html and place it at `inputs/boundaries/gadm41_CHN_2.json` for Figure 1 and the methane allocation.
- Figure 1 is included as an academic map; the boundary file itself is not redistributed.

## Environment and execution

Use Python 3.12 and the versions in `requirements.txt`. Install dependencies in an isolated Python environment with `python -m pip install -r requirements.txt`. Run commands from this directory:

```text
python reproduce.py --stage figures
python reproduce.py --stage metrics
python reproduce.py --stage scenarios
python reproduce.py --stage fair
python reproduce.py --stage weather
python reproduce.py --stage air-quality
python reproduce.py --stage scale-fair
```

The `figures` stage needs the GADM file for Figure 1, and `air-quality` needs the CNEMC city-day series; both are described under “Not included in this deposit”.

The default `figures` stage uses the supplied outputs to regenerate Figures 1–5 and Figures S1–S4 (Figure S5 is generated only when the CNEMC-derived summary is supplied locally; Figure 1 requires the external GADM file, see above). The `metrics` stage recomputes the main uncertainty summaries, regional supplementary temperature series and Table S7 from the supplied ensemble temperatures. `scenarios` reconstructs the provincial energy and emission pathways. `fair` reruns the 61 main and 29 supplementary emission cases for all 841 calibrated FaIR members. `weather` regenerates the 128-member heatwave-response ensemble from daily input temperatures and is substantially more expensive than plotting. `air-quality` recomputes the regional annual summary used for Figure S5 from user-supplied CNEMC city-day series; neither the summary nor Figure S5 is redistributed here.

To rebuild the primary chain from prepared inputs, place the GADM file (see above), then run `scenarios`, `fair`, `weather`, `metrics`, then `figures`. If regenerating only climate metrics, the supplied weather-response curves can be reused. Generated results are written inside this directory; keep a separate extracted copy if retaining the supplied results is important.

## Directory guide

| Location | Contents |
| --- | --- |
| `scripts` | Main figure and table generation, methane allocation and supplementary FaIR calculations |
| `models/main` | Main NPI/NPR/NPLS pathway construction, FaIR runs and metrics |
| `models/scales` | Shared DPEC scale calculations and regional noise parameters |
| `models/heatwave` | Weather generator, response curves, ensemble members and observed metrics |
| `models/air_quality` | City heatwave series from ERA5-Land, grid-to-city mapping and the CNEMC processing script (CNEMC data and derived outputs not included) |
| `models/exposure/results` | Scale and committed-exposure summaries used in the figure calculations |
| `weather/scripts/common.py` | Heatwave definitions and whole-year weather-resampling functions |
| `inputs/fair` | FaIR species, 841 parameter sets, harmonised SSP emissions, solar and volcanic forcing, regional coefficients |
| `inputs/emissions` | DPEC paths, the provincial SSP2-com totals of Zhong et al. (`THUCMA_provincial_totals_from_sector_table.csv`; THUCMA is the file label used in this package) and the CEADs historical series |
| `inputs/gcam` | GCAM-China-v6.0 projections for 2025–2060 (low, mid and high scenarios, four sectors), connected to CEADs emissions in 2020 |
| `inputs/EDGAR_v8` | EDGAR v8.0 country totals and the 2020 methane grid used for allocation shares |
| `inputs/climate` | CMIP6 regional monthly changes used in the 2035 background |
| `inputs/weather` | Prepared ERA5-Land Tmax/Tmin arrays and grid metadata |
| `inputs/boundaries` | Expected location for the GADM 4.1 geometry (not included; see above) |
| `data` | Main figure data, table values (`tableS1C_pathways.csv`, `tableS5B_sensitivity_summary.csv`, `tableS7_methane_sensitivity.csv`, `tableS8_experiment_protocol.csv`, `tableS9_full_results.csv`; `table2_tableS5B_S9_values.json` for main Table 2), methane shares and supplementary ensemble results |
| `data/si` | Tables S1–S6, S10–S12 and weather-generator summaries |
| `figures` | Regenerated PNG, vector PDF and TIFF artwork for Figures 1–5 and S1–S4 (Figure S5 not redistributed) |
| `diagrams` | Editable draw.io source for Figure 1C; its native node geometry, labels and connections are read by the figure script |

## Variables and units

- Main scenario CSVs: `region` is Sichuan or Chongqing; `scenario` identifies the NPI construction, NPR, or NPLS pathway; `year` is calendar year; `gdp_index` is real GDP relative to 1978 = 100; `energy_10k_tce` is 10,000 tonnes of coal equivalent; `co2_mt` is Mt CO₂. The low-GDP file retains actual and official growth through 2025 and reduces growth by one percentage point from 2026.
- `models/main/fair/global_T.npz`: `T` has axes year × case × FaIR member; `years`, `cases`, and `members` label those axes. Temperature is in °C. Main reported warming subtracts each member's NPR 2015–2024 mean before regional scaling; it is not an absolute local temperature.
- `data/fair_supplement/global_temperature.npz` has the same labelled axes. Supplementary regional changes are case minus SSP2-4.5 within each member, multiplied by three CMIP6 coefficients. Positive `dTmax` means warming. Figures reporting *cooling* reverse this sign.
- `timeseries_metrics_1990_2060.csv`: `median`, `p05`, and `p95` summarize 841 FaIR × 10 regional-coefficient combinations. W and A are °C; P, C and G are percentages; dHWF is days per year. `view=own` changes only the indicated province.
- `tableS7_methane_sensitivity.csv`: temperature columns are °C; `SC_share_warming_2035` is a fraction; crossover and emergence columns are calendar years. Empty year fields mean the criterion was not reached by 2100.
- `models/heatwave/response_members.npz`: `values` axes are case × weather member × region × threshold definition × heatwave kind × weighting × metric; labels are included in the archive. HWF is days and HWN is events per year.
- `inputs/weather/tmax_float64.npy` and `tmin_float64.npy`: float64 arrays, 12,784 consecutive days from 1990-01-01 to 2024-12-31 × 5,326 grid cells, °C. Column order matches `grid_metadata_population.csv`; coordinates are degrees and `population2020` is persons. `area_weight` is the normalized regional area weight.
- Population-weighted heatwave counts are calculated after cell-level identification, using strict threshold exceedance for at least three consecutive days. The primary threshold is each cell's annual 90th percentile in 1990–2005.
- Some code comments, region labels and column names are in Chinese: 四川 = Sichuan, 重庆 = Chongqing; in `inputs/emissions/ceads_history.csv`, 年份 = year and CO2总排放_百万吨 = total CO₂ emissions (Mt). City names in the air-quality files are the official Chinese names.
- The air-quality files use µg m⁻³ for PM₂.₅ and ozone MDA8.

## Experiment details

The main climate background is SSP2-4.5. The harmonised SSPs are identical through 2022. Global SSP1-1.9 and SSP1-2.6 benchmarks diverge from 2023; their timing variants hold SSP2-4.5 through 2025 and delay each species' emission difference by three years. These timing variants do not establish pathway feasibility.

The main FaIR experiment uses 841 members and 10 regional coefficients; supplementary scale experiments use 841 members and three coefficients. Weather resampling uses 128 members and the original deterministic seed. Table S8 defines perturbation onset and extension rules. The policy-pathway and DPEC benchmark experiments have different constructions and are compared only in sign and order of magnitude.

The package starts from prepared input data and the GCAM-China projections for 2025–2060. The spatial and temporal extracts used by the analysis are supplied, apart from the external inputs listed above; their original sources are listed in `DATA_SOURCES.md`.

## Reproduction checks

Scenario reconstruction and the main climate metrics reproduce the source values within numerical precision. The 61 main and 29 supplementary emission-input cases were rebuilt. The main 841-member FaIR ensemble was rerun with 61 cases (no-policy constructions and single-province runs); the 24 cases common to both runs matched to numerical precision. The figure scripts use package-relative paths; the external inputs are noted above. The supplementary FaIR ensemble and the weather ensemble were not rerun; their precomputed outputs are included.

## Sources and reuse

Use the article and the source citations in `DATA_SOURCES.md` when citing these materials. Third-party datasets and software retain their original terms and attribution requirements. The scripts use package-relative paths and export publication-quality figures; the scientific definitions and numerical inputs are those reported in the manuscript.

The additional `scale-fair` stage regenerates the broader DPEC scale ensemble. Its precomputed ensemble is also supplied, and `metrics` recalculates the exposure and noise summaries from it. This auxiliary experiment uses a CO2-share proxy for methane; the methane sensitivity and Figure 5 use the separate 29-case supplementary ensemble.
