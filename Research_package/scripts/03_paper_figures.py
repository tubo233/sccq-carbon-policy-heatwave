"""图 1–5（英文标签、白底、300 dpi PNG + 矢量 PDF）。

所有数字来自主情景/补充实验输出，以及 01_paper_numbers.py 写到 data/ 的表。
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import FancyBboxPatch
import matplotlib.patheffects as pe
from shapely.geometry import shape
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
EXP = HERE / 'models'
D = EXP / 'main'
DATA = HERE / 'data'
FIG = HERE / 'figures'; FIG.mkdir(exist_ok=True)

rcParams.update({'font.family': ['Arial', 'Segoe UI'],  # Segoe UI 只补 Arial 缺的下标/上标字形（₂ ₄ ⁻）
                 'font.size': 7.5, 'axes.unicode_minus': True, 'axes.linewidth': 0.6,
                 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
COL = {'NPI': '#eb6834', 'NPR': '#2a78d6', 'NPLS': '#1baf7a'}
COOL, WARM = '#2a78d6', '#e34948'
ORD = ['#86b6ef', '#3987e5', '#1c5cab', '#0d366b']
BLUE_SEQ = LinearSegmentedColormap.from_list('blue', ['#e6f0fc', '#9ec5f4', '#3987e5', '#1c5cab', '#0d366b'])
ORANGE_SEQ = LinearSegmentedColormap.from_list('orange', ['#fdeee6', '#f6b894', '#eb6834', '#b8461a', '#6e2408'])
REGN = {'sichuan': 'Sichuan', 'chongqing': 'Chongqing'}
CN = {'sichuan': '四川', 'chongqing': '重庆'}


def style(ax, grid_axis='y'):
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=7)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, lw=0.5); ax.set_axisbelow(True)


def tag(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s.upper(), transform=ax.transAxes, fontsize=9, fontweight='bold', color=INK, ha='right', va='bottom')


def save(fig, name):
    from matplotlib.text import Text
    from PIL import Image
    if not name.startswith('FigureS'):
        frontiers_style(fig)
    for t in fig.findobj(match=Text):
        if t.get_text():
            t.set_fontsize(max(8, t.get_fontsize()))
    fig.savefig(FIG / f'{name}.png', dpi=300, facecolor='white')
    fig.savefig(FIG / f'{name}.pdf', facecolor='white')
    Image.open(FIG / f'{name}.png').convert('RGB').save(FIG / f'{name}.tiff', compression='tiff_lzw', dpi=(300,300))
    plt.close(fig)


def frontiers_style(fig):
    """只改版式：180 mm 宽，保留所有科学序列和轴范围。"""
    from matplotlib.text import Text
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.collections import Collection
    from matplotlib.ticker import FuncFormatter, NullLocator
    import re
    superscript = str.maketrans('0123456789-+', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺')
    def log_label(value, _position):
        if value <= 0:
            return ''
        power = int(np.round(np.log10(value)))
        return ('10' + str(power).translate(superscript)) if np.isclose(value, 10.**power) else f'{value:g}'
    for ax in fig.axes:
        ax.tick_params(axis='both', which='major', labelsize=9, width=2, length=3, pad=4)
        ax.tick_params(axis='both', which='minor', width=2, length=0)
        if ax.get_xscale() == 'log':
            ax.xaxis.set_major_formatter(FuncFormatter(log_label))
            ax.xaxis.set_minor_locator(NullLocator())
        if ax.get_yscale() == 'log':
            ax.yaxis.set_major_formatter(FuncFormatter(log_label))
            ax.yaxis.set_minor_locator(NullLocator())
        for spine in ax.spines.values():
            spine.set_linewidth(2)
        for gridline in ax.get_xgridlines() + ax.get_ygridlines():
            gridline.set_color('#ededeb')
            gridline.set_linewidth(2)
    # 先固定刻度再设字号，否则导出时懒创建刻度会把 6.5/7 pt 的旧设置带回来。
    fig.canvas.draw()
    for text in fig.findobj(match=Text):
        label = text.get_text()
        if not label:
            continue
        label = label.replace('CO$_2$', 'CO₂').replace('CH$_4$', 'CH₄').replace('λ$_r$', 'λᵣ')
        label = re.sub(r'\$\^\{([+-]?\d+)\}\$', lambda m: m.group(1).translate(superscript), label)
        text.set_text(label)
        text.set_fontsize(max(9, text.get_fontsize()))
        if label in ('A', 'B', 'C', 'D'):
            text.set_fontsize(11)
        if '$' not in label:
            text.set_parse_math(False)
    for line in fig.findobj(match=Line2D):
        if line.get_linewidth() > 0:
            line.set_linewidth(max(2, line.get_linewidth()))
        if line.get_marker() not in (None, 'None', '', ' '):
            line.set_markersize(max(5, line.get_markersize()))
            line.set_markeredgewidth(2)
    for patch in fig.findobj(match=Patch):
        if patch.get_linewidth() > 0:
            patch.set_linewidth(max(2, patch.get_linewidth()))
    for collection in fig.findobj(match=Collection):
        widths = collection.get_linewidths()
        if len(widths):
            collection.set_linewidths(np.where(np.asarray(widths) > 0, np.maximum(widths, 2), widths))
    fig.canvas.draw()


# ---------------------------------------------------------------- 图 1：研究区与框架
def fig1():
    g = pd.read_csv(DATA / 'grid_hwf_population.csv')
    gj = json.load(open(HERE / 'inputs/boundaries/gadm41_CHN_2.json', encoding='utf-8'))
    prov = {}
    for f in gj['features']:
        prov.setdefault(f['properties']['NAME_1'], []).append(shape(f['geometry']))
    prov = {k: unary_union(v) for k, v in prov.items()}
    lat = np.round(g.latitude, 1); lon = np.round(g.longitude, 1)
    la = np.round(np.arange(lat.min(), lat.max() + 0.05, 0.1), 1); lo = np.round(np.arange(lon.min(), lon.max() + 0.05, 0.1), 1)
    ii = np.searchsorted(la, lat); jj = np.searchsorted(lo, lon)
    cell_km2 = (0.1 * 111.32) ** 2 * np.cos(np.deg2rad(lat))
    dens = np.full((len(la), len(lo)), np.nan); hw = dens.copy()
    dens[ii, jj] = np.maximum(g.population2020 / cell_km2, 0.1); hw[ii, jj] = g.hwf_2015_2024
    fig = plt.figure(figsize=(180 / 25.4, 4.7))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.0, 0.8], hspace=0.32, wspace=0.24, left=0.075, right=0.98, top=0.95, bottom=0.025)
    edges = (np.r_[lo - 0.05, lo[-1] + 0.05], np.r_[la - 0.05, la[-1] + 0.05])
    for k, (arr, cmap, norm, lab) in enumerate([(dens, BLUE_SEQ, LogNorm(1, 5000), 'Population density (persons km$^{-2}$)'),
                                                (hw, ORANGE_SEQ, plt.Normalize(20, 70), 'Heatwave days (days yr$^{-1}$)')]):
        ax = fig.add_subplot(gs[0, k])
        for name, geom in prov.items():
            if name in ('Sichuan', 'Chongqing'):
                continue
            if geom.intersects(unary_union([prov['Sichuan'], prov['Chongqing']]).buffer(3)):
                for poly in getattr(geom, 'geoms', [geom]):
                    ax.fill(*poly.exterior.xy, color='#f3f2ee', lw=0.3, ec='#d9d8d2', zorder=0)
        m = ax.pcolormesh(*edges, arr, cmap=cmap, norm=norm, shading='flat', zorder=1, rasterized=True)
        for name in ('Sichuan', 'Chongqing'):
            for poly in getattr(prov[name], 'geoms', [prov[name]]):
                ax.plot(*poly.exterior.xy, color=INK, lw=0.7, zorder=2)
        ax.set_xlim(97.2, 110.4); ax.set_ylim(25.9, 34.5); ax.set_aspect(1 / np.cos(np.deg2rad(30)))
        ax.set_xticks([98, 102, 106, 110]); ax.set_yticks([26, 28, 30, 32, 34])
        ax.set_xticklabels([f'{v}°E' for v in ax.get_xticks().astype(int)]); ax.set_yticklabels([f'{v}°N' for v in ax.get_yticks().astype(int)])
        ax.tick_params(labelsize=6.5, colors=INK2, length=2)
        for s in ax.spines.values():
            s.set_color(AXIS); s.set_linewidth(0.6)
        halo = [pe.withStroke(linewidth=2.2, foreground='white')]
        ax.text(100.6, 31.6, 'Sichuan', fontsize=8, color=INK, ha='center', fontweight='bold', path_effects=halo, zorder=4)
        ax.text(108.3, 30.75, 'Chongqing', fontsize=8, color=INK, ha='center', fontweight='bold', path_effects=halo, zorder=4)
        if k == 0:
            ax.plot(104.07, 30.57, 'o', ms=3, mfc='white', mec=INK, mew=0.8, zorder=4)
            ax.text(103.9, 30.75, 'Chengdu', fontsize=6.5, color=INK, ha='right', path_effects=halo, zorder=4)
        cb = fig.colorbar(m, ax=ax, orientation='horizontal', fraction=0.05, pad=0.08, aspect=35)
        cb.set_label(lab, fontsize=7, color=INK2); cb.ax.tick_params(labelsize=6.5, colors=INK2, length=2); cb.outline.set_edgecolor(AXIS)
        tag(ax, 'ab'[k], x=0.0, y=1.01)
    # 把 .drawio 里的框和标签画成矢量，改图后这里会跟着变。
    import xml.etree.ElementTree as ET
    import re, html
    from matplotlib.patches import Rectangle
    root = ET.parse(HERE / 'diagrams/figure1-framework.drawio').getroot()
    ax = fig.add_subplot(gs[1, :]); ax.axis('off')
    ax.set_xlim(0, 902); ax.set_ylim(122, 10)
    tag(ax, 'c', x=0.0, y=1.01)
    nodes = {}
    for cell in root.findall('.//mxCell[@vertex="1"]'):
        geom = cell.find('mxGeometry')
        x, y, w, h = [float(geom.get(k)) for k in ('x','y','width','height')]
        nodes[cell.get('id')] = (x,y,w,h)
        ax.add_patch(Rectangle((x,y),w,h,fc='white',ec='#49708A',lw=0.8))
        parts = re.split(r'<br\s*/?>',cell.get('value'))
        labels = [html.unescape(re.sub('<[^>]+>','',v)) for v in parts]
        ax.text(x+w/2,y+h*.35,labels[0].replace('CO₂',r'CO$_2$'),ha='center',va='center',fontsize=8.6,fontweight='bold',color=INK)
        ax.text(x+w/2,y+h*.69,labels[1],ha='center',va='center',fontsize=8,color=INK2)
    for cell in root.findall('.//mxCell[@edge="1"]'):
        x,y,w,h=nodes[cell.get('source')]; tx,ty,tw,th=nodes[cell.get('target')]
        ax.annotate('',xy=(tx-3,ty+th/2),xytext=(x+w+3,y+h/2),arrowprops=dict(arrowstyle='-|>',color='#52616B',lw=2,mutation_scale=10))
    save(fig, 'Figure1_study_area_framework')


# ---------------------------------------------------------------- 图 2：能源与 CO2 路径
def fig2():
    from matplotlib.lines import Line2D
    sp = pd.read_csv(D / 'scenario_paths_energy_co2.csv')
    fig, axs = plt.subplots(2, 2, figsize=(180 / 25.4, 5.3), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.975, top=0.875, bottom=0.09, hspace=0.28, wspace=0.30)
    yy = np.arange(2005, 2036)
    keys={'NPI':'NPI_aeei1.0','NPR':'NPR','NPLS':'NPLS_gcam_mid'}
    for j,r in enumerate(['sichuan','chongqing']):
        s=sp[sp.region.eq(r)]
        for i,(var,unit,scale) in enumerate([('energy_10k_tce','Energy use (Mtce)',0.01),('co2_mt','CO$_2$ emissions (Mt yr$^{-1}$)',1)]):
            ax=axs[i,j]; style(ax)
            p=s.pivot_table(index='year',columns='scenario',values=var).loc[yy]*scale
            ax.axvspan(2005,2022,color='#f4f5f6',lw=0,zorder=0)
            ax.axvline(2022,color=AXIS,lw=0.6,ls=':')
            for k in keys:
                years=yy if k=='NPI' else yy[yy>=2022]
                ax.plot(years,p.loc[years,keys[k]],color=COL[k],lw=1.8,
                        ls=(0,(4,2)) if k=='NPLS' else '-',zorder=4 if k=='NPLS' else 3)
            ax.plot(yy[yy<=2022],p.loc[yy[yy<=2022],'NPR'],color=INK2,lw=1.3,zorder=5)
            top=p[list(keys.values())].to_numpy().max()*1.12
            ax.set_ylim(0,top); ax.set_xlim(2005,2035)
            ax.set_xticks([2005,2015,2025,2035])
            if j==0: ax.set_ylabel(unit,color=INK2)
            if i==0:
                ax.set_title(REGN[r],fontsize=9,fontweight='bold',pad=10)
                ax.text(2006,top*.94,'Observed',color=MUTED,va='top')
                ax.text(2023,top*.94,'Projected',color=MUTED,va='top')
            tag(ax,'abcd'[i*2+j],x=0,y=1.02)
    h=[Line2D([],[],color=INK2,lw=1.3,label='Observed')]
    h += [Line2D([],[],color=COL[k],lw=1.8,ls=(0,(4,2)) if k=='NPLS' else '-',label=k) for k in keys]
    fig.legend(handles=h,ncol=4,loc='upper center',bbox_to_anchor=(.53,.995),frameon=False,columnspacing=2.8)
    save(fig,'Figure2_energy_CO2_pathways')


# ---------------------------------------------------------------- 图 3：增温与避免增温
def fig3():
    ts = pd.read_csv(DATA / 'timeseries_metrics_1990_2060.csv')
    t = ts[(ts.view == 'combined')]
    get = lambda r, k: t[(t.region == r) & (t.metric == k)].set_index('year').loc[2006:2060]
    fig, axs = plt.subplots(2, 2, figsize=(180 / 25.4, 5.7))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.93, bottom=0.08, hspace=0.35, wspace=0.32)
    ext = lambda ax: (ax.axvspan(2035.5, 2060, color='#f6f6f3', lw=0, zorder=0), ax.axvline(2022.5, color=AXIS, lw=0.6, zorder=0))
    # (a) 增温 W
    ax = axs[0, 0]; style(ax); ext(ax)
    for k in ['NPI', 'NPR', 'NPLS']:
        d = get('sichuan', f'W_{k}')
        if k == 'NPI':
            ax.fill_between(d.index, d.p05, d.p95, color=COL[k], alpha=0.12, lw=0)
        ax.plot(d.index, d['median'], color=COL[k], lw=1.6 if k != 'NPLS' else 1.2, ls='-' if k != 'NPLS' else (0, (3, 1.5)), label=k)
    ax.axhline(0, color=AXIS, lw=0.6)
    ax.set_xlim(2006, 2060); ax.set_ylabel('Warming W relative to 2015–2024 (°C)', fontsize=7.5, color=INK2)
    ax.legend(frameon=False, fontsize=7, loc='upper left'); tag(ax, 'a')
    ax.text(2059, -0.5, 'Sichuan; band: NPI 5–95%\n(FaIR × λ$_r$); NPR and NPLS\noverlap NPI at this scale', fontsize=6.3, color=INK2, ha='right')
    ax.set_ylim(-0.6, 2.4)
    # (b) 避免增温 A
    ax = axs[0, 1]; style(ax); ext(ax)
    for k in ['NPR', 'NPLS']:
        d = get('sichuan', f'A_NPI_{k}')
        ax.fill_between(d.index, d.p05, d.p95, color=COL[k], alpha=0.12, lw=0)
        ax.plot(d.index, d['median'], color=COL[k], lw=1.6, label=f'NPI − {k}')
    d = get('sichuan', 'A_NPR_NPLS')
    ax.plot(d.index, d['median'], color=INK2, lw=1.0, ls=(0, (1, 1.2)), label='NPR − NPLS (local increment)')
    ax.set_xlim(2006, 2060); ax.set_ylim(0, None); ax.set_ylabel('Avoided warming A (°C)', fontsize=7.5, color=INK2)
    ax.legend(frameon=False, fontsize=7, loc='upper left'); tag(ax, 'b')
    for y in [2022, 2035]:
        v = get('sichuan', 'A_NPI_NPLS').loc[y, 'median']
        ax.plot(y, v, 'o', ms=3.5, color=COL['NPLS'], mec='white', mew=0.8, zorder=4)
        ax.annotate(f'{v:.4f} °C', (y, v), xytext=(-4, 9), textcoords='offset points', fontsize=6.5, color=INK, ha='right')
    # (c) 比例 P
    ax = axs[1, 0]; style(ax); ext(ax)
    for k in ['NPR', 'NPLS']:
        d = get('sichuan', f'P_{k}').loc[2030:]
        ax.fill_between(d.index, d.p05, d.p95, color=COL[k], alpha=0.12, lw=0)
        ax.plot(d.index, d['median'], color=COL[k], lw=1.6, label=f'P, {k}')
    ax.set_xlim(2006, 2060); ax.set_ylim(0, 5); ax.set_ylabel('Share of NPI warming avoided, P (%)', fontsize=7.5, color=INK2)
    ax.text(2008, 0.35, 'P is shown from 2030: near the\n2015–2024 reference period its\ndenominator approaches zero', fontsize=6.3, color=INK2)
    ax.legend(frameon=False, fontsize=7, loc='upper left'); tag(ax, 'c')
    # (d) 贡献比例 C
    ax = axs[1, 1]; style(ax); ext(ax)
    for k in ['NPR', 'NPLS']:
        d = get('sichuan', f'C_{k}').loc[2008:]
        ax.fill_between(d.index, d.p05, d.p95, color=COL[k], alpha=0.12, lw=0)
        ax.plot(d.index, d['median'], color=COL[k], lw=1.6, label=f'C, {k}')
    ax.set_xlim(2006, 2060); ax.set_ylim(0, 100); ax.set_ylabel('Own-emission warming avoided, C (%)', fontsize=7.5, color=INK2)
    for y in [2035, 2060]:
        v = get('sichuan', 'C_NPLS').loc[y, 'median']
        ax.plot(y, v, 'o', ms=3.5, color=COL['NPLS'], mec='white', mew=0.8, zorder=4)
        ax.annotate(f'{v:.0f}%', (y, v), xytext=(0, 7), textcoords='offset points', fontsize=6.5, color=INK, ha='center')
    ax.legend(frameon=False, fontsize=7, loc='lower right'); tag(ax, 'd')
    for ax in axs.flat:
        lo, hi = ax.get_ylim()
        ax.text(2059.5, lo + 0.985 * (hi - lo) if ax is not axs[0, 0] else 2.35, '2036–2060 extension', fontsize=6, color=MUTED, va='top', ha='right')
    save(fig, 'Figure3_warming_avoided')


# ---------------------------------------------------------------- 图 4：敏感性（单因素 + 联合）
def sensitivity_table():
    """图 4 和表 3 用的行：气候行来自 models/main/04_metrics.py，构造行每次只变一个因素，
    其余保持主设定，取 FaIR 和 λ 的中位响应。
    """
    m = pd.read_csv(D / 'results/metrics_FaIR_lambda_sensitivity.csv', encoding='utf-8-sig')
    cs = pd.read_csv(D / 'results/metrics_construction_sensitivity.csv', encoding='utf-8-sig')
    rows = []
    for r in ['sichuan', 'chongqing']:
        k = cs[(cs.region == CN[r]) & (cs.year == 2035) & ~cs.NPI.str.contains('lowGDP')]
        npls_var = k[k.NPI.eq('NPI_aeei1.0')]
        aeei = k[k.NPI.str.fullmatch(r'NPI_aeei\d\.\d') & k.NPLS.eq('NPLS_gcam_mid')]   # AEEI 0.5–5%/年，碳强度因子冻结
        cfobs = k[k.NPI.str.contains('cfobs') & k.NPLS.eq('NPLS_gcam_mid')]             # 碳强度因子取观测到 2022
        cqcf = k[k.NPI.str.contains('cqcf') & k.NPLS.eq('NPLS_gcam_mid')]               # 重庆 2005 碳强度因子 ±10%
        allnpi = k[k.NPLS.eq('NPLS_gcam_mid')]
        low = cs[(cs.region == CN[r]) & (cs.year == 2035) & cs.NPI.str.contains('lowGDP')].iloc[0]
        central = k[k.NPI.eq('NPI_aeei1.0') & k.NPLS.eq('NPLS_gcam_mid')].iloc[0]
        for metric in ['A_NPI_NPLS', 'A_NPR_NPLS', 'P_NPLS', 'C_NPLS']:
            d = m[(m.region == CN[r]) & (m.year == 2035) & (m.metric == metric)].iloc[0]
            lam = (np.nan, np.nan) if metric in ('P_NPLS', 'C_NPLS') else (d.lambda_only_min, d.lambda_only_max)
            den = (np.nan, np.nan)
            if metric == 'P_NPLS':
                c = m[(m.region == CN[r]) & (m.year == 2035) & (m.metric == 'P_NPLS_cmipdenom')].iloc[0]
                den = (c.p05, c.p95)
            indep_npi = metric == 'A_NPR_NPLS'
            rows += [dict(region=r, metric=metric, source='central (FaIR × λ median)', lo=d['median'], hi=d['median']),
                     dict(region=r, metric=metric, source='FaIR members (5–95%)', lo=d.fair_only_p05, hi=d.fair_only_p95),
                     dict(region=r, metric=metric, source='λ values (min–max)', lo=lam[0], hi=lam[1]),
                     dict(region=r, metric=metric, source='CMIP6 regional-warming denominator (5–95%)', lo=den[0], hi=den[1]),
                     dict(region=r, metric=metric, source='FaIR × λ (5–95%)', lo=d.p05, hi=d.p95),
                     dict(region=r, metric=metric, source='central (median FaIR, median CMIP6 λ)', lo=central[metric], hi=central[metric]),
                     dict(region=r, metric=metric, source='NPLS: GCAM low–high and provincial targets', lo=npls_var[metric].min(), hi=npls_var[metric].max()),
                     dict(region=r, metric=metric, source='NPI: AEEI 0.5–5%', lo=np.nan if indep_npi else aeei[metric].min(), hi=np.nan if indep_npi else aeei[metric].max()),
                     dict(region=r, metric=metric, source='NPI: carbon factor observed to 2022', lo=np.nan if indep_npi else cfobs[metric].min(), hi=np.nan if indep_npi else cfobs[metric].max()),
                     dict(region=r, metric=metric, source='NPI: CQ 2005 carbon factor ±10%', lo=np.nan if indep_npi else cqcf[metric].min(), hi=np.nan if indep_npi else cqcf[metric].max()),
                     dict(region=r, metric=metric, source='NPI: all constructions', lo=np.nan if indep_npi else allnpi[metric].min(), hi=np.nan if indep_npi else allnpi[metric].max()),
                     dict(region=r, metric=metric, source='GDP growth −1 percentage point', lo=low[metric], hi=low[metric])]
    t = pd.DataFrame(rows)
    t.to_csv(DATA / 'sensitivity_table_2035.csv', index=False)
    return t


def fig4():
    t = sensitivity_table()
    labels = [('FaIR members (5–95%)', 'FaIR members (5–95%)'),
              ('λ values (min–max)', 'λᵣ values\n(min–max, n = 10)'),
              ('FaIR × λ (5–95%)', 'FaIR × λᵣ jointly (5–95%)'),
              ('NPLS: GCAM low–high and provincial targets', 'NPLS: GCAM low–high,\nprovincial targets'),
              ('NPI: AEEI 0.5–5%', 'NPI: AEEI 0.5–5% yr⁻¹'),
              ('NPI: carbon factor observed to 2022', 'NPI: carbon factor observed\nto 2022 (AEEI 1%, 4%)'),
              ('NPI: all constructions', 'NPI: all constructions')]
    fig, axs = plt.subplots(1, 3, figsize=(180 / 25.4, 5.0))
    fig.subplots_adjust(left=0.32, right=0.975, top=0.84, bottom=0.16, wspace=0.25)
    shades = {'sichuan': '#0072B2', 'chongqing': '#D55E00'}
    panels = [('A_NPI_NPLS', 'Avoided warming (°C)\nNPI − NPLS', 'log', (2e-4, 3e-2)),
              ('A_NPR_NPLS', 'Local increment (°C)\nNPR − NPLS', 'log', (2e-5, 1e-3)),
              ('P_NPLS', 'P, NPLS (%)\nof NPI warming', 'lin', (0, 3.3))]
    n = len(labels)
    for p, (metric, xl, scale, xlim) in enumerate(panels):
        ax = axs[p]; style(ax, 'x')
        ax.axhspan(n - 3.5, n - 0.5, color='#eeede7', lw=0, zorder=0)
        for off, r in [(0.15, 'sichuan'), (-0.15, 'chongqing')]:
            s = t[(t.region == r) & (t.metric == metric)].set_index('source')
            cen = s.loc['central (FaIR × λ median)', 'lo']
            for i, (key, _) in enumerate(labels):
                y = n - 1 - i + off
                lo, hi = s.loc[key, 'lo'], s.loc[key, 'hi']
                note = None
                if metric == 'P_NPLS' and key == 'λ values (min–max)':
                    lo, hi = s.loc['CMIP6 regional-warming denominator (5–95%)', ['lo', 'hi']]
                    note = 'λ$_r$ cancels; bar:\nCMIP6 denominator'
                if np.isnan(lo):
                    if r == 'sichuan':
                        ax.text(np.sqrt(xlim[0] * xlim[1]) if scale == 'log' else sum(xlim) / 2, n - 1 - i, 'independent\nof NPI',
                                fontsize=6.3, color=MUTED, ha='center', va='center')
                    continue
                ax.plot([lo, hi], [y, y], color=shades[r], lw=2.6, solid_capstyle='round')
                if i < 3:
                    ax.plot(cen, y, 'o', ms=3.2, mfc='white', mec=shades[r], mew=0.9, zorder=3)
        if scale == 'log':
            ax.set_xscale('log')
        ax.set_xlim(*xlim)
        ax.set_yticks(range(n)); ax.set_yticklabels([x[1] for x in labels][::-1] if p == 0 else [], fontsize=6.8, color=INK)
        ax.set_ylim(-0.6, n - 0.4); ax.tick_params(axis='y', length=0)
        ax.set_xlabel(xl, fontsize=7.2, color=INK2)
        tag(ax, 'abc'[p], x=0.0, y=1.02)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=shades['sichuan'], lw=2.6, label='Sichuan'), Line2D([], [], color=shades['chongqing'], lw=2.6, label='Chongqing'),
         Line2D([], [], marker='o', ls='', ms=3.2, mfc='white', mec=INK, label='Central estimate'),
         plt.Rectangle((0, 0), 1, 1, fc='#eeede7', ec='none', label='climate-response rows')]
    fig.legend(handles=h, ncol=2, loc='upper left', bbox_to_anchor=(0.21, 1.0), frameon=False, fontsize=6.8, columnspacing=1.5)
    save(fig, 'Figure4_sensitivity')


# ---------------------------------------------------------------- 图 5：热浪背景、尺度与物种、显现
def fig5():
    from matplotlib.lines import Line2D
    a2 = pd.read_csv(EXP / 'exposure/results/A2_committed_exposure_2035.csv', encoding='utf-8-sig')
    sup = pd.read_csv(DATA / 'fair_supplement/regional_timeseries_2023_2100.csv')
    b1 = sup[sup.region == 'sichuan']
    m = pd.read_csv(D / 'results/metrics_FaIR_lambda_sensitivity.csv', encoding='utf-8-sig')
    fig = plt.figure(figsize=(180 / 25.4, 7.4))
    # (a) 已锁定 vs 避免的热浪天数
    ax = fig.add_axes([0.255, 0.6, 0.22, 0.36]); style(ax, 'x')
    items = []
    for kind, lab in [('day', 'daytime'), ('night', 'night-time'), ('compound', 'compound')]:
        items.append((f'Committed, {lab}', {r: a2[(a2.region == CN[r]) & (a2.kind == kind)].dHWF for r in ['sichuan', 'chongqing']}, 'range'))
    for met, lab in [('dHWF_NPI_NPR', 'Avoided, NPI − NPR'), ('dHWF_NPI_NPLS', 'Avoided, NPI − NPLS'), ('dHWF_NPR_NPLS', 'Local, NPR − NPLS')]:
        items.append((lab, {r: m[(m.region == CN[r]) & (m.year == 2035) & (m.metric == met)].iloc[0] for r in ['sichuan', 'chongqing']}, 'q'))
    n = len(items)
    ax.axhspan(n - 3.5, n - 0.5, color='#f6f6f3', lw=0, zorder=0)
    for i, (lab, vals, kind) in enumerate(items):
        for off, r, c in [(0.15, 'sichuan', '#0072B2'), (-0.15, 'chongqing', '#D55E00')]:
            v = vals[r]; y = n - 1 - i + off
            lo, hi, med = (v.min(), v.max(), v.median()) if kind == 'range' else (v.p05, v.p95, v['median'])
            ax.plot([lo, hi], [y, y], color=c, lw=2.6, solid_capstyle='round')
            ax.plot(med, y, 'o', ms=3.2, mfc='white', mec=c, mew=0.9, zorder=3)
    ax.set_xscale('log'); ax.set_xlim(1e-3, 30)
    ax.set_yticks(range(n)); ax.set_yticklabels([x[0] for x in items][::-1], fontsize=6.8, color=INK); ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel('Heatwave days (days yr$^{-1}$)', fontsize=7.2, color=INK2)
    ax.legend(handles=[Line2D([], [], color='#0072B2', lw=2.6, label='Sichuan'), Line2D([], [], color='#D55E00', lw=2.6, label='Chongqing')],
              frameon=False, fontsize=6.6, loc='lower right')
    tag(ax, 'a', x=0.0)
    # (b) 尺度与物种，四川 2035
    ax = fig.add_axes([0.74, 0.6, 0.24, 0.36]); style(ax, 'x')
    dn = m[(m.region == '四川') & (m.year == 2035) & (m.metric == 'A_NPI_NPLS')].iloc[0]
    sel = [('This study: NPLS − NPI', -dn['median'], -dn.p95, -dn.p05),
           ('SC–CQ share: CO₂', 'own_share__OTNZ__co2only'),
           ('SC–CQ share:\n+ pollutants', 'own_share__OTNZ__all'),
           ('SC–CQ share:\n+ pollutants + CH₄', 'own_share__OTNZ__allCH4_V1'),
           ('China: CO₂', 'national__OTNZ__co2only'),
           ('China:\n+ pollutants + CH₄', 'national__OTNZ__allCH4_V1'),
           ('World: CO₂', 'universal__OTNZ__co2only'),
           ('World:\n+ pollutants + CH₄', 'universal__OTNZ__allCH4_V1')]
    n = len(sel)
    for i, s in enumerate(sel):
        if len(s) == 2:
            row = b1[(b1.case == s[1]) & (b1.year == 2035)].iloc[0]
            med, lo, hi = row.dTmax_p50, row.dTmax_p05, row.dTmax_p95
        else:
            med, lo, hi = s[1], s[2], s[3]
        y = n - 1 - i; c = COOL if med < 0 else WARM
        ax.plot([lo, hi], [y, y], color=c, lw=2.6, solid_capstyle='round', alpha=0.9)
        ax.plot(med, y, 'o', ms=3.2, mfc='white', mec=c, mew=0.9, zorder=3)
    ax.set_xscale('symlog', linthresh=1e-4, linscale=0.4); ax.set_xlim(-0.2, 0.7)
    ax.xaxis.set_minor_locator(plt.NullLocator()); ax.set_xticks([-0.1, -0.001, 0, 0.001, 0.1]); ax.set_xticklabels(['−0.1', '−0.001', '0', '0.001', '0.1'], fontsize=6.3)
    ax.axvline(0, color=AXIS, lw=0.6)
    ax.set_yticks(range(n)); ax.set_yticklabels([s[0] for s in sel][::-1], fontsize=6.8, color=INK); ax.tick_params(axis='y', length=0)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel('Sichuan Tmax change (°C)', fontsize=7.2, color=INK2)
    ax.legend(handles=[Line2D([], [], color=COOL, lw=2.6, label='median cooling'), Line2D([], [], color=WARM, lw=2.6, label='median warming')],
              frameon=False, fontsize=6.4, loc='upper right')
    tag(ax, 'b', x=0.0)
    # (c) 显现
    ax = fig.add_axes([0.085, 0.07, 0.75, 0.36]); style(ax)
    series = [('SC–CQ: CO$_2$', 'own_share__OTNZ__co2only', ':'), ('China: CO$_2$', 'national__OTNZ__co2only', '-.'),
              ('World: CO$_2$', 'universal__OTNZ__co2only', '-'), ('World: CO$_2$ +\npollutants + CH$_4$', 'universal__OTNZ__allCH4_V1', (0, (5, 2)))]
    top = 0
    for (lab, case, ls), c in zip(series, ['#009E73','#D55E00','#0072B2','#9B3A97']):
        d = b1[(b1.case == case) & (b1.year >= 2025)].sort_values('year')
        if case == 'universal__OTNZ__co2only':
            ax.fill_between(d.year, -d.dTmax_p95, -d.dTmax_p05, color=c, alpha=0.12, lw=0)
            top = max(top, (-d.dTmax_p05).max())
        ax.plot(d.year, -d.dTmax_p50, color=c, lw=1.9, ls=ls, label=lab)
        ax.annotate(lab, (2100,-d.dTmax_p50.iloc[-1]),xytext=(6,0),textcoords='offset points',ha='left',va='center',color=c,fontsize=8,annotation_clip=False)
    thr20, thr1 = 0.2087, 0.4666
    ax.axhline(thr20, color=INK2, lw=0.7, ls=(0, (3, 2))); ax.axhline(thr1, color=INK2, lw=0.7, ls=(0, (1, 1.5)))
    ax.text(2026, thr20 + 0.02, '20-year threshold (0.21 °C)', fontsize=6.3, color=INK2, ha='left', va='bottom')
    ax.text(2026, thr1 + 0.02, 'Single-year threshold (0.47 °C)', fontsize=6.3, color=INK2, ha='left', va='bottom')
    ax.axhline(0, color=AXIS, lw=0.6)
    ax.axvspan(2060.5, 2100, color='#f6f6f3', lw=0, zorder=0)
    ax.text(2061.5, -0.33, 'Post-2060 extension', fontsize=6, color=MUTED)
    ax.set_xlim(2025, 2100); ax.set_ylim(-0.4, top + 0.05)
    ax.set_ylabel('Sichuan Tmax cooling vs SSP2-4.5 (°C)', fontsize=7.5, color=INK2)
    ax.set_xticks([2025,2040,2060,2080,2100])
    tag(ax, 'c', x=0.0)
    save(fig, 'Figure5_heatwave_scale_emergence')


# ---------------------------------------------------------------- 图 S4：全球基准，2023 vs 2026 起点
def figS4():
    sup = pd.read_csv(DATA / 'fair_supplement/regional_timeseries_2023_2100.csv')
    thr = json.load(open(EXP / 'scales/analysis_inputs.json', encoding='utf-8'))['noise']
    fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.0), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.975, bottom=0.14, top=0.9, wspace=0.08)
    for ax, r in zip(axs, ['sichuan', 'chongqing']):
        style(ax)
        for case, c, ls, lab in [('ssp119__start2023', '#0072B2', '-', 'SSP1-1.9, from 2023'), ('ssp119__start2026', '#0072B2', (0, (4, 2)), 'SSP1-1.9, from 2026'),
                                 ('ssp126__start2023', '#D55E00', '-', 'SSP1-2.6, from 2023'), ('ssp126__start2026', '#D55E00', (0, (4, 2)), 'SSP1-2.6, from 2026')]:
            d = sup[(sup.region == r) & (sup.case == case) & (sup.year >= 2023)].sort_values('year')
            if ls == '-':
                ax.fill_between(d.year, -d.dTmax_p95, -d.dTmax_p05, color=c, alpha=0.08, lw=0)
            ax.plot(d.year, -d.dTmax_p50, color=c, lw=1.4, ls=ls, label=lab)
        s = thr[r]['t']
        ax.axhline(2 * s / np.sqrt(20), color=INK2, lw=0.7, ls=(0, (3, 2))); ax.axhline(s, color=INK2, lw=0.7, ls=(0, (1, 1.5)))
        ax.axhline(0, color=AXIS, lw=0.6)
        ax.set_xlim(2023, 2100); ax.set_ylim(-0.3, 2.2)
        ax.set_title(('A  ' if r == 'sichuan' else 'B  ') + REGN[r], fontsize=8, color=INK, loc='left')
    axs[0].set_ylabel('Tmax cooling vs SSP2-4.5 (°C)', fontsize=7.5, color=INK2)
    axs[0].text(2098, 2 * thr['sichuan']['t'] / np.sqrt(20) + 0.03, '2σ/√20', fontsize=6.3, color=INK2, ha='right')
    axs[0].text(2098, thr['sichuan']['t'] + 0.03, '1σ', fontsize=6.3, color=INK2, ha='right')
    axs[1].legend(frameon=False, fontsize=6.6, loc='upper left')
    save(fig, 'FigureS4_global_benchmarks_start_year')


if __name__ == '__main__':
    import sys
    todo = sys.argv[1:] or ['fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'figS4']
    for f in todo:
        globals()[f]()
    print('saved to', FIG)
