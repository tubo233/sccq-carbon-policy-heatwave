"""用现有结果画图 S4 的尺度和物种延续面板。"""
from pathlib import Path
import csv,json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parents[1]
with (HERE/'data/fair_supplement/regional_timeseries_2023_2100.csv').open(encoding='utf-8') as f: times=list(csv.DictReader(f))
noise=json.loads((HERE/'models/scales/analysis_inputs.json').read_text(encoding='utf-8'))['noise']
# 补图 S4 时间序列的配套面板。
plt.rcParams.update({'font.family':'Arial', 'font.size':9, 'pdf.fonttype':42, 'axes.spines.top':False, 'axes.spines.right':False})
fig, axs = plt.subplots(1,2,figsize=(7.1,4.4),sharey=True)
colors = {'own_share':'#c75d21','national':'#258552','universal':'#175a9e'}
labels = {'own_share':'SC–CQ share','national':'China','universal':'World'}
for ax, region, panel in zip(axs,['sichuan','chongqing'],['C','D']):
    for scale in colors:
        for variant, ls in [('co2only','-'),('all','--'),('allCH4_V1',':')]:
            rr = sorted([r for r in times if r['case']==f'{scale}__OTNZ__{variant}' and r['region']==region], key=lambda r:int(r['year']))
            ax.plot([int(r['year']) for r in rr],[-float(r['dTmax_p50']) for r in rr], color=colors[scale],ls=ls,lw=1.5)
    ax.axhline(0,color='#777777',lw=.6)
    ax.axhline(noise[region]['t'],color='#777777',ls=(0,(1,2)),lw=.8)
    ax.axhline(2*noise[region]['t']/np.sqrt(20),color='#777777',ls=(0,(4,3)),lw=.8)
    ax.set_title(f'{panel}  {region.title()}',loc='left',fontsize=11)
    ax.set_xlim(2023,2100); ax.set_xlabel('Year')
    ax.grid(axis='y',alpha=.2)
    ax.set_yscale('symlog',linthresh=.01,linscale=.7)
    ax.set_yticks([-.3,-.1,-.01,0,.01,.1,.3,1,2])
    ax.set_yticklabels(['−0.3','−0.1','−0.01','0','0.01','0.1','0.3','1','2'])
axs[0].set_ylabel('Tmax cooling relative to SSP2-4.5 (°C; symlog)')
from matplotlib.lines import Line2D
handles = [Line2D([],[],color=colors[s],lw=1.5,label=labels[s]) for s in colors]
handles += [Line2D([],[],color='#333333',lw=1.5,ls=ls,label=lab) for ls,lab in [('-','CO$_2$ only'),('--','DPEC species'),(':','DPEC species + CH$_4$ (V1)')]]
fig.legend(handles=handles,ncol=3,loc='lower center',frameon=False,bbox_to_anchor=(.52,-.01),fontsize=8)
fig.subplots_adjust(left=.09,right=.98,top=.91,bottom=.25,wspace=.08)
fig.savefig(HERE/'figures/FigureS4_scale_species_timeseries.png',dpi=300,facecolor='white')
fig.savefig(HERE/'figures/FigureS4_scale_species_timeseries.pdf',facecolor='white')
plt.close(fig)

