"""热浪实验共用的数值定义（带版本）。"""
from pathlib import Path
import os,json
from datetime import datetime,timezone
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'2')
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parents[1]
ROOT=HERE.parent
REGIONS=['sichuan','chongqing']
DEFINITIONS=['annual90','annual95','warm90']
KINDS=['day','night','compound']
SEED=20260921

def save_json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=str),encoding='utf-8')

def grid_meta():
    df=pd.read_csv(ROOT/'inputs/weather/grid_metadata_population.csv')
    assert len(df)==5326 and (df.population2020>=0).all()
    return df


def load_cache():
    meta=grid_meta();dates=pd.date_range('1990-01-01','2024-12-31')
    cache=ROOT/'inputs/weather'
    tx=np.load(cache/'tmax_float64.npy',mmap_mode='r');tn=np.load(cache/'tmin_float64.npy',mmap_mode='r')
    assert tx.shape==tn.shape==(len(dates),len(meta))
    return dates,meta,tx,tn

def design(dates,trend=True):
    d=pd.DatetimeIndex(dates);phase=(d.dayofyear.to_numpy()-1)/np.where(d.is_leap_year,366.,365.)
    terms=[np.ones(len(d))]
    if trend:terms.append(d.year.to_numpy()+phase-1990)
    terms.extend(v for k in (1,2,3) for v in (np.sin(2*np.pi*k*phase),np.cos(2*np.pi*k*phase)))
    return np.column_stack(terms)

def represent(tx,tn):
    assert np.isfinite(tx).all() and np.isfinite(tn).all() and np.all(tx>tn)
    return np.stack(((tx+tn)*.5,np.log(tx-tn)),axis=-1)

def restore(values):
    half=np.exp(values[...,1])*.5
    tx=values[...,0]+half;tn=values[...,0]-half
    assert np.isfinite(tx).all() and np.isfinite(tn).all() and np.all(tn<=tx)
    return tx,tn

def mark(hot):
    """先按时间排序；符合条件的连续段可能跨年。"""
    hot=np.asarray(hot,dtype=bool);out=np.zeros_like(hot)
    triples=hot[:-2]&hot[1:-1]&hot[2:]
    for i in range(3):out[i:i+len(triples)]|=triples
    return out

def annual_metrics(hot,dates):
    d=pd.DatetimeIndex(dates)
    assert len(d)==len(hot) and np.all(np.diff(d.values)==np.timedelta64(1,'D'))
    marked=mark(hot);start=marked.copy();start[1:]&=~marked[:-1]
    years,index=np.unique(d.year,return_index=True)
    return years,np.add.reduceat(marked,index,axis=0,dtype=np.int32),np.add.reduceat(start,index,axis=0,dtype=np.int32)

def thresholds(tx,tn,dates,definition='annual90'):
    d=pd.DatetimeIndex(dates);train=np.asarray((d.year>=1990)&(d.year<=2005))
    if definition=='warm90':train&=np.asarray((d.month>=5)&(d.month<=9))
    q=.95 if definition=='annual95' else .9
    return np.quantile(tx[train],q,axis=0),np.quantile(tn[train],q,axis=0)

def donor_plan(train_dates,target_dates,candidate,members=500,seed=SEED):
    """所有格点和两个温度分量共用的供体索引。

    数组形状为 目标日 × 成员；供体年没有 2 月 29 日时在 M/logDTR 空间插值；
    30 天块要求供体有完整 30 天。
    """
    rd=pd.DatetimeIndex(train_dates);td=pd.DatetimeIndex(target_dates)
    rng=np.random.default_rng(seed)
    a=np.empty((len(td),members),dtype=np.int32);b=np.empty_like(a)
    w=np.zeros(a.shape,dtype=np.float32)
    if candidate=='block30':
        rp=(rd.dayofyear.to_numpy()-1)*365.25/np.where(rd.is_leap_year,366.,365.)
        tp=(td.dayofyear.to_numpy()-1)*365.25/np.where(td.is_leap_year,366.,365.)
        valid=np.arange(len(rd)-29)
        for offset in range(0,len(td),30):
            length=min(30,len(td)-offset);dist=np.abs(rp[valid]-tp[offset]);dist=np.minimum(dist,365.25-dist)
            chosen=rng.choice(valid[dist<=15],size=members)
            a[offset:offset+length]=chosen[None,:]+np.arange(length)[:,None]
        b[:]=a
    elif candidate in ('whole_year','seasonal_climatology'):
        donors=np.unique(rd.year)
        lookup={int(yy):{int(rd[i].month*100+rd[i].day):int(i) for i in np.flatnonzero(rd.year==yy)} for yy in donors}
        for yy in np.unique(td.year):
            pos=np.flatnonzero(td.year==yy);chosen=rng.choice(donors,size=members)
            for donor in np.unique(chosen):
                mm=np.flatnonzero(chosen==donor);la=[];lb=[];lw=[]
                for date in td[pos]:
                    md=date.month*100+date.day
                    if md in lookup[donor]:la.append(lookup[donor][md]);lb.append(lookup[donor][md]);lw.append(0.)
                    else:
                        assert md==229
                        la.append(lookup[donor][228]);lb.append(lookup[donor][301]);lw.append(.5)
                a[np.ix_(pos,mm)]=np.array(la)[:,None];b[np.ix_(pos,mm)]=np.array(lb)[:,None];w[np.ix_(pos,mm)]=np.array(lw)[:,None]
    else:raise ValueError(candidate)
    assert a.min()>=0 and b.max()<len(rd)
    return a,b,w

def draw_residuals(residuals,plan,member_slice=slice(None)):
    a,b,w=(z[:,member_slice] for z in plan)
    ww=w.reshape(w.shape+(1,)*(residuals.ndim-1))
    return residuals[a]+ww*(residuals[b]-residuals[a])

def crps(samples,obs,axis=-1):
    s=np.moveaxis(np.asarray(samples,dtype=float),axis,-1);n=s.shape[-1]
    ranks=2*np.arange(1,n+1)-n-1
    return np.mean(abs(s-np.asarray(obs)[...,None]),axis=-1)-np.sum(np.sort(s,axis=-1)*ranks,axis=-1)/(n*n)

def summary(values):
    a=np.asarray(values,dtype=float)
    return dict(mean=float(a.mean()),median=float(np.median(a)),q025=float(np.quantile(a,.025)),q975=float(np.quantile(a,.975)))
