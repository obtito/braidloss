"""Controlled parameter figures, created only from saved evaluated records."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from optimize_design import design_id

ROOT=Path(__file__).resolve().parents[1]

def plot():
    records={p.stem:json.loads(p.read_text(encoding='utf-8')) for p in (ROOT/'data/search_64').glob('*.json')}
    if len(records)!=129:raise RuntimeError('The predeclared search is not complete')
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,
        'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(14,4.8),facecolor='#f7fafc')
    fig.subplots_adjust(top=.73,bottom=.2,left=.06,right=.97,wspace=.31)
    pitch=[16,32,64]
    reference_K=records[design_id((4,4,4),(32,32,96))]['ROM']['K']
    for level,label,color in [(0,'只改变内级节距','#147d92'),(1,'只改变中级节距','#b66a29')]:
        y=[]
        for value in pitch:
            p=[32,32,96];p[level]=value
            r=records[design_id((4,4,4),p)]
            y.append(100*(r.get('ROM',{}).get('K',np.nan)/reference_K-1))
        axes[0].plot(pitch,y,'o-',color=color,label=label)
    axes[0].set(xlabel='改变的节距 / mm',ylabel='相对参考 K 的变化 / %',ylim=(-1,1),title='① 单因素节距对照')
    axes[0].axhline(0,color='#8696a7',lw=.8,ls='--')
    axes[0].set_xticks(pitch);axes[0].legend(fontsize=9);axes[0].grid(alpha=.15)
    labels=[];y=[]
    for a,b in [(1,1),(-1,1),(1,-1),(-1,-1)]:
        r=records[design_id((4,4,4),(32*a,32*b,96))]
        labels.append(('+' if a>0 else '−')+'/'+('+' if b>0 else '−'))
        y.append(r.get('ROM',{}).get('K',np.nan))
    axes[1].bar(labels,y,color='#577fa1',width=.6)
    axes[1].set(xlabel='内级 / 中级绞向；外级固定 +',ylabel='Rac/Rdc',title='② 仅改变绞向')
    for i,v in enumerate(y):
        if np.isfinite(v):axes[1].text(i,v,f'{v:.3f}',ha='center',va='bottom',fontsize=9)
    axes[1].grid(axis='y',alpha=.15)
    for levels,label,color in [(1,'单级','#b66a29'),(2,'两级','#147d92'),(3,'三级','#6957a4')]:
        selected=[r for r in records.values() if len(r['factors'])==levels and 'ROM' in r]
        axes[2].scatter([100*(r['geometry']['growth_mean']-1) for r in selected],
                         [r['ROM']['K'] for r in selected],s=23,alpha=.65,label=label,color=color)
    axes[2].set(xlabel='平均实际丝长增加 / %',ylabel='Rac/Rdc',title='③ 同时报告几何代价')
    axes[2].legend(fontsize=9);axes[2].grid(alpha=.15)
    fig.text(.06,.925,'Python 参数搜索：受控比较，保留全部候选',fontsize=21,weight='bold',color='#142a43')
    fig.text(.06,.84,'固定 64 股 · 6 mm² · 200 kHz · 20 A；左、中图采用 4×4×4，参考节距 32 / 32 / 96 mm。',fontsize=11,color='#526478')
    fig.text(.06,.062,'Python 近似预测；最终方案以 COMSOL 复算为准。左图不足 1% 的差异不能据此认定优劣。缺失点表示未满足预设约束。',fontsize=9,color='#526478')
    out=ROOT/'figures/parameter_controls.png';fig.savefig(out,dpi=180,facecolor=fig.get_facecolor());plt.close(fig)
    return out

if __name__=='__main__':print(plot())
