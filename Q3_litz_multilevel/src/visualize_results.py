"""Scientific figures from saved COMSOL field samples and numerical evidence."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from summarize_pilots import load_result

ROOT=Path(__file__).resolve().parents[1]


def field_dashboard(case,final=False):
    folder=ROOT/'data/raw'/case
    cfg=json.loads((folder/'inputs.json').read_text(encoding='utf-8'))
    sample=np.load(folder/'sampling_points.npz')
    surf=sample['surface_m'];cross=sample['cross_m'];coords=sample['points_m']
    raw=(folder/'stdout.log').read_text(encoding='utf-8',errors='replace')
    values=np.array([float(l.split('|')[-1]) for l in raw.splitlines() if l.startswith('FIELD|')])
    if len(values)!=len(coords) or not np.isfinite(values).all():raise ValueError('Incomplete field values')
    result=load_result(folder);m=result['METRIC'];current=result['STRAND'][::3,-1]
    u=values[np.prod(surf.shape[:-1]):].reshape(cross.shape[:-1])/1e6
    n=len(cross);total=abs(current.sum());reference=total/n
    pos=cross[:,0,0,:2]*1000;radius=np.linalg.norm(pos,axis=1)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,
                         'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(13,7),facecolor='#f7fafc')
    ax=fig.add_axes([.055,.15,.47,.66],facecolor='white')
    norm=Normalize(0,35);cmap=plt.get_cmap('turbo')
    for j in range(n):
        xy=np.vstack([cross[j,0,:1,:2],cross[j,1:,:,:2].reshape(-1,2)])*1000
        val=np.r_[u[j,0,0],u[j,1:].ravel()]
        tri=mtri.Triangulation(xy[:,0],xy[:,1])
        ax.tripcolor(tri,val,shading='gouraud',cmap=cmap,norm=norm,rasterized=True)
    ax.set_aspect('equal');ax.set(xlabel='x / mm',ylabel='y / mm',title='铜截面内的电流密度有效值')
    lim=np.max(abs(cross[:,:,:,:2]))*1000+.2;ax.set_xlim(-lim,lim);ax.set_ylim(-lim,lim)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=ax,fraction=.045,pad=.035)
    cb.set_label('|J| RMS / (A/mm²)')
    bx=fig.add_axes([.65,.47,.30,.31],facecolor='white')
    bx.scatter(radius,abs(current),s=25,color='#147d92',alpha=.8,edgecolors='none')
    bx.axhline(reference,color='#c05730',ls='--',lw=1.2,label='总电流 / 股数，仅作参考')
    bx.set(xlabel='股线中心距总轴线 / mm',ylabel='各股电流幅值 RMS / A',title='各股电流由场方程求解')
    bx.legend(fontsize=8);bx.grid(alpha=.15)
    fig.text(.65,.37,'按轴向每米、20 A RMS 评估',fontsize=13,weight='bold',color='#142a43')
    lines=[f'直流极限电阻：{m[4,0].real*1000:.4f} mΩ/m',
           f'交流电阻：{m[4,-1].real*1000:.4f} mΩ/m',
           f'Rac/Rdc = {m[4,-1].real/m[4,0].real:.4f}',
           f'铜损：{m[6,-1].real:.4f} W/m']
    for k,line in enumerate(lines):fig.text(.65,.315-.046*k,line,color='#334e68',fontsize=11)
    status='通过指定验证的模型' if final else '初算结果，正式验收仍在进行'
    fig.text(.055,.935,f'6 mm² Litz：{n} 股 · 200 kHz · 20 A',fontsize=22,weight='bold',color='#142a43')
    fig.text(.055,.88,f'真实三维 COMSOL 场的截面采样 | {status}',fontsize=11,color='#526478')
    fig.text(.055,.07,'颜色是局部电流密度；右图是每股净电流。高频下，同一根丝内部也会出现不均匀分布。',fontsize=10,color='#526478')
    fig.text(.055,.033,f'数据：{case}；取单元中间截面；绘图采样至铜半径的 98%。各股未预设为 I/N。',fontsize=8.5,color='#526478')
    out=ROOT/'figures'/(case+'_field.png');out.parent.mkdir(exist_ok=True)
    fig.savefig(out,dpi=170,facecolor=fig.get_facecolor());plt.close(fig)
    return out


if __name__=='__main__':
    import sys
    print(field_dashboard(sys.argv[1]))
