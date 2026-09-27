"""Controlled structural comparison; drawings are geometry, not field results."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from recursive_geometry import RecursiveCable
from visual_geometry import tube_surface

ROOT=Path(__file__).resolve().parents[1]

def draw():
    cases=[('单级：64 股', (64,), (96.,)),
           ('两级：4 × 16', (4,16), (16.,96.)),
           ('三级：4 × 4 × 4', (4,4,4), (32.,32.,96.))]
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,
        'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(13.8,10),facecolor='#f7fafc')
    gs=fig.add_gridspec(3,3,left=.06,right=.97,top=.86,bottom=.14,
                        height_ratios=[1.1,1.05,.9],wspace=.32,hspace=.47)
    for column,(name,factors,pitches) in enumerate(cases):
        c=RecursiveCable(factors,pitches,gap_um=40.)
        z=np.linspace(0,.096,385);xyz,_,_=c.at_physical_z(z)
        surface=tube_surface(c,[0.],ntheta=32)[:,0,:,:2]*1000
        group=np.arange(64)//int(np.prod(factors[:-1])) if len(factors)>1 else np.arange(64)//8
        palette=plt.get_cmap('tab20')(group/ max(group.max(),1))
        ax=fig.add_subplot(gs[0,column]);ax.set_facecolor('white')
        for j in range(64):ax.add_patch(Polygon(surface[j],facecolor=palette[j],edgecolor='white',lw=.35))
        ax.set_aspect('equal');ax.set(xlim=(-3,3),ylim=(-3,3),xlabel='x / mm',ylabel='y / mm',title=name+'：起始截面')
        ax.grid(alpha=.12)
        bx=fig.add_subplot(gs[1,column],projection='3d');bx.set_facecolor('#f7fafc')
        p=xyz[:,:129]*1000
        for j in range(64):bx.plot(p[j,:,2],p[j,:,0],p[j,:,1],color=palette[j],lw=.5,alpha=.2)
        tracked=[0,17,63]
        colors=['#117a91','#b76a27','#69459b']
        for j,color in zip(tracked,colors):bx.plot(p[j,:,2],p[j,:,0],p[j,:,1],color=color,lw=2.)
        bx.set(xlabel='z / mm',ylabel='x / mm',zlabel='y / mm',title='前 32 mm 的股线中心线')
        bx.set_box_aspect((2.5,1,1));bx.view_init(22,-68);bx.tick_params(labelsize=8)
        cx=fig.add_subplot(gs[2,column]);cx.set_facecolor('white')
        for j,color in zip(tracked,colors):
            cx.plot(z*1000,np.linalg.norm(xyz[j,:,:2],axis=1)*1000,color=color,lw=1.4,label=f'第 {j+1} 股')
        cx.set(xlabel='轴向 z / mm',ylabel='距总轴线 / mm',ylim=(0,3),title='追踪同一根丝的位置变化')
        cx.legend(fontsize=8,loc='upper right',ncol=3);cx.grid(alpha=.13)
        label=' / '.join(f'{x:g}' for x in pitches)
        cx.text(.5,-.38,f'相对节距（内→外）：{label} mm',transform=cx.transAxes,ha='center',fontsize=9,color='#526478')
    fig.text(.06,.95,'同样 64 股，改变逐级绞合方式',fontsize=22,weight='bold',color='#142a43')
    fig.text(.06,.9,'固定铜面积 6 mm²、单丝直径 0.3455 mm、名义间距 40 μm；多级绞合使单丝相对总轴线产生径向移动。',fontsize=11,color='#526478')
    fig.text(.06,.028,'这是几何示意，不是电流云图。截面使用同一尺度；三维中心线图缩短了轴向显示比例。分组变化也会改变外径和空隙率。',fontsize=9,color='#526478')
    out=ROOT/'figures/controlled_structure_geometry.png'
    fig.savefig(out,dpi=170,facecolor=fig.get_facecolor());plt.close(fig)
    return out

if __name__=='__main__':print(draw())
