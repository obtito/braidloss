"""Explain the recommended recursive stranding with its actual trajectories."""
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
    c=RecursiveCable((4,16),(-16.,128.),gap_um=40.)
    cross=tube_surface(c,[0.],ntheta=40)[:,0]*1000
    xyz,_,_=c.at_physical_z(np.linspace(0,.032,321));xyz*=1000
    colors=plt.get_cmap('tab20')(np.arange(16)/15)
    strandcolors=['#167d9a','#cf7930','#75569e','#36977b']
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,
        'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(13,8.7),facecolor='#f7fafc')
    gs=fig.add_gridspec(2,2,left=.07,right=.95,bottom=.14,top=.8,wspace=.34,hspace=.49)
    ax=fig.add_subplot(gs[0,0]);ax.set_facecolor('white')
    mid=cross[:4].mean(axis=(0,1))
    for j in range(4):
        ax.add_patch(Polygon(cross[j,:,:2]-mid[:2],facecolor=strandcolors[j],edgecolor='white',lw=.8))
        ax.text(*(cross[j,:,:2].mean(0)-mid[:2]),str(j+1),ha='center',va='center',color='white',weight='bold')
    ax.set_aspect('equal');ax.set(xlim=(-.65,.65),ylim=(-.55,.55),xlabel='x / mm',ylabel='y / mm',title='① 4 根绝缘圆铜丝，绞成一个子束')
    ax.text(.5,-.36,'单丝直径 0.3455 mm；内级相对节距 −16 mm',transform=ax.transAxes,ha='center',color='#526478',fontsize=9)
    ax=fig.add_subplot(gs[0,1]);ax.set_facecolor('white')
    for j in range(64):ax.add_patch(Polygon(cross[j,:,:2],facecolor=colors[j//4],edgecolor='white',lw=.45))
    ax.set_aspect('equal');ax.set(xlim=(-2.5,2.5),ylim=(-2.5,2.5),xlabel='x / mm',ylabel='y / mm',title='② 16 个子束，绞成 64 股总束')
    ax.text(.5,-.36,'同色属于同一子束；外级相对节距 +128 mm',transform=ax.transAxes,ha='center',color='#526478',fontsize=9)
    ax=fig.add_subplot(gs[1,0],projection='3d');ax.set_facecolor('#f7fafc')
    for j in range(64):ax.plot(xyz[j,:,2],xyz[j,:,0],xyz[j,:,1],lw=.45,color=colors[j//4],alpha=.13)
    selected=[24,25,26,27]
    for j,color in zip(selected,strandcolors):
        ax.plot(xyz[j,:,2],xyz[j,:,0],xyz[j,:,1],color=color,lw=1.9)
    ax.set(xlabel='轴向 z / mm',ylabel='x / mm',zlabel='y / mm',title='③ 子束内旋转，同时绕总束轴线旋转')
    ax.set_box_aspect((2.8,1,1));ax.view_init(22,-68);ax.tick_params(labelsize=8)
    ax=fig.add_subplot(gs[1,1]);ax.set_facecolor('white')
    for j,color in zip(selected,strandcolors):
        ax.plot(xyz[j,:,2],np.linalg.norm(xyz[j,:,:2],axis=1),color=color,lw=1.5,label=f'股 {j+1}')
    ax.set(xlabel='轴向 z / mm',ylabel='距总轴线 / mm',title='④ 同一根丝的位置随长度周期变化')
    ax.grid(alpha=.15);ax.legend(ncol=4,fontsize=8,loc='upper right')
    fig.text(.07,.945,'推荐结构怎样实现：4 股 → 子束 → 64 股',fontsize=22,weight='bold',color='#142a43')
    fig.text(.07,.892,'内、外级采用相反绞向。子丝位于父束轴线的法向平面内，各级轨迹递归生成，再扫掠成圆铜实体。',fontsize=11,color='#526478')
    fig.text(.07,.054,'这张图解释几何，颜色用于区分股线；电流密度另见实际场图。多级绞合产生径向移动，但不保证每根丝遍历整个截面。',fontsize=9,color='#526478')
    fig.text(.07,.024,'节距按全局轴向参数在父束材料坐标系中相对转一圈定义；实际局部绞角和丝长由轨迹计算。三维图调整了轴向显示比例。',fontsize=8.6,color='#526478')
    out=ROOT/'figures/recommended_structure.png'
    fig.savefig(out,dpi=170,facecolor=fig.get_facecolor());plt.close(fig)
    return out


if __name__=='__main__':print(draw())
