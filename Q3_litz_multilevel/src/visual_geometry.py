"""Scientific geometry views and points for actual COMSOL field sampling."""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline
from recursive_geometry import RecursiveCable

ROOT=Path(__file__).resolve().parents[1]


def tube_surface(cable, zlevels, ntheta=24, radius_fraction=1.):
    """Exact normal-circle tube evaluated at common physical-z cut planes."""
    phi=np.arange(ntheta)*2*np.pi/ntheta
    zgrid,ang=np.meshgrid(np.asarray(zlevels),phi,indexing='ij')
    z=zgrid.ravel();ang=ang.ravel()
    a=cable.a*radius_fraction
    result=[]
    for j in range(cable.n):
        cs=CubicSpline(cable.t,cable.xyz[j]-cable.t[:,None]*[0,0,1],bc_type='periodic',extrapolate='periodic')
        t=z.copy()
        for _ in range(8):
            c=cs(t);c[:,2]+=t
            d=cs(t,1)+[0,0,1];dd=cs(t,2)
            speed=np.linalg.norm(d,axis=1)
            tangent=d/speed[:,None]
            dtan=(dd-tangent*np.sum(tangent*dd,axis=1)[:,None])/speed[:,None]
            u=np.cross([0,1,0],tangent);du=np.cross([0,1,0],dtan)
            un=np.linalg.norm(u,axis=1);e1=u/un[:,None]
            de1=(du-e1*np.sum(e1*du,axis=1)[:,None])/un[:,None]
            e2=np.cross(tangent,e1);de2=np.cross(dtan,e1)+np.cross(tangent,de1)
            q=c+a*(np.cos(ang)[:,None]*e1+np.sin(ang)[:,None]*e2)
            dq=d+a*(np.cos(ang)[:,None]*de1+np.sin(ang)[:,None]*de2)
            step=(q[:,2]-z)/dq[:,2]
            t-=step
            if abs(step).max()<1e-13:break
        if abs(q[:,2]-z).max()>1e-10:raise ValueError('Tube cross-section inverse did not converge')
        result.append(q.reshape(len(zlevels),ntheta,3))
    return np.array(result)


def overview():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei'],'axes.unicode_minus':False,'font.size':11,
                         'axes.spines.top':False,'axes.spines.right':False})
    c=RecursiveCable((7,7,7),(10.,30.,90.))
    xyz,tan,_=c.at_physical_z(np.linspace(0,c.period,361))
    colors=plt.get_cmap('tab10')(np.arange(7))
    fig=plt.figure(figsize=(13.2,8.5),facecolor='#f7fafc')
    gs=fig.add_gridspec(2,2,left=.065,right=.98,top=.86,bottom=.09,hspace=.35,wspace=.28)
    ax=fig.add_subplot(gs[0,0]);ax.set_facecolor('#fff')
    for j in range(c.n):
        color=colors[j//49]
        ax.add_patch(Circle(xyz[j,0,:2]*1000,c.a*1000,facecolor=color,edgecolor='white',lw=.25,alpha=.85))
    ax.set(xlim=(-2.6,2.6),ylim=(-2.6,2.6),xlabel='x / mm',ylabel='y / mm',title='① 分组：7 股 → 49 股子束 → 343 股总束')
    ax.set_aspect('equal');ax.grid(alpha=.1)
    ax.text(.01,.01,'颜色区分最外一级的 7 个子束；圆符号标记股线',transform=ax.transAxes,fontsize=9)
    a3=fig.add_subplot(gs[0,1],projection='3d');a3.set_facecolor('#f7fafc')
    keep=xyz[:, :121]*1000
    for j in range(c.n):
        a3.plot(keep[j,:,2],keep[j,:,0],keep[j,:,1],color=colors[j//49],lw=.4,alpha=.18)
    for j in [0,1,7,50,57,99,342]:
        a3.plot(keep[j,:,2],keep[j,:,0],keep[j,:,1],color=colors[j//49],lw=1.5,alpha=1)
    a3.set(xlabel='轴向 z / mm',ylabel='x / mm',zlabel='y / mm',title='② 真实递归中心线，显示前 30 mm')
    a3.set_box_aspect((2.8,1,1));a3.view_init(20,-70)
    ax=fig.add_subplot(gs[1,0]);ax.set_facecolor('#fff')
    for j,label in [(0,'总束中心股'),(1,'中心子束内的外圈股'),(7,'中心子束内的次级中心股'),(57,'外层子束中的股线')]:
        ax.plot(xyz[j,:,2]*1000,np.linalg.norm(xyz[j,:,:2],axis=1)*1000,label=label,lw=1.6)
    ax.set(xlabel='轴向 z / mm',ylabel='距总束轴线的半径 / mm',title='③ 有径向移动，但并非每股都遍历内外层')
    ax.legend(fontsize=9,ncol=1);ax.grid(alpha=.17)
    ax=fig.add_subplot(gs[1,1]);ax.axis('off')
    ax.text(0,1,'④ 当前候选参数与验证范围',fontsize=14,weight='bold',va='top',color='#142a43')
    items=[('固定工况','6 mm² · 200 kHz · 20 A RMS · 20°C'),
           ('几何候选','7×7×7；丝径 149.24 μm'),
           ('相对节距','内→外：10 / 30 / 90 mm'),
           ('几何代价','平均丝长 +2.32%；外径上界约 4.56 mm'),
           ('完整重复','90 mm；5 mm + 旋转 20° 的几何映射已检查'),
           ('当前证据','单丝与 3×3 小规模实体场试算；正式基准未冻结')]
    for k,(key,val) in enumerate(items):
        y=.8-k*.128
        ax.text(0,y,key,weight='bold',fontsize=11,color='#53677b')
        ax.text(.24,y,val,fontsize=10.3,color='#142a43')
    fig.text(.065,.955,'传统多级绞合 Litz：从结构到验证',fontsize=22,weight='bold',color='#142a43')
    fig.text(.065,.905,'本图展示几何与研究进度。它不是 343 股已求解的三维电磁场，也不是最优参数结果。',fontsize=11,color='#526478')
    fig.text(.065,.025,'“各级节距”相对于父束材料坐标系定义；绞角和实际丝长由递归轨迹推导。三维图为便于观察调整了坐标轴显示比例。',fontsize=9,color='#526478')
    folder=ROOT/'figures';folder.mkdir(exist_ok=True)
    fig.savefig(folder/'recursive_geometry_overview.png',dpi=180,facecolor=fig.get_facecolor())
    plt.close(fig)
    return folder/'recursive_geometry_overview.png'


if __name__=='__main__':print(overview())
