"""Match XY sections, full 3D path and a finite-radius local piece by strand ID."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

out=Path(__file__).resolve().parent
actual=json.loads((out/'actual_mph_strands.json').read_text())
data=json.loads((out.parent/'render_data.json').read_text())
a=np.array(actual['43']['points_m'])*1000
b=np.array(actual['21']['points_m'])*1000
c=np.array(data['coeff484']);q=data['parameters']['station_advance_mm']
d=data['parameters']['diameter_um']/1000
grid=c[:,0,:]
blue='#2478b4';orange='#d6812e';ink='#243647';grey='#b9c1ca'
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei'],'axes.unicode_minus':False,'font.size':12})

fig=plt.figure(figsize=(15,11),dpi=180,facecolor='white')
gs=GridSpec(2,5,figure=fig,height_ratios=[1,2.4],top=.86,bottom=.14,left=.035,right=.95,hspace=.28,wspace=.33)
fig.suptitle('二维截面的蓝点，沿 z 方向连起来，就是三维中的同一根线',fontsize=22,y=.98,color=ink)
fig.text(.5,.938,'统一追踪 Strand 44（path43）｜原始 COMSOL .mph 坐标重绘，非软件截图',ha='center',fontsize=14,color='#536273')
fig.text(.5,.897,'① 从同一根线的完整 500 mm 周期中，取出五张横截面',ha='center',fontsize=15,color=ink)
indices=[0,242,484,726,968]
for i,k in enumerate(indices):
    ax=fig.add_subplot(gs[0,i]);pt=a[k]
    ax.scatter(grid[:,0],grid[:,1],s=5,color=grey,linewidth=0)
    ax.scatter([pt[0]],[pt[1]],s=210,color=blue,edgecolor='white',linewidth=1.5,zorder=5)
    ax.text(pt[0],pt[1],'A',ha='center',va='center',color='white',weight='bold',fontsize=12,zorder=6)
    ax.set_title(f'S{i}  |  z = {pt[2]:.0f} mm',fontsize=13,pad=7,color=ink)
    ax.set_xlim(-2.65,2.65);ax.set_ylim(-2.65,2.65);ax.set_aspect('equal')
    ax.set_xticks([-2,0,2]);ax.set_yticks([-2,0,2]);ax.tick_params(labelsize=10,length=2)
    ax.set_xlabel('x / mm',fontsize=10);ax.set_ylabel('y / mm',fontsize=10)
    for spine in ax.spines.values():spine.set_color('#d6dce3')

ax=fig.add_subplot(gs[1,:3],projection='3d')
ax.set_title('② 保留每个蓝点的 z 坐标，连成完整三维中心线',fontsize=15,pad=16,color=ink)
ax.plot(a[:,0],a[:,1],a[:,2],color=blue,lw=1.9)
ax.plot(a[:9,0],a[:9,1],a[:9,2],color=orange,lw=4.8)
for i,k in enumerate(indices):
    z=a[k,2];corners=[[-2.48,-2.48,z],[2.48,-2.48,z],[2.48,2.48,z],[-2.48,2.48,z]]
    ax.add_collection3d(Poly3DCollection([corners],facecolors='#eaf1f7',edgecolors='#99aabc',alpha=.08,linewidth=.7))
    ax.scatter(*a[k],s=45,color=blue,depthshade=False)
    ax.text(-2.9,-2.8,z,f'S{i}',color=ink,fontsize=11)
ax.text(2.4,-2.4,12,'橙色：前 4.132 mm',color=orange,fontsize=11)
ax.set_xlim(-2.8,2.8);ax.set_ylim(-2.8,2.8);ax.set_zlim(0,500)
ax.set_xticks([-2,0,2]);ax.set_yticks([-2,0,2]);ax.set_zticks([0,125,250,375,500])
ax.set_xlabel('x / mm',labelpad=2);ax.set_ylabel('y / mm',labelpad=2);ax.set_zlabel('z / mm',labelpad=11)
ax.set_box_aspect([1,1,1.8]);ax.view_init(elev=17,azim=-62)
ax.grid(False)
fig.text(.31,.095,'完整周期 500 mm；轴向显示比例已压缩，以便看清横向运动。',ha='center',fontsize=10.5,color='#536273')

ax2=fig.add_subplot(gs[1,3:],projection='3d')
ax2.set_title('③ 放大橙色段，加入铜线直径',fontsize=14,pad=16,color=ink)
def tube(loopid):
    z=np.linspace(0,4*q,180);u=(loopid+z/q)%484;k=np.floor(u).astype(int);t=(u-k)[:,None];cc=c[k]
    xy=((cc[:,3]*t+cc[:,2])*t+cc[:,1])*t+cc[:,0]
    v=((3*cc[:,3]*t+2*cc[:,2])*t+cc[:,1])/q
    T=np.column_stack([v,np.ones(len(v))]);T/=np.linalg.norm(T,axis=1)[:,None]
    n1=np.column_stack([np.ones(len(v)),np.zeros(len(v)),-v[:,0]]);n1/=np.linalg.norm(n1,axis=1)[:,None];n2=np.cross(T,n1)
    ang=np.linspace(0,2*np.pi,25)
    vec=n1[:,None,:]*np.cos(ang)[None,:,None]+n2[:,None,:]*np.sin(ang)[None,:,None]
    ctr=np.column_stack([xy,z]);return ctr[:,None,:]+d/2*vec,ctr
allpts=[]
for loopid,col in [(463,'#bcc5ce'),(462,blue)]:
    surf,ctr=tube(loopid);allpts.append(surf.reshape(-1,3))
    ax2.plot_surface(surf[:,:,0],surf[:,:,1],surf[:,:,2],rstride=2,cstride=1,color=col,linewidth=0,antialiased=True,shade=True,alpha=.95)
    if loopid==462:
        ax2.plot(ctr[:,0],ctr[:,1],ctr[:,2],lw=1,color='#0b3555',alpha=.8)
        for j in [0,2,4,6,8]:ax2.scatter(*a[j],s=22,color=orange,depthshade=False)
        ax2.text(*ctr[-1], '  A / Strand 44',color=ink,fontsize=11)
xyz=np.vstack(allpts);lo=xyz.min(axis=0);hi=xyz.max(axis=0)
ax2.set_xlim(lo[0]-.05,hi[0]+.05);ax2.set_ylim(lo[1]-.05,hi[1]+.05);ax2.set_zlim(lo[2],hi[2])
ax2.set_box_aspect(hi-lo+np.array([.1,.1,0]));ax2.view_init(elev=12,azim=-57);ax2.grid(False)
ax2.set_xlabel('');ax2.set_ylabel('');ax2.set_zlabel('z / mm',labelpad=8)
ax2.set_xticks([]);ax2.set_yticks([]);ax2.set_zticks([0,1,2,3,4]);ax2.tick_params(labelsize=10)
fig.text(.785,.095,'各轴等比例；铜径 125.63 μm\n灰线为相邻 Strand 22；管状外形是可视化示意。',ha='center',fontsize=10.5,color='#536273',linespacing=1.5)
fig.text(.5,.035,'同一截面处的蓝点 = 蓝色空间曲线与该截面的交点。前 4.132 mm 只占完整周期的 0.83%。',ha='center',fontsize=13,color=ink)
fig.savefig(out/'01-section-to-comsol-path.png',dpi=180)
fig.savefig(out/'01-section-to-comsol-path.pdf')
plt.close(fig)

# A clear plan projection of the same MPH-defined curve, with actual slice points.
fig,ax=plt.subplots(figsize=(9,8),dpi=180)
ax.plot(a[:,0],a[:,1],color=blue,lw=1.15,alpha=.8)
ax.plot(a[:9,0],a[:9,1],color=orange,lw=3,zorder=4)
ax.scatter(grid[:,0],grid[:,1],s=8,color='#bfc7cf',zorder=1)
for k in [28,100,164,244,390,620,890]:
    if np.min(np.ptp(a[k:k+8,:2],axis=0)) < .01:
        ax.annotate('',xy=a[k+7,:2],xytext=a[k,:2],arrowprops={'arrowstyle':'->','color':blue,'lw':1.2},zorder=3)
for i,k in enumerate(indices[:-1]):
    pt=a[k];ax.scatter(pt[0],pt[1],s=125,color=blue,edgecolors='white',zorder=5)
    dx=-.18 if pt[0]>1 else .14;dy=.18
    ax.annotate(f'S{i}: {pt[2]:.0f} mm',xy=pt[:2],xytext=(pt[0]+dx,pt[1]+dy),fontsize=13,ha='right' if dx<0 else 'left',color=ink)
ax.set(xlim=(-2.62,2.62),ylim=(-2.62,2.62),xlabel='x / mm',ylabel='y / mm');ax.set_aspect('equal')
ax.set_title('同一根 Strand 44 的全周期 x–y 投影',fontsize=18,pad=30)
fig.text(.5,.908,'原始 COMSOL 建模点重绘；箭头表示 z 增大的方向；橙色对应前 4.132 mm',ha='center',fontsize=10.5,color='#536273')
fig.text(.5,.023,'投影中各段位于不同 z；在任一横截面上，这根单丝只有一个中心点。S4 回到 S0 的位置。',ha='center',fontsize=10.5,color='#536273')
fig.subplots_adjust(top=.89,bottom=.10,left=.10,right=.96)
fig.savefig(out/'02-original-mph-xy-projection.png',dpi=180)
plt.close(fig)
print('Rendered matched sections, full path, finite-radius local model, and full XY projection.')
