from pathlib import Path
import json
import numpy as np
from scipy.special import jv
from scipy.optimize import brentq
from scipy.integrate import simpson
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT=Path(__file__).resolve().parent.parent
plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
sigma=5.8e7; mu=4*np.pi*1e-7; a=1e-3; I=20.; L=1.; f=2e5
delta=np.sqrt(1/(np.pi*f*mu*sigma)); rdc=L/(sigma*np.pi*a*a)
def theory(freq, r):
    de=np.sqrt(1/(np.pi*freq*mu*sigma)); k=(1-1j)/de
    J=I*k/(2*np.pi*a)*jv(0,k*r)/jv(1,k*a)
    Z=L*k/(2*np.pi*a*sigma)*jv(0,k*a)/jv(1,k*a)
    return J,Z
raw=np.genfromtxt(OUT/'radial_200kHz.csv',delimiter=',',names=True)
mesh=np.genfromtxt(OUT/'mesh_convergence.csv',delimiter=',',names=True)
sweep=np.genfromtxt(OUT/'frequency_sweep.csv',delimiter=',',names=True)
axial=np.genfromtxt(OUT/'axial_mesh_check.csv',delimiter=',',names=True)
r=raw['r_m']; j=raw['J_rms_A_m2']; jc=raw['Jz_real_rms_A_m2']+1j*raw['Jz_imag_rms_A_m2']
ja,z=theory(f,r); rac=float(mesh['Rac_ohm'][-1]); ratio=float(mesh['Rac_over_Rdc'][-1]); loss=float(mesh['loss_W'][-1])
dfem=a-np.interp(j[-1]/np.e,j,r)
dexact=a-brentq(lambda rr:abs(theory(f,rr)[0])/abs(ja[-1])-1/np.e,0,a,xtol=1e-15)
rowdiff=100*np.abs(np.diff(mesh['Rac_ohm'])/mesh['Rac_ohm'][1:])
mesherr=100*np.abs(mesh['Rac_ohm']/z.real-1)
freqtheory=np.array([theory(ff,a)[1].real for ff in sweep['frequency_Hz']])
current=complex(mesh['I_real_rms_A'][-1],mesh['I_imag_rms_A'][-1])
summary=dict(frequency_Hz=f,delta_um=delta*1e6,Rdc_mohm=rdc*1000,Rac_FEM_mohm=rac*1000,Rac_exact_mohm=z.real*1000,ratio_FEM=ratio,loss_W=loss,skin_layer_efold_FEM_um=dfem*1e6,skin_layer_efold_Bessel_um=dexact*1e6,Rac_error_percent=float(mesherr[-1]),last_mesh_change_percent=float(rowdiff[-1]),I_real_rms_A=current.real,I_imag_rms_A=current.imag,radial_integral_current_real_A=float(simpson(2*np.pi*r*jc,x=r).real),surface_Jrms_A_mm2=float(j[-1]*1e-6),center_Jrms_A_mm2=float(j[0]*1e-6),max_profile_abs_error_normalized_to_surface_percent=float(100*np.max(np.abs(jc-ja))/abs(ja[-1])),axial_refinement_change_percent=float(100*abs(axial['Rac_ohm']/rac-1)),max_frequency_sweep_Rac_error_percent=float(100*np.max(np.abs(sweep['Rac_ohm']/freqtheory-1))))
assert abs(current-20)<1e-5, current
assert mesherr[-1]<0.01
assert summary['max_profile_abs_error_normalized_to_surface_percent']<0.1
(OUT/'results_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
np.savetxt(OUT/'radial_theory_comparison.csv',np.column_stack([r,j,abs(ja),jc.real,jc.imag,ja.real,ja.imag]),delimiter=',',header='r_m,Jrms_FEM_A_m2,Jrms_Bessel_A_m2,Jz_real_FEM,Jz_imag_FEM,Jz_real_Bessel,Jz_imag_Bessel',comments='',fmt='%.15g')
np.savetxt(OUT/'frequency_theory_comparison.csv',np.column_stack([sweep['frequency_Hz'],sweep['Rac_ohm'],freqtheory,100*(sweep['Rac_ohm']/freqtheory-1)]),delimiter=',',header='frequency_Hz,Rac_FEM_ohm,Rac_Bessel_ohm,relative_error_percent',comments='',fmt='%.15g')

# Reconstruct the circular section by rotational symmetry, using FEM radial samples only.
xy=np.linspace(-a,a,801); X,Y=np.meshgrid(xy,xy); R=np.hypot(X,Y)
JJ=np.ma.array(np.interp(R,r,j)/1e6,mask=R>a)
fig,ax=plt.subplots(figsize=(6.6,5.6),layout='constrained')
im=ax.pcolormesh(X*1e3,Y*1e3,JJ,cmap='inferno',shading='auto',rasterized=True,vmin=0,vmax=j[-1]/1e6)
ax.set(aspect='equal',xlabel='x / mm',ylabel='y / mm',title='200 kHz 铜导体截面电流密度')
ax.text(.02,.02,'COMSOL 有限元数据按轴对称旋转重建\n20 A 有效值 · D = 2 mm',transform=ax.transAxes,color='black',fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
fig.colorbar(im,ax=ax,label=r'$|J|_{\rm rms}$ / (A/mm²)');fig.savefig(OUT/'cross_section_FEM.png',dpi=240);plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
axs[0].plot(r*1e3,abs(ja)/1e6,color='#185c9e',lw=2,label='圆柱导体 Bessel 解析解')
axs[0].plot(r[::40]*1e3,j[::40]/1e6,'o',color='#e46a35',ms=4,mfc='white',label='COMSOL 有限元')
axs[0].axhline(I/(np.pi*a*a)/1e6,color='#738393',ls='--',label='直流均匀电流密度')
axs[0].set(xlabel='半径 r / mm',ylabel='电流密度有效值 / (A/mm²)',title='径向分布与独立解析解');axs[0].legend(fontsize=9);axs[0].grid(alpha=.2)
depth=(a-r)*1e6
axs[1].plot(depth,j/j[-1],color='#e46a35',lw=2,label='COMSOL：圆柱导体')
axs[1].plot(depth,np.exp(-(a-r)/delta),'--',color='#185c9e',label='平面近似 exp(−深度/δ)')
axs[1].axhline(1/np.e,color='#778899',ls=':');axs[1].axvline(delta*1e6,color='#185c9e',ls=':',label=f'理论 δ = {delta*1e6:.2f} μm')
axs[1].axvline(dfem*1e6,color='#e46a35',ls=':',label=f'圆柱 1/e 层厚 = {dfem*1e6:.2f} μm')
axs[1].set(xlim=(0,600),ylim=(0,1.05),xlabel='距表面的深度 / μm',ylabel='相对于表面电流密度的幅值',title='趋肤深度与曲率修正');axs[1].legend(fontsize=9);axs[1].grid(alpha=.2)
fig.savefig(OUT/'radial_validation.png',dpi=220);plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(11,4.4),layout='constrained')
axs[0].plot(mesh['radial_elements'],mesh['Rac_ohm']*1e3,'o-',color='#185c9e')
axs[0].axhline(z.real*1e3,color='#e46a35',ls='--',label='Bessel 解析解')
axs[0].set(xscale='log',xlabel='径向单元数',ylabel='1 m 交流电阻 / mΩ',title='五档网格收敛');axs[0].set_xticks(mesh['radial_elements'],[str(int(n)) for n in mesh['radial_elements']]);axs[0].legend();axs[0].grid(alpha=.2)
axs[1].loglog(mesh['radial_elements'],mesherr,'o-',color='#e46a35');axs[1].set(xlabel='径向单元数',ylabel='相对解析解的电阻误差 / %',title='离散误差随网格加密下降');axs[1].grid(alpha=.2,which='both')
fig.savefig(OUT/'mesh_convergence.png',dpi=220);plt.close(fig)

fig,ax=plt.subplots(figsize=(7.3,4.5),layout='constrained')
ff=np.geomspace(5e4,1e6,200); rr=[theory(fi,a)[1].real/rdc for fi in ff]
ax.plot(ff/1e3,rr,color='#185c9e',label='Bessel 解析解');ax.plot(sweep['frequency_Hz']/1e3,sweep['Rac_over_Rdc'],'o',color='#e46a35',label='COMSOL 有限元')
ax.set(xlabel='频率 / kHz',ylabel=r'$R_{ac}/R_{dc}$',title='50 kHz–1 MHz 频率扫描');ax.legend();ax.grid(alpha=.2)
fig.savefig(OUT/'frequency_sweep.png',dpi=220);plt.close(fig)

print(json.dumps(summary,ensure_ascii=False,indent=2))
