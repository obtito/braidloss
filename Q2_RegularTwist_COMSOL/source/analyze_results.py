"""Reproduce numerical checks and figures from exported COMSOL data."""
from pathlib import Path
import csv, json, sys, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

sys.stdout.reconfigure(encoding='utf-8')
OUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                     'axes.unicode_minus':False,'font.size':11,'figure.dpi':160})
def row(path):
    with open(path,encoding='utf-8-sig',newline='') as f:
        return {k:float(v) for k,v in next(csv.DictReader(f)).items()}
def rows(path):
    with open(path,encoding='utf-8-sig',newline='') as f:
        return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def write_csv(name,items):
    with open(OUT/name,'w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(items[0]));w.writeheader();w.writerows(items)

main=row(OUT/'regular_twist_200kHz_results.csv')
low=row(OUT/'verification/low_frequency_10Hz.csv')
zero=row(OUT/'verification/zero_twist_limit.csv')
boundary=row(OUT/'verification/outer_boundary_40mm.csv')
mesh=[row(OUT/f'verification/mesh_{k}.csv') for k in (3,5,8)]
for r in [main,low,zero,boundary,*mesh]:
    assert r['air_regularization_S_per_m']==0, 'Final validation must use exact insulation.'
    assert abs(r['I_real_rms_A']-20)<2e-4
    assert abs(r['I_imag_rms_A'])<2e-4
    assert abs(r['power_balance_relative_residual'])<1e-6
strand=rows(OUT/'strand_currents_200kHz.csv')
assert len(strand)==331
geom=json.loads((OUT/'geometry_checks.json').read_text(encoding='utf-8'))
ref=json.loads((OUT/'references/circular_results_summary.json').read_text(encoding='utf-8'))
current=np.array([r['I_real_rms_A']+1j*r['I_imag_rms_A'] for r in strand])
loss=np.array([r['loss_W_per_m'] for r in strand])
ring=np.array([int(r['radial_ring']) for r in strand])
assert abs(current.sum()-20)<2e-4
assert abs(loss.sum()-main['loss_W'])<1e-7
amp=np.abs(current)
angle=np.angle(current,deg=True)
ideal=20/331
ringdata=[]
for k in range(11):
    ix=ring==k
    ringdata.append(dict(ring=k,radius_mm=.15*k,count=int(ix.sum()),
        current_mean_rms_mA=float(amp[ix].mean()*1000),
        current_min_rms_mA=float(amp[ix].min()*1000),
        current_max_rms_mA=float(amp[ix].max()*1000),
        mean_real_current_mA=float(current[ix].real.mean()*1000),
        mean_imag_current_mA=float(current[ix].imag.mean()*1000),
        copper_loss_W=float(loss[ix].sum()),
        copper_loss_share_percent=float(loss[ix].sum()/loss.sum()*100),
        length_factor=math.sqrt(1+(2*math.pi/.04*k*.00015)**2)))
write_csv('ring_statistics.csv',ringdata)
analysis=[]
for r,a,p in zip(strand,amp,angle):
    analysis.append({**r,'I_magnitude_rms_mA':float(a*1000),'phase_deg':float(p),
                     'length_factor':math.sqrt(1+(2*math.pi/.04*r['radius_mm']*.001)**2)})
write_csv('strand_analysis.csv',analysis)
meshdata=[]
for k,r in zip((3,5,8),mesh):
    meshdata.append({'mesh_div':k,'strand_max_element_um':140/k,**r})
write_csv('mesh_convergence.csv',meshdata)
comparisons=[
 dict(model='等名义铜截面实心圆线',Rdc_mohm_per_m=ref['Rdc_mohm'],Rac_mohm_per_m=ref['solid_Rac_mohm'],
      Rac_over_Rdc=ref['solid_Rac_over_Rdc'],copper_loss_W_per_m=ref['solid_loss_W_per_m']),
 dict(model='圆形简单成束',Rdc_mohm_per_m=ref['Rdc_mohm'],Rac_mohm_per_m=ref['Rac_mohm'],
      Rac_over_Rdc=ref['Rac_over_Rdc'],copper_loss_W_per_m=ref['loss_W_per_m']),
 dict(model='规则绞合，绞距40 mm',Rdc_mohm_per_m=main['Rdc_ohm']*1000,
      Rac_mohm_per_m=main['Rac_ohm']*1000,Rac_over_Rdc=main['Rac_over_Rdc'],copper_loss_W_per_m=main['loss_W'])]
write_csv('comparison.csv',comparisons)
summary={**geom,
    'frequency_Hz':200000,'current_rms_A':20,
    'Rdc_mohm_per_axial_m':main['Rdc_ohm']*1000,
    'Rac_mohm_per_axial_m':main['Rac_ohm']*1000,
    'Rac_over_Rdc':main['Rac_over_Rdc'],
    'copper_loss_W_per_axial_m':main['loss_W'],
    'Rac_change_vs_circular_percent':100*(main['Rac_ohm']*1000/ref['Rac_mohm']-1),
    'Rac_change_vs_solid_percent':100*(main['Rac_ohm']*1000/ref['solid_Rac_mohm']-1),
    'center_current_rms_mA':float(amp[0]*1000),
    'outer_mean_current_rms_mA':float(amp[ring==10].mean()*1000),
    'outer_min_current_rms_mA':float(amp[ring==10].min()*1000),
    'outer_max_current_rms_mA':float(amp[ring==10].max()*1000),
    'outer_loss_share_percent':float(loss[ring==10].sum()/loss.sum()*100),
    'current_sum_real_A':float(current.real.sum()),
    'current_sum_imag_A':float(current.imag.sum()),
    'max_phasor_deviation_from_equal_current_percent':float(np.max(np.abs(current-ideal))/ideal*100),
    'last_mesh_Rac_change_percent':abs(mesh[-1]['Rac_ohm']/mesh[-2]['Rac_ohm']-1)*100,
    'outer_boundary_Rac_change_percent':abs(boundary['Rac_ohm']/main['Rac_ohm']-1)*100,
    'zero_twist_reference_Rac_change_percent':abs(zero['Rac_ohm']/(ref['Rac_mohm']/1000)-1)*100,
    'low_frequency_Rac_over_Rdc':low['Rac_over_Rdc'],
    'main_power_balance_relative_residual':main['power_balance_relative_residual'],
    'normal_area_Jdc_nominal_A_per_mm2':20/geom['normal_copper_area_mm2'],
    'slender_DC_center_J_A_per_mm2':20*main['Rdc_ohm']*5.8e7/1e6,
    'final_insulation_conductivity_S_per_m':0,
    'method':'helicoidal reduction, 3-component magnetic vector potential, exact insulating gaps',
    'reference_note':'Same normal wire diameters and strand counts; twisted copper volume per axial meter is 1.507% greater.'
}
(OUT/'results_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')

fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
names=['等铜截面实心线','圆形简单成束','规则绞合 40 mm']
colors=['#8c99a5','#3688b4','#d77c28']
for ax,key,ylabel in zip(axs,['Rac_mohm_per_m','Rac_over_Rdc'],['交流电阻 / (mΩ/m)','交流 / 直流电阻比']):
    vals=[r[key] for r in comparisons]
    ax.bar(names,vals,color=colors,width=.55)
    for i,v in enumerate(vals):ax.text(i,v+.03,f'{v:.4f}',ha='center',va='bottom')
    ax.set_ylabel(ylabel);ax.set_ylim(0,max(vals)*1.15);ax.spines[['right','top']].set_visible(False)
fig.suptitle('200 kHz · 20 A RMS · 按线束轴向 1 m 评估')
fig.text(.5,-.025,'331 股 × 0.14 mm：规则绞合用铜量增加 1.507%，各股保持原径向层。',ha='center',fontsize=10)
fig.savefig(OUT/'resistance_comparison.png',bbox_inches='tight');plt.close(fig)

straight_strand=rows(OUT/'references/straight_strand_currents_200kHz.csv')
straight_i=np.array([abs(r['I_real_rms_A']+1j*r['I_imag_rms_A'])*1000 for r in straight_strand])
straight_ring=np.array([int(r['radial_ring']) for r in straight_strand])
fig,axs=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
x=np.arange(11)*.15
axs[0].plot(x,[straight_i[straight_ring==k].mean() for k in range(11)],'o-',color='#3688b4',label='简单成束')
axs[0].plot(x,[r['current_mean_rms_mA'] for r in ringdata],'s-',color='#d77c28',label='规则绞合 40 mm')
axs[0].axhline(ideal*1000,ls='--',color='#6d6d6d',label='理想等分电流')
axs[0].set(xlabel='股线螺旋半径 / mm',ylabel='每股电流有效值幅值 / mA')
axs[0].legend();axs[0].grid(alpha=.15)
axs[1].bar(x,[r['copper_loss_share_percent'] for r in ringdata],width=.095,color='#d77c28')
axs[1].set(xlabel='股线螺旋半径 / mm',ylabel='该层占总铜损 / %')
axs[1].text(.04,.92,f'最外层 60 股承担 {summary["outer_loss_share_percent"]:.2f}% 铜损',transform=axs[1].transAxes)
for ax in axs:ax.spines[['right','top']].set_visible(False)
fig.suptitle('规则绞合仍存在明显的内外层电流差异')
fig.savefig(OUT/'radial_current_loss_comparison.png',bbox_inches='tight');plt.close(fig)

# Scientific visualization of the actual helical centerline equations.
fig=plt.figure(figsize=(11,5.6),layout='constrained')
ax=fig.add_subplot(121,projection='3d')
z=np.linspace(0,10,100)
for r in strand:
    k=int(r['radial_ring']);n=1 if k==0 else 6*k
    theta=2*np.pi*r['ring_index']/n+2*np.pi*z/40
    rr=r['radius_mm']
    highlight=(k in (0,5,10) and int(r['ring_index'])==0)
    col={0:'#333333',5:'#3688b4',10:'#d77c28'}.get(k,'#999999') if highlight else '#a7b1bb'
    ax.plot(rr*np.cos(theta),rr*np.sin(theta),z,color=col,lw=2.2 if highlight else .38,alpha=1 if highlight else .35)
ax.set(xlabel='x / mm',ylabel='y / mm',zlabel='轴向 z / mm',title='真实螺旋中心线：前 10 mm 局部')
ax.set_box_aspect((1,1,2.2));ax.view_init(elev=22,azim=35)
ax.set_xticks([-1.5,0,1.5]);ax.set_yticks([-1.5,0,1.5])
ax2=fig.add_subplot(122)
zz=np.linspace(0,40,300)
for rr,col in [(0,'#333333'),(.75,'#3688b4'),(1.5,'#d77c28')]:
    ax2.plot(zz,np.full_like(zz,rr),color=col,lw=2,label=f'初始半径 {rr:g} mm')
ax2.set(xlabel='轴向 z / mm',ylabel='距线束中心的半径 / mm',title='一个完整绞距中，各股不换径向层',ylim=(-.1,1.8))
ax2.legend(loc='lower right',bbox_to_anchor=(.99,.08));ax2.spines[['right','top']].set_visible(False)
fig.suptitle('规则绞合：各股绕中心转动，保持原有半径')
fig.savefig(OUT/'helical_geometry_explanation.png',bbox_inches='tight');plt.close(fig)
print(json.dumps(summary,ensure_ascii=False,indent=2))
