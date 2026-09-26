from pathlib import Path
import sys, json, re
import numpy as np
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection

O=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent.parent
plt.rcParams.update({'font.family':'Microsoft YaHei','font.size':11,'axes.unicode_minus':False,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
def data(name):return np.genfromtxt(O/name,delimiter=',',names=True)
b=data('bundle_200kHz_results.csv'); ref=data('equal_area_solid_results.csv')
st=data('strand_currents_200kHz.csv'); low=data('low_frequency_check.csv'); lowst=data('strand_currents_10Hz.csv')
mesh=data('mesh_convergence.csv'); edge=data('outer_boundary_check.csv'); scan=data('frequency_sweep.csv')
old=json.loads((O/'reference_hexagonal'/'results_summary.json').read_text(encoding='utf-8'))
oldst=data('reference_hexagonal/strand_currents_200kHz.csv')
oldscan=data('reference_hexagonal/frequency_sweep.csv')
N=len(st); d=.14; sigma=5.8e7; mu=4*np.pi*1e-7; Irms=20.
I=st['I_real_rms_A']+1j*st['I_imag_rms_A']; Imag=abs(I); avg=Irms/N
Il=lowst['I_real_rms_A']+1j*lowst['I_imag_rms_A']
P=st['loss_W_per_m']; rings=st['radial_ring'].astype(int); outer=rings==10
A=N*np.pi*d*d/4; spacing=pdist(np.column_stack([st['x_mm'],st['y_mm']])).min()
Rb=1/(sigma*st['area_m2']); Ptransport=float(sum(Rb*Imag**2)); Pdc=float(400*b['Rdc_ohm'])
s={
 'N':N,'strand_diameter_mm':d,'copper_area_mm2':A,'Jdc_A_mm2':Irms/A,
 'delta_um':float(np.sqrt(1/(np.pi*2e5*mu*sigma))*1e6),
 'minimum_center_spacing_mm':float(spacing),'minimum_copper_gap_um':float((spacing-d)*1000),
 'bare_copper_envelope_diameter_mm':float(2*max(st['radius_mm'])+d),
 'Rdc_mohm':float(b['Rdc_ohm']*1e3),'Rac_mohm':float(b['Rac_ohm']*1e3),
 'Rac_over_Rdc':float(b['Rac_over_Rdc']),'loss_W_per_m':float(b['loss_W']),
 'solid_Rac_mohm':float(ref['Rac_ohm']*1e3),'solid_Rac_over_Rdc':float(ref['Rac_over_Rdc']),'solid_loss_W_per_m':float(ref['loss_W']),
 'hex_Rac_mohm':old['bundle_Rac_mohm'],'hex_Rac_over_Rdc':old['bundle_ratio'],'hex_loss_W_per_m':old['bundle_loss_W_m'],
 'Rac_change_vs_hex_percent':float(100*(b['Rac_ohm']*1e3/old['bundle_Rac_mohm']-1)),
 'Rac_change_vs_solid_percent':float(100*(b['Rac_ohm']/ref['Rac_ohm']-1)),
 'current_sum_real_A':float(sum(I).real),'current_sum_imag_A':float(sum(I).imag),
 'ideal_strand_current_mA':avg*1000,'center_current_mA':float(Imag[0]*1000),
 'minimum_current_mA':float(min(Imag)*1000),'maximum_current_mA':float(max(Imag)*1000),
 'outer_mean_current_mA':float(np.mean(Imag[outer])*1000),
 'outer_min_current_mA':float(np.min(Imag[outer])*1000),'outer_max_current_mA':float(np.max(Imag[outer])*1000),
 'outer_loss_percent':float(100*sum(P[outer])/sum(P)),
 'eta_complex_percent':float(100*max(abs(I-avg))/avg),'eta_magnitude_percent':float(100*max(abs(Imag-avg))/avg),
 'ideal_dc_loss_W_per_m':Pdc,'strand_imbalance_excess_loss_W_per_m':Ptransport-Pdc,'within_strand_excess_loss_W_per_m':float(sum(P)-Ptransport),
 'last_mesh_Rac_change_percent':float(100*abs(mesh['Rac_ohm'][-1]/mesh['Rac_ohm'][-2]-1)),
 'outer_boundary_Rac_change_percent':float(100*abs(edge['Rac_ohm']/b['Rac_ohm']-1)),
 'power_balance_relative_error':float(abs(b['Rac_voltage_ohm']/b['Rac_ohm']-1)),
 'low_frequency_Rac_Rdc':float(low['Rac_over_Rdc']),
 'low_frequency_current_magnitude_deviation_percent':float(100*max(abs(abs(Il)/avg-1))),
 'low_frequency_complex_current_deviation_percent':float(100*max(abs(Il/avg-1))),
}
assert N==331 and A>=5 and spacing>.15-1e-9
assert abs(sum(I)-20)<1e-7 and abs(sum(P)-b['loss_W'])<1e-8
assert s['power_balance_relative_error']<1e-7
assert s['last_mesh_Rac_change_percent']<.01 and s['outer_boundary_Rac_change_percent']<.01
assert abs(low['Rac_over_Rdc']-1)<1e-4 and s['low_frequency_current_magnitude_deviation_percent']<.001
assert max(abs(scan['I_real_rms_A']-20))<1e-7
mesh_counts=[]
log=(O/'verification'/'logs'/'Q2CircularBundle.log').read_text(encoding='utf-8',errors='replace')
for line in log.splitlines():
    if '自由度' in line:
        numbers=[int(x) for x in re.findall(r'\d+',line) if int(x)>10000]
        for count in numbers:
            if count not in mesh_counts:mesh_counts.append(count)
s['mesh_dof_counts']=mesh_counts
(O/'results_summary.json').write_text(json.dumps(s,ensure_ascii=False,indent=2),encoding='utf-8')
rows=[]
for k in range(11):
    mask=rings==k
    rows.append([k,sum(mask),np.mean(Imag[mask]),np.min(Imag[mask]),np.max(Imag[mask]),np.mean(Imag[mask]/(st['area_m2'][mask]*1e6)),sum(P[mask]),100*sum(P[mask])/sum(P)])
np.savetxt(O/'ring_statistics.csv',rows,delimiter=',',header='radial_ring,count,mean_Irms_A,min_Irms_A,max_Irms_A,mean_net_J_magnitude_A_mm2,loss_W_per_m,loss_fraction_percent',comments='',fmt='%.15g')
np.savetxt(O/'strand_analysis.csv',np.column_stack([st['strand_id'],rings,st['radius_mm'],Imag,np.angle(I,deg=True),Imag/(st['area_m2']*1e6),P,Rb*Imag**2,P-Rb*Imag**2]),delimiter=',',header='strand_id,radial_ring,radius_mm,I_rms_magnitude_A,I_phase_deg,net_J_magnitude_A_mm2,loss_W_per_m,uniform_same_current_loss_W_per_m,within_strand_excess_loss_W_per_m',comments='',fmt='%.15g')
np.savetxt(O/'strand_geometry.csv',np.column_stack([st['strand_id'],rings,st['ring_index'],st['x_mm'],st['y_mm'],np.full(N,.07),np.full(N,1000)]),delimiter=',',header='strand_id,radial_ring,ring_index,x_mm,y_mm,copper_radius_mm,length_mm',comments='',fmt='%.15g')
comparison=np.array([[s['solid_Rac_mohm'],s['solid_Rac_over_Rdc'],s['solid_loss_W_per_m']],[s['hex_Rac_mohm'],s['hex_Rac_over_Rdc'],s['hex_loss_W_per_m']],[s['Rac_mohm'],s['Rac_over_Rdc'],s['loss_W_per_m']]])
with (O/'shape_comparison.csv').open('w',encoding='utf-8') as f:
    f.write('configuration,copper_area_mm2,Rdc_mohm_per_m,Rac_mohm_per_m,Rac_over_Rdc,loss_W_per_m\n')
    for label,row in zip(['equal_area_solid','hexagonal_parallel_bundle','circular_parallel_bundle'],comparison):
        f.write(','.join([label]+[f'{v:.15g}' for v in [A,s['Rdc_mohm'],*row]])+'\n')

fig,axs=plt.subplots(1,2,figsize=(11.5,5.6),layout='constrained')
oldI=abs(oldst['I_real_rms_A']+1j*oldst['I_imag_rms_A'])*1000
vmax=max(max(oldI),max(Imag*1000))
for ax,ds,vals,title in [(axs[0],oldst,oldI,'原六边形排布'),(axs[1],st,Imag*1000,'圆形同心环排布')]:
    pc=PatchCollection([Circle((x,y),.07) for x,y in zip(ds['x_mm'],ds['y_mm'])],cmap='inferno',edgecolor='#46515b',linewidth=.3)
    pc.set_array(vals);pc.set_clim(0,vmax);ax.add_collection(pc)
    ax.set(aspect='equal',xlim=(-1.7,1.7),ylim=(-1.7,1.7),xlabel='x / mm',ylabel='y / mm',title=title)
fig.colorbar(pc,ax=axs,label='逐股净电流有效值 / mA',shrink=.85)
fig.suptitle('相同 331 股、相同丝径 · 200 kHz · 20 A RMS',fontsize=15)
fig.savefig(O/'current_comparison.png',dpi=210);plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
labels=['等铜截面\n实心线','六边形\n简单成束','圆形\n简单成束']
colors=['#758492','#d18541','#247e95']
for ax,index,title,ylabel in [(axs[0],0,'交流电阻比较','交流电阻 / (mΩ/m)'),(axs[1],2,'每米铜损比较','铜损 / (W/m)')]:
    vals=comparison[:,index];ax.bar(labels,vals,color=colors,width=.62)
    ax.set(title=title,ylabel=ylabel,ylim=(0,max(vals)*1.18));ax.grid(axis='y',alpha=.15)
    for j,v in enumerate(vals):ax.text(j,v+max(vals)*.025,f'{v:.4f}',ha='center')
fig.suptitle('均为 5.095349 mm² 铜截面、20 A 有效值、200 kHz',fontsize=14)
fig.savefig(O/'resistance_loss_comparison.png',dpi=210);plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(11,4.8),layout='constrained')
axs[0].scatter(oldst['radius_mm'],oldI,s=17,color='#d18541',label='六边形',alpha=.5)
axs[0].scatter(st['radius_mm'],Imag*1000,s=25,color='#247e95',label='圆形',alpha=.65)
axs[0].axhline(avg*1000,color='#58656c',ls='--',label=f'理想均分 {avg*1000:.3f} mA')
axs[0].set(xlabel='股线中心距束心 / mm',ylabel='股电流有效值 / mA',title='内外层电流仍有差别');axs[0].legend(fontsize=9);axs[0].grid(alpha=.2)
rr=np.array(rows);axs[1].bar(rr[:,0],rr[:,7],color='#247e95')
axs[1].set(xlabel='同心圆环编号（0 为中心）',ylabel='总铜损占比 / %',title=f'圆形外层 60 股承担 {s["outer_loss_percent"]:.2f}% 铜损')
axs[1].set_xticks(range(11));axs[1].grid(axis='y',alpha=.2)
fig.savefig(O/'radial_current_and_loss.png',dpi=210);plt.close(fig)

fig,ax=plt.subplots(figsize=(7.6,4.7),layout='constrained')
ax.plot(oldscan['frequency_Hz']/1000,oldscan['Rac_over_Rdc'],'s--',color='#d18541',label='六边形简单成束')
ax.plot(scan['frequency_Hz']/1000,scan['Rac_over_Rdc'],'o-',color='#247e95',label='圆形简单成束')
ax.set(xlabel='频率 / kHz',ylabel='Rac/Rdc',title='两种排布的 COMSOL 频率扫描')
ax.legend();ax.grid(alpha=.2)
fig.savefig(O/'frequency_comparison.png',dpi=210);plt.close(fig)
print(json.dumps(s,ensure_ascii=False,indent=2))
