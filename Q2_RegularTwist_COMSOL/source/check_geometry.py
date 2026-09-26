from pathlib import Path
import numpy as np,json,sys
sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parent.parent
O=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else ROOT;(O/'source').mkdir(parents=True,exist_ok=True);(O/'verification/logs').mkdir(parents=True,exist_ok=True)
p=.15e-3; a=.07e-3; P=.04; alpha=2*np.pi/P;sigma=5.8e7
rows=[]
for k in range(11):
    n=1 if k==0 else 6*k
    for j in range(n):
        th=2*np.pi*j/n;r=k*p;s=np.sqrt(1+(alpha*r)**2)
        rows.append([len(rows)+1,k,j,r,th,r*np.cos(th),r*np.sin(th),s,np.arctan(alpha*r)*180/np.pi])
g=np.array(rows)
ii,jj=np.triu_indices(len(g),1)
r1=g[ii,3];r2=g[jj,3];theta=g[ii,4]-g[jj,4];dz=np.zeros(len(ii))
for _ in range(15):
    f=dz-alpha*r1*r2*np.sin(theta-alpha*dz)
    fp=1+alpha**2*r1*r2*np.cos(theta-alpha*dz)
    dz-=f/fp
assert alpha**2*max(g[:,3])**2 < 1, 'Squared distances must be strictly convex in axial offset.'
distance=np.sqrt(r1*r1+r2*r2-2*r1*r2*np.cos(theta-alpha*dz)+dz*dz)
gap=distance-2*a
assert min(gap)>0
A0=np.pi*a*a
summary={'N':len(g),'normal_wire_diameter_mm':2*a*1000,'twist_pitch_mm':P*1000,'twist_rate_per_m':alpha,'turns_per_axial_meter':1/P,'outer_helix_angle_deg':float(max(g[:,8])),'outer_wire_length_per_axial_meter_m':float(max(g[:,7])),'normal_copper_area_mm2':len(g)*A0*1e6,'total_horizontal_copper_area_mm2':float(sum(g[:,7])*A0*1e6),'copper_volume_per_axial_meter_cm3':float(sum(g[:,7])*A0*1e6),'copper_length_volume_increase_percent':float((np.mean(g[:,7])-1)*100),'Rdc_slender_wire_parallel_ohm_per_m':float(1/(sigma*A0*np.sum(1/g[:,7]))),'minimum_physical_copper_gap_um':float(min(gap)*1e6),'nearest_strand_pair':[int(ii[np.argmin(gap)]+1),int(jj[np.argmin(gap)]+1)],'nearest_pair_axial_offset_um':float(dz[np.argmin(gap)]*1e6)}
np.savetxt(O/'strand_geometry.csv',g,delimiter=',',header='strand_id,radial_ring,ring_index,helix_radius_m,phase_rad,x_at_z0_m,y_at_z0_m,strand_length_per_axial_meter_m,helix_angle_deg',comments='',fmt='%.15g')
(O/'geometry_checks.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
