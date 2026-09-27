"""Read the existing braid; independently bound continuous finite-radius clearance.

Outputs are confined to this explanatory folder. No COMSOL model is modified.
All calculations below use millimetres. Bernstein convex-hull bounds cover
whole polynomial intervals (not just points); arithmetic is ordinary float64,
with explicit small margins, not a formal interval-arithmetic certificate.
"""
from pathlib import Path
import sys, json, math
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'litz_q3'
sys.path.insert(0, str(SOURCE / 'src'))
from geometry import Braid, hamiltonian_cycle

cfg = json.loads((SOURCE/'data/final_results.json').read_text(encoding='utf-8'))['recommended_within_0p2_percent_band']
b = Braid('full_exchange', 22, cfg['half_width_m'], cfg['period_m'])
c = b.splines[0].c[::-1].transpose(1, 0, 2) * 1000  # [484, ascending coefficient, xy]
n = b.n
q = cfg['period_m']*1000/n
d = cfg['diameter_um']/1000

def bernstein(a):
    deg = a.shape[-1]-1
    out = np.zeros_like(a)
    for k in range(deg+1):
        for i in range(k+1):
            out[...,k] += a[...,i]*math.comb(k,i)/math.comb(deg,i)
    return out

def split(be):
    deg=be.shape[-1]-1
    tmp=be.copy(); left=np.empty_like(be); right=np.empty_like(be)
    left[...,0]=tmp[...,0]; right[...,deg]=tmp[...,deg]
    for j in range(1,deg+1):
        tmp=(tmp[...,:-1]+tmp[...,1:])/2
        left[...,j]=tmp[...,0]; right[...,deg-j]=tmp[...,-1]
    return np.stack((left,right),axis=-2).reshape(-1,deg+1)

def square_norm(a):
    # last axes are ascending polynomial coefficients and x/y
    out=np.zeros(a.shape[:-2]+(2*a.shape[-2]-1,))
    for i in range(a.shape[-2]):
        for j in range(a.shape[-2]):
            out[...,i+j]+=np.sum(a[...,i,:]*a[...,j,:],axis=-1)
    return out

def global_bounds(be,depth):
    for _ in range(depth):
        be=split(be)
    return float(be.min()),float(be.max())

# Continuous speed and acceleration upper bounds from polynomial convex hulls.
v=c[:,1:,:]*np.arange(1,4)[None,:,None]/q
a=v[:,1:,:]*np.arange(1,3)[None,:,None]/q
_, speed2_hi=global_bounds(bernstein(square_norm(v)),7)
M=math.sqrt(speed2_hi)+1e-10
_, acc2_hi=global_bounds(bernstein(square_norm(a)),7)
acc_hi=math.sqrt(acc2_hi)+1e-10

# Eliminate distant pairs using displacement <= max speed * one axial step.
# This is valid throughout the complete phase interval, including endpoints.
ii,jj=np.triu_indices(n,1)
dist0=np.linalg.norm(c[ii,0,:]-c[jj,0,:],axis=1)
cut=0.18+2*M*q
near=dist0<=cut
ni,nj=ii[near],jj[near]
sq=square_norm(c[ni]-c[nj])
be=bernstein(sq)
# Remove any interval whose convex-hull lower bound exceeds a sampled upper
# bound on the global minimum; refine only intervals that can still matter.
sample_upper=min(float(np.min(np.polynomial.polynomial.polyval(t,sq.T))) for t in np.linspace(0,1,17))
active=be[be.min(axis=1)<=sample_upper+1e-12]
for _ in range(12):
    active=split(active)
    active=active[active.min(axis=1)<=sample_upper+1e-12]
delta_lo=math.sqrt(max(0,float(active.min())))-1e-8
far_lo=float(np.min(dist0[~near]-2*M*q)) if np.any(~near) else float('inf')
delta_lo=min(delta_lo,far_lo)

# Numerical stationary-point search identifies a useful close-up, independently
# of the above conservative all-pair continuous lower bound.
mid=np.polynomial.polynomial.polyval(.5,sq.T)
ids=np.argsort(mid)[:160]
best=(float('inf'),None)
for k in ids:
    der=sq[k,1:]*np.arange(1,7)
    roots=np.polynomial.polynomial.polyroots(der)
    ts=[0.,1.]+[float(z.real) for z in roots if abs(z.imag)<1e-8 and 0<z.real<1]
    for t in ts:
        val=float(np.polynomial.polynomial.polyval(t,sq[k]))
        if val<best[0]: best=(val,(int(ni[k]),int(nj[k]),t))
i,j,t=best[1]
spline=b.splines[0]
def closest_fun(x):
    tt,eta=x
    xy1=spline((i+tt)%n)*1000
    xy2=spline((j+tt+eta)%n)*1000
    return float(np.dot(xy1-xy2,xy1-xy2)+(q*eta)**2)
opt=minimize(closest_fun,[t,0],method='Nelder-Mead',options={'xatol':1e-13,'fatol':1e-16,'maxiter':2000})
center_bound=delta_lo/math.sqrt(1+M*M)

# Continuous box and radius bounds (slightly conservative Bernstein hulls).
maxabs=0
for axis in [0,1]:
    lo,hi=global_bounds(bernstein(c[:,:,axis]),8)
    maxabs=max(maxabs,abs(lo),abs(hi))
_, r2hi=global_bounds(bernstein(square_norm(c)),8)
boxside=2*maxabs+d
envelopeD=2*math.sqrt(r2hi)+d
length=cfg['mean_length_factor']
rows=[]
for enamel in [0,4,6,8,10,12]:
    rows.append({'enamel_radial_um':enamel,'outer_diameter_um':d*1000+2*enamel,
                 'clearance_bound_um':center_bound*1000-d*1000-2*enamel})
scales=[]
for enamel in [0,8]:
    for gap in [0,5]:
        need=d+(2*enamel+gap)/1000
        scale=need/math.sqrt(delta_lo**2-(need*M)**2)
        scales.append({'enamel_um':enamel,'gap_um':gap,'minimum_xy_scale_for_bound':scale})
inspect=b.inspect(samples=2048)
result={
    'scope':'Analytical tube-envelope screen using continuous polynomial interval bounds; no new COMSOL solid/EM solve',
    'arithmetic':'float64 Bernstein hulls with 0.01 nm distance margin, not formal interval arithmetic',
    'N':n,'period_mm':q*n,'station_advance_mm':q,'diameter_um':d*1000,
    'grid_spacing_um':(2*cfg['half_width_m']/21)*1e6,
    'all_strand_pairs':len(ii),'retained_candidate_pairs':int(near.sum()),
    'continuous_same_z_center_lower_mm':delta_lo,
    'continuous_max_slope_upper':M,'continuous_max_angle_upper_deg':math.degrees(math.atan(M)),
    'continuous_3D_center_distance_lower_mm':center_bound,
    'continuous_bare_surface_gap_lower_um':(center_bound-d)*1000,
    'continuous_local_bend_radius_lower_mm':1/acc_hi,
    'sampled_min_bend_radius_mm':inspect['min_bend_radius_m']*1000,
    'geometric_outer_copper_bend_strain_estimate':d/(2*inspect['min_bend_radius_m']*1000),
    'mean_length_factor':length,
    'box_side_mm':boxside,'round_envelope_D_mm':envelopeD,
    'nominal_copper_fill_square':6/boxside**2,
    'mean_copper_volume_fill_square':6*length/boxside**2,
    'nominal_copper_fill_round':6/(math.pi*(envelopeD/2)**2),
    'mean_copper_volume_fill_round':6*length/(math.pi*(envelopeD/2)**2),
    'insulation_cases':rows,'uniform_xy_scaling_sufficient_conditions':scales,
    'closeup':{'loop_indices':[i,j],'same_z_phase':t,'same_z_center_mm':math.sqrt(best[0]),
               'numerically_minimized_phase_eta':opt.x.tolist(),'numerically_minimized_center_mm':math.sqrt(opt.fun)},
    'not_verified':['specific enamel specification','manufacturing tolerance','solid CAD boolean validation',
                    'contact/friction/tension','manufacturing guide trajectories','three-dimensional volume FEM']
}
(HERE/'geometry_analysis.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
# Low-volume exact polynomial data for explanatory rendering.
small=Braid('full_exchange',4,.0015,16*q/1000)
payload={'parameters':result,'coeff484':np.round(c,12).tolist(),
         'coeff16':np.round(small.splines[0].c[::-1].transpose(1,0,2)*1000,12).tolist(),
         'loop_slot_ids':hamiltonian_cycle(22).tolist()}
(HERE/'render_data.json').write_text(json.dumps(payload,separators=(',',':')),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
