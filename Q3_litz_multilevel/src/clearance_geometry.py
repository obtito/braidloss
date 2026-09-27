"""Conservative inter-strand tube clearance for the periodic cubic-spline axes.

The certificate applies to the Python spline geometry, not to a subsequently
approximated CAD sweep. A CAD interpolation/meshing tolerance must be subtracted
separately. Swept circular tubes are contained in radius-a neighborhoods of axes.
"""
import numpy as np
from scipy.spatial import cKDTree


def segment_distances(p, u, q, v):
    """Exact distances between batches of closed straight line segments."""
    w=p-q
    aa=np.einsum('ij,ij->i',u,u);bb=np.einsum('ij,ij->i',u,v)
    cc=np.einsum('ij,ij->i',v,v)
    dd=np.einsum('ij,ij->i',u,w);ee=np.einsum('ij,ij->i',v,w)
    determinant=aa*cc-bb*bb
    valid=determinant>1e-12*aa*cc
    s=np.zeros(len(p));t=s.copy()
    s[valid]=(bb[valid]*ee[valid]-cc[valid]*dd[valid])/determinant[valid]
    t[valid]=(aa[valid]*ee[valid]-bb[valid]*dd[valid])/determinant[valid]
    interior=valid&(s>=0)&(s<=1)&(t>=0)&(t<=1)
    delta=w+s[:,None]*u-t[:,None]*v
    best=np.where(interior,np.einsum('ij,ij->i',delta,delta),np.inf)
    for point,start,direction,norm2 in [(p,q,v,cc),(p+u,q,v,cc),(q,p,u,aa),(q+v,p,u,aa)]:
        f=np.clip(np.einsum('ij,ij->i',point-start,direction)/norm2,0,1)
        delta=point-start-f[:,None]*direction
        best=np.minimum(best,np.einsum('ij,ij->i',delta,delta))
    return np.sqrt(np.maximum(best,0))


def periodic_curve_clearance(c, search_extra_gap_um=None):
    """Bound the minimum distance of different strand axes over a full repeat."""
    if c.n<2:return {'clearance_lower_bound_um':float('inf')}
    dt=c.t[1]-c.t[0]
    gap=(1.5*c.gap_um if search_extra_gap_um is None else search_extra_gap_um)*1e-6
    threshold=2*c.a+gap
    dzmin=c.evaluate(c.t,1)[...,2].min()
    if dzmin<=0:raise ValueError('Axial reversal')
    acceleration=np.linalg.norm(c.evaluate(c.t,2),axis=2).max()
    chord_error=acceleration*dt*dt/8
    # Different normal-plane offsets can shift two axes axially by different
    # amounts. Include that range when padding translated neighboring repeats.
    axial_offset_range=np.ptp(c.xyz[...,2]-c.t[None,:])+2*chord_error
    pad=int(np.ceil((threshold+axial_offset_range)/dt))+2
    times=np.arange(-pad,len(c.t)+pad)*dt
    pts=c.evaluate(times)
    start=pts[:,:-1].reshape(-1,3)
    vector=np.diff(pts,axis=1).reshape(-1,3)
    mids=start+vector/2
    strand=np.repeat(np.arange(c.n),len(times)-1)
    axial=np.tile(times[:-1],c.n)
    halfmax=np.linalg.norm(vector,axis=1).max()/2
    # For cubic splines C'' is linear on each interval. Its norm is convex;
    # the maximum at knot endpoints bounds the chord interpolation error.
    pairs=cKDTree(mids).query_pairs(threshold+2*halfmax,output_type='ndarray')
    keep=(strand[pairs[:,0]]!=strand[pairs[:,1]])&(
       ((axial[pairs[:,0]]>=0)&(axial[pairs[:,0]]<c.period))|
       ((axial[pairs[:,1]]>=0)&(axial[pairs[:,1]]<c.period)))
    pairs=pairs[keep]
    minimum=np.inf;best=None
    for k in range(0,len(pairs),100000):
        ij=pairs[k:k+100000];i,j=ij.T
        ds=segment_distances(start[i],vector[i],start[j],vector[j])
        b=ds.argmin()
        if ds[b]<minimum:minimum=float(ds[b]);best=ij[b]
    lower=min(minimum,threshold)-2*chord_error-2*c.a
    return dict(clearance_lower_bound_um=float(lower*1e6),
       closest_polygon_clearance_um=float((minimum-2*c.a)*1e6),
       spline_to_chord_bound_um=float(chord_error*1e6),
       closest_strands=None if best is None else strand[best].tolist(),
       candidate_segment_pairs=len(pairs),
       scope='All different spline axes, including adjacent periodic copies; excludes CAD approximation error and same-strand self-intersection')


if __name__=='__main__':
    from recursive_geometry import RecursiveCable
    import json
    for factors,pitches in [((7,7),(12.,36.)),((4,4,4),(16.,16.,96.)),((5,5,5),(20.,40.,120.))]:
        c=RecursiveCable(factors,pitches)
        print(json.dumps({'factors':factors,'pitches':pitches,**periodic_curve_clearance(c)}),flush=True)
