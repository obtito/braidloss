"""Experimental round-wire impedance ROM; NOT accepted until 3D validation.

Exact cylindrical skin impedance and a truncated 2D harmonic scattering solve
are recomputed at the actual recursive geometry's cross sections. Currents are
unknown coupled circuit quantities. There is no square-grid impedance cache,
equal-current prescription, or empirically fitted pitch penalty.

Version 2 includes the exact cylindrical response to a magnetic field parallel
to the wire axis. Version 1 omitted this loss mechanism in its local curved-wire
extension. It is added with its analytic coefficient, without empirical fitting.

The extension to curved wires is an approximation. In particular the full
nonlocal 3D inductance is not represented by a 2D Green function. Its usefulness
must be decided by independent 3D A-phi comparisons, including rank agreement.
"""
import numpy as np
from math import comb
from scipy.special import jv,jvp
from scipy.linalg import solve
from threadpoolctl import threadpool_limits
from build_recursive_periodic import cycles

MU=4*np.pi*1e-7
SIGMA=5.8e7


def wire_response(a,frequency,modes=3):
    w=2*np.pi*frequency;k=(1-1j)*np.sqrt(w*MU*SIGMA/2)
    x=k*a
    zi=k*jv(0,x)/(2*np.pi*a*SIGMA*jv(1,x))
    n=np.arange(1,modes+1)
    q=x*jvp(n,x)/jv(n,x)
    scatter=(n-q)/(n+q)
    return zi,scatter


@threadpool_limits.wrap(limits=1)
def round_impedance(xy,a,frequency=200000.,outer_radius=.006,modes=3):
    """Parallel round conductors; Dirichlet outer circular magnetic boundary."""
    xy=np.asarray(xy);n=len(xy);q=n*modes
    z=xy[:,0]+1j*xy[:,1];d=z[:,None]-z[None,:]
    if n>1 and np.min(abs(d+np.eye(n)*1e30))<=2*a:
        raise ValueError('Overlapping circles in ROM slice')
    if max(abs(z))+a>=outer_radius:raise ValueError('Outside return boundary')
    d=d.copy();np.fill_diagonal(d,np.inf)
    ratio=a/d
    zi,tn=wire_response(a,frequency,modes)
    order=np.arange(1,modes+1)
    T=np.tile(tn,n)
    p=np.empty((n,modes,n),complex)
    u=np.empty((n,modes,n,modes),complex)
    image=np.zeros_like(u)
    R2=outer_radius**2
    den=R2-z[:,None]*z.conj()[None,:]
    # a/(z_j - R^2/conj(z_l)), including an image at infinity for z_l=0.
    ri=-a*z.conj()[None,:]/den
    for i,nn in enumerate(order):
        coef=(-1)**(nn+1)/nn
        p[:,i,:]=-MU/(4*np.pi)*coef*(ratio**nn-ri**nn)
        for h,mm in enumerate(order):
            u[:,i,:,h]=(-1)**nn*comb(mm+nn-1,nn)*ratio**(mm+nn)
            # Taylor coefficients of the reflected outgoing multipoles.
            for h0 in range(min(mm,nn)+1):
                image[:,i,:,h]-=(a**(mm+nn)*comb(mm,h0)*z[:,None]**(mm-h0)
                   /den**mm*comb(mm+nn-h0-1,nn-h0)*(z.conj()[None,:]/den)**(nn-h0))
    p=p.reshape(q,n);u=u.reshape(q,q);image=image.reshape(q,q)
    eye=np.eye(q,dtype=complex)
    matrix=np.block([[eye-T[:,None]*image.conj(),-T[:,None]*u.conj()],
                     [-T[:,None]*u,eye-T[:,None]*image]])
    rhs=np.vstack([T[:,None]*p.conj(),T[:,None]*p])
    s=solve(matrix,rhs,assume_a='gen',check_finite=False)
    kp=np.empty((n,n,modes),complex)
    for h,mm in enumerate(order):
        kp[:,:,h]=ratio**mm-(a*z.conj()[:,None]/(R2-z[None,:]*z.conj()[:,None]))**mm
    kp=kp.reshape(n,q)
    delta=kp@s[:q]+kp.conj()@s[q:]
    distance=abs(d);np.fill_diagonal(distance,a)
    L=MU/(2*np.pi)*np.log(abs(den)/(outer_radius*distance))
    Z=np.eye(n)*zi+1j*2*np.pi*frequency*(L+delta)
    reciprocity=np.max(abs(Z-Z.T))/max(abs(Z).max(),1e-30)
    if reciprocity>1e-7:raise ValueError('Multipole reciprocity failure: '+str(reciprocity))
    return .5*(Z+Z.T)


@threadpool_limits.wrap(limits=1)
def local_dipole_impedance(xyz,tangent,growth,a,frequency=200000.):
    """Reciprocal local straight-line dipole approximation in 3D.

    Length weighting comes from the local loss integral per axial length. The
    reciprocal average approximates mutual dipole interactions of skew axes;
    this is an explicit ROM approximation, not an exact 3D Green function.
    """
    n=len(xyz);t=np.asarray(tangent);g=np.asarray(growth)
    P=np.eye(3)[None]-t[:,:,None]*t[:,None,:]
    r=xyz[:,None]-xyz[None,:]
    rp=r-np.einsum('ijk,jk->ij',r,t)[:,:,None]*t[None]
    rho2=np.einsum('ijk,ijk->ij',rp,rp);np.fill_diagonal(rho2,np.inf)
    H=np.cross(t[None],rp)/(2*np.pi*rho2[:,:,None])
    H=np.einsum('iab,ijb->ija',P,H)
    H=H.transpose(0,2,1).reshape(3*n,n)*np.repeat(np.sqrt(g),3)[:,None]
    D=(2*rp[:,:,:,None]*rp[:,:,None,:]/rho2[:,:,None,None]-P[None])/(2*np.pi*rho2[:,:,None,None])
    D=np.einsum('iab,ijbc,jcd->ijad',P,D,P)
    D=D.transpose(0,2,1,3).reshape(3*n,3*n)
    weights=np.repeat(np.sqrt(g),3)
    D=D*weights[:,None]/weights[None,:]
    D=(D+D.T)/2
    alpha=2*np.pi*a*a*wire_response(a,frequency,1)[1][0]
    response=solve(np.eye(3*n)-alpha*D,alpha*H,check_finite=False)
    return 1j*2*np.pi*frequency*MU*(H.T@response)


def axial_polarizability(a,frequency=200000.):
    """Infinite round cylinder in uniform axial H; magnetic moment per length.

    H_inside(r)=H_applied*J0(k*r)/J0(k*a). Integrating the excess axial
    flux gives alpha_parallel = pi*a^2*(2*J1(x)/(x*J0(x))-1).
    With the exp(+i*omega*t) convention, -Im(alpha) is nonnegative.
    """
    x=(1-1j)*np.sqrt(np.pi*frequency*MU*SIGMA)*a
    return np.pi*a*a*(2*jv(1,x)/(x*jv(0,x))-1)


def local_axial_impedance(xyz,tangent,growth,a,frequency=200000.):
    """Leading local axial-field loss/reactance with exact round-wire response.

    The source H is the field of the locally straight filament tangents. Axial
    eddy-current interactions and the self-field of distant parts of the same
    curved strand are not resolved here; independent 3D validation is required.
    """
    t=np.asarray(tangent)
    r=xyz[:,None]-xyz[None,:]
    rp=r-np.einsum('ijk,jk->ij',r,t)[:,:,None]*t[None]
    rho2=np.einsum('ijk,ijk->ij',rp,rp);np.fill_diagonal(rho2,np.inf)
    H=np.cross(t[None],rp)/(2*np.pi*rho2[:,:,None])
    parallel=np.einsum('ijk,ik->ij',H,t)*np.sqrt(growth)[:,None]
    return 1j*2*np.pi*frequency*MU*axial_polarizability(a,frequency)*(parallel.T@parallel)


@threadpool_limits.wrap(limits=1)
def evaluate(c,cell_length,frequency=200000.,outer_radius=.006,slices=12,modes=3,
             orientation_correction=True,axial_correction=True):
    mapping=c.screw_mapping(cell_length,2*np.pi*cell_length/(c.pitches_mm[-1]*1e-3))
    if not mapping['geometry_mapping_pass']:raise ValueError('ROM periodic map failed')
    groups=cycles(mapping['strand_destination_to_source'])
    B=np.zeros((c.n,len(groups)))
    for j,group in enumerate(groups):B[group,j]=1
    count=B.sum(0)
    z=(np.arange(slices)+.5)/slices*cell_length
    xyz,t,_=c.at_physical_z(z)
    zi=wire_response(c.a,frequency,modes)[0]
    avg=np.zeros((c.n,c.n),complex)
    for k in range(slices):
        pos=xyz[:,k];tan=t[:,k];growth=1/tan[:,2]
        mat=round_impedance(pos[:,:2],c.a,frequency,outer_radius,modes)
        mat+=np.diag(zi*(growth-1))
        if orientation_correction:
            vertical=np.zeros_like(tan);vertical[:,2]=1
            mat+=local_dipole_impedance(pos,tan,growth,c.a,frequency)-local_dipole_impedance(pos,vertical,np.ones(c.n),c.a,frequency)
        if axial_correction:
            mat+=local_axial_impedance(pos,tan,growth,c.a,frequency)
        avg+=mat/slices
    reduced=B.T@avg@B
    currents=solve(reduced,count,assume_a='sym')
    total=count@currents
    strand=B@currents/total
    Z=1/total
    dc=c.diagnostics()['Rdc_length_parallel_ohm_per_m']
    return dict(Rac=float(Z.real),Rdc=float(dc),K=float(Z.real/dc),Z=Z,
       currents_fraction=strand,imbalance_complex=float(abs(strand*c.n-1).max()),
       passivity_min_eigenvalue=float(np.linalg.eigvalsh(avg.real).min()),
       status='experimental_ROM_requires_independent_3D_validation',slices=slices,modes=modes,
       version=2,axial_cylinder_response=axial_correction)


if __name__=='__main__':
    from recursive_geometry import RecursiveCable
    for f,p,area,ell,R in [((3,3),(4.,12.),6*9/343,4e-3/3,.001),((64,),(96.,),6.,.0005,.006)]:
        c=RecursiveCable(f,p,copper_area_mm2=area,gap_um=20. if len(f)==2 else 40.)
        ans=evaluate(c,ell,outer_radius=R,slices=12 if len(f)>1 else 1)
        print(f,{k:v for k,v in ans.items() if k!='currents_fraction'},flush=True)
