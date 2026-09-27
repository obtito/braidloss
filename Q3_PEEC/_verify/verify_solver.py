"""
Independent numpy reimplementation of Question/Q3_PEEC/q3_solver.m (frozen v3).
Purpose: verify (do not modify core) the T2' calibration numbers from first
principles, reproducing the frozen regression anchors to machine precision.
Mirrors q3_solver.m line-for-line for scheme in {straight, twist, braid,
counter_braid}, with optional analytic intra-strand proximity (cfg.intra).
"""
import numpy as np
from cmath import sqrt as csqrt

def besselj0(x):
    x = complex(x); term = 1.0 + 0j; s = 1.0 + 0j
    x2 = -x*x/4.0; k = 1
    while abs(term) > 1e-18*max(1.0, abs(s)):
        term *= x2/(k*k); s += term; k += 1
        if k > 200: break
    return s

def besselj1(x):
    x = complex(x); half = x/2.0
    term = half; s = half; k = 1
    x2 = -half*half
    while abs(term) > 1e-18*max(1.0, abs(s)):
        term *= x2/(k*(k+1)); s += term; k += 1
        if k > 200: break
    return s

def q3_solver(cfg):
    sigma=cfg['sigma']; mu=4*np.pi*1e-7; f=cfg['f']; w=2*np.pi*f; Itot=cfg['I']
    a=0.5*cfg['d_mm']*1e-3; As=np.pi*a*a
    delta=np.sqrt(2/(w*mu*sigma)); kc=(1-1j)/delta
    Zint=kc/(2*np.pi*a*sigma)*besselj0(kc*a)/besselj1(kc*a)   # ohm/m complex
    Rdc_strand=1/(sigma*As)
    intra=cfg.get('intra', False)

    # ---- Q2 circular layout (1,6,...,60 ; r_k = p*k) ----
    K0=cfg['ringmax']; cap=6*np.arange(0,K0+1); cap[0]=1
    ring=[]; th0=[]
    for k in range(0,K0+1):
        c=int(cap[k]); j=np.arange(0,c)
        ring += [k]*c
        th0 += list(2*np.pi*j/c)
    ring=np.array(ring); th0=np.array(th0,dtype=float)
    N=ring.size; p=cfg['p_mm']*1e-3; rr=p*np.arange(0,K0+1)
    r0=p*ring; x=r0*np.cos(th0); y=r0*np.sin(th0)

    # ---- trajectory ----
    sch=cfg['scheme']
    if sch=='straight':
        K=1; X=x[:,None]; Y=y[:,None]; VX=np.zeros((N,1)); VY=np.zeros((N,1)); S=np.ones((N,1))
    elif sch=='twist':
        K=1; om=2*np.pi/(cfg['pitch_mm']*1e-3)
        X=x[:,None]; Y=y[:,None]; VX=(-om*y)[:,None]; VY=(om*x)[:,None]
        S=np.sqrt(1+(om*r0)**2)[:,None]
    elif sch in ('braid','counter_braid'):
        K=max(4,round(cfg['K'])); Lam=cfg['Lambda_mm']*1e-3
        tanA=np.tan(np.deg2rad(cfg['alpha_deg'])); beta=tanA/rr[-1]
        Fv=np.cumsum(cap)/N
        phi=((np.arange(1,N+1))-0.5)/N
        fam=(2*(ring%2)-1) if sch=='counter_braid' else np.ones(N)
        zm=((np.arange(0,K))+0.5)*(Lam/K)
        X=np.zeros((N,K)); Y=np.zeros((N,K)); VX=np.zeros((N,K)); VY=np.zeros((N,K)); S=np.zeros((N,K))
        drdzmag=(rr[-1]-rr[0])/(Lam/2)
        for m in range(K):
            u=np.mod(phi+zm[m]/Lam,1.0)
            tri=2*np.minimum(u,1-u)                 # rising then falling, one cycle/Lam
            lay=layerq(tri,Fv)
            cnt=np.bincount(lay, minlength=K0+1)     # occupancy conservation assert
            assert np.array_equal(cnt, cap), f"occupancy {cnt} != caps at section {m}"
            r=p*lay
            th=th0+fam*beta*zm[m]
            for k in range(0,K0+1):
                idx=np.where(lay==k)[0]
                o=np.argsort(th[idx], kind='stable')
                th[idx[o]]=2*np.pi*np.arange(idx.size)/cap[k]   # music-chairs slots
            drdz=drdzmag*(1-2*(u>=0.5).astype(float)); dthdz=fam*beta
            X[:,m]=r*np.cos(th); Y[:,m]=r*np.sin(th)
            VX[:,m]=drdz*np.cos(th)-r*np.sin(th)*dthdz
            VY[:,m]=drdz*np.sin(th)+r*np.cos(th)*dthdz
            S[:,m]=np.sqrt(1+drdz**2+(r*dthdz)**2)
    else:
        raise ValueError(f"unknown scheme {sch}")

    # ---- mutual matrix, per-section averaging; assert >= 2a (v3: no clamp) ----
    gap=2*a
    Macc=np.zeros((N,N)); sacc=np.zeros(N); mind2=np.inf
    for kk in range(K):
        xk=X[:,kk]; yk=Y[:,kk]; vx=VX[:,kk]; vy=VY[:,kk]; sk=S[:,kk]
        dx=xk[:,None]-xk[None,:]; dy=yk[:,None]-yk[None,:]; d2=dx*dx+dy*dy
        dtmp=d2.copy(); np.fill_diagonal(dtmp,np.inf)
        mind2=min(mind2, dtmp.min())
        if mind2 < gap*gap:
            raise RuntimeError(f"spacing {np.sqrt(mind2)*1e6:.1f}um < 2a={gap*1e6:.1f} at section {kk}")
        d2=np.maximum(d2, gap*gap); np.fill_diagonal(d2, a*a)
        Mseg=-0.5*np.log(d2)
        tdot=(1+vx[:,None]*vx[None,:]+vy[:,None]*vy[None,:])/(sk[:,None]*sk[None,:])
        Macc=Macc+tdot*Mseg; sacc=sacc+sk
    sbar=sacc/K; Mavg=Macc/K

    # ---- parallel common-voltage PEEC solve (Q2 convention) ----
    Zmat=np.diag(Zint*sbar)+1j*(w*mu/(2*np.pi))*Mavg
    wi=1.0/sbar
    zsol=np.linalg.solve(Zmat, wi)
    Vdrive=Itot/(wi@zsol)
    Iwire=Vdrive*zsol
    Iz=Iwire/sbar
    Ibar=Itot/N
    eta_c=np.max(np.abs(Iz-Ibar))/Ibar*100
    eta_m=np.max(np.abs(np.abs(Iz)-Ibar))/Ibar*100

    P_peec=np.real(Zint)*np.sum(sbar*np.abs(Iwire)**2)
    Rac_peec=P_peec/Itot**2
    Rdc_bundle=Rdc_strand/np.sum(1.0/sbar)
    J_peec=Rac_peec/Rdc_bundle

    # ---- optional analytic intra-strand proximity term (uniform-field, /4, RMS) ----
    P_intra=0.0
    if intra:
        Cintra=sigma*np.pi*a**4*w**2/4
        B2=np.zeros(N)
        for kk in range(K):
            xk=X[:,kk]; yk=Y[:,kk]
            dx=xk[:,None]-xk[None,:]; dy=yk[:,None]-yk[None,:]
            d2=np.maximum(dx*dx+dy*dy, gap*gap); np.fill_diagonal(d2, np.inf)
            B=(mu/(2*np.pi))*(((-dy+1j*dx)/d2)@Iwire)
            B2=B2+np.abs(B)**2
        P_intra=Cintra*np.sum(B2/K)
    Rac_full=(P_peec+P_intra)/Itot**2
    J_full=Rac_full/Rdc_bundle

    isOuter=(ring==K0)
    P_ring=np.real(Zint)*sbar*np.abs(Iwire)**2
    ic=np.where(ring==0)[0]
    ic=ic[0] if ic.size>0 else 0

    return dict(scheme=sch, J_peec=J_peec, J_full=J_full,
        Rac_peec_mohm=Rac_peec*1e3, Rac_full_mohm=Rac_full*1e3, Rdc_mohm=Rdc_bundle*1e3,
        P_peec_W_m=P_peec, P_intra_W_m=P_intra, P_total_W_m=P_peec+P_intra,
        eta_c_percent=eta_c, eta_m_percent=eta_m,
        sbar_mean=float(np.mean(sbar)), sbar_max=float(np.max(sbar)),
        center_current_mA=float(abs(Iz[ic])*1e3),
        outer_loss_share_percent=float(100*np.sum(P_ring[isOuter])/np.sum(P_ring)),
        min_center_dist_um=float(np.sqrt(mind2)*1e6),
        Zint_ratio=float(np.real(Zint)/Rdc_strand), N=int(N), K=K)

def layerq(u, Fv):
    lay=np.zeros(u.shape, dtype=int)
    for j in range(len(Fv)-1):
        lay += (u>Fv[j])
    return lay

# frozen regression anchors
BASE=dict(ringmax=10, p_mm=0.15, d_mm=0.14, sigma=5.8e7, f=200e3, I=20,
          K=16, Lambda_mm=40, alpha_deg=25, pitch_mm=40)

def cfg(**kw):
    c=dict(BASE); c.update(kw); return c

if __name__=='__main__':
    print("== frozen regression anchors (verify vs protocol_frozen regression_anchor) ==")
    rS=q3_solver(cfg(scheme='straight', intra=True))
    rT=q3_solver(cfg(scheme='twist',    intra=True))
    print(f"straight Rac={rS['Rac_peec_mohm']:.12f} (anchor 14.357229796331)  dR={rS['Rac_peec_mohm']-14.357229796331:+.2e}")
    print(f"straight J   ={rS['J_peec']:.9f}      (anchor 4.242996)")
    print(f"straight Zint={rS['Zint_ratio']:.12f} (anchor 1.001048146746)")
    print(f"straight P_intra={rS['P_intra_W_m']:.6f} (expect ~0.4016 ; Q2 COMSOL 0.395802)")
    print(f"twist    Rac={rT['Rac_peec_mohm']:.12f} (anchor 14.0359392794269) dR={rT['Rac_peec_mohm']-14.0359392794269:+.2e}")
    print(f"twist    J   ={rT['J_peec']:.9f}      (anchor 4.086757)")
    print(f"twist    sbar={rT['sbar_mean']:.12f} (anchor 1.015070108124)")
    print("\n== G1 (straight / twist) full metrics ==")
    for nm,r in [('straight',rS),('twist',rT)]:
        print(f" {nm}: J_peec={r['J_peec']:.6f} J_full={r['J_full']:.6f} Rac={r['Rac_peec_mohm']:.6f} "
              f"Rdc={r['Rdc_mohm']:.6f} P_peec={r['P_peec_W_m']:.6f} P_intra={r['P_intra_W_m']:.6f} "
              f"P_tot={r['P_total_W_m']:.6f} etaC={r['eta_c_percent']:.6f} etaM={r['eta_m_percent']:.6f} "
              f"center={r['center_current_mA']:.6f} outer={r['outer_loss_share_percent']:.6f} "
              f"sbar={r['sbar_mean']:.6f} minUm={r['min_center_dist_um']:.2f}")

    print("\n== K sweep braid(alpha=25,Lambda=40) ==")
    prev=None
    for K in [8,16,32,64,128,256]:
        r=q3_solver(cfg(scheme='braid', K=K, alpha_deg=25, Lambda_mm=40, intra=True))
        adj = np.nan if prev is None else (r['Rac_peec_mohm']-prev)/prev*100
        print(f" braid K={K:3d}(eff {r['K']}): Rac={r['Rac_peec_mohm']:.6f} J_peec={r['J_peec']:.6f} "
              f"J_full={r['J_full']:.6f} P={r['P_total_W_m']:.6f} etaI(mag)={r['eta_m_percent']:.4f} "
              f"etaI(cplx)={r['eta_c_percent']:.4f} clamp=0 minUm={r['min_center_dist_um']:.1f} "
              f"sbar={r['sbar_mean']:.6f} center={r['center_current_mA']:.3f} outer={r['outer_loss_share_percent']:.4f} "
              f"dRac_adj={'n/a' if np.isnan(adj) else format(adj,'+.6f')+'%'}")
        prev=r['Rac_peec_mohm']
