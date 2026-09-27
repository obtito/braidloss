function res = q3_solver(cfg)
%Q3_SOLVER  Lean Q3 PEEC evaluator, four schemes.
%   cfg.scheme = 'straight' | 'twist' | 'braid' | 'counter_braid'
%   Layout: Q2 circular rings (capacities 1,6,...,60; r_k = p_mm*k).
%   Radial transposition: per-strand triangle path with a golden-ratio phase.
%   Per-strand phase (instead of ring-locked phase) keeps strands apart when a
%   ring moves inward, so the geometry stays physical and no distance clamp is
%   ever applied: min strand spacing is ASSERTED to be >= 2a.
%   Parallel common-voltage PEEC; straight/twist reproduce the Q2 PEEC anchors.
%
%   cfg: scheme, ringmax, p_mm, d_mm, sigma, f, I, K,
%        Lambda_mm, alpha_deg (braid/counter), pitch_mm (twist),
%        intra (optional, default false: add analytic intra-strand proximity)
%
%   res: J_peec, J_full, Rac_peec_mohm, Rac_full_mohm, Rdc_mohm, P_peec_W_m,
%        P_intra_W_m, P_total_W_m, eta_c_percent, eta_m_percent, sbar_mean,
%        center_current_mA, outer_loss_share_percent, min_center_dist_um,
%        Zint_ratio, I_z, ring, r0

sigma=cfg.sigma; mu=4*pi*1e-7; f=cfg.f; w=2*pi*f; Itot=cfg.I;
a=0.5*cfg.d_mm*1e-3; As=pi*a^2;
delta=sqrt(2/(w*mu*sigma)); kc=(1-1i)/delta;
Zint=kc/(2*pi*a*sigma)*besselj(0,kc*a)/besselj(1,kc*a);   % ohm/m
Rdc_strand=1/(sigma*As);
if ~isfield(cfg,'intra'), cfg.intra=false; end

% ---- Q2 circular layout ----
K0=cfg.ringmax; cap=6*(0:K0); cap(1)=1;
ring=zeros(sum(cap),1); th0=zeros(sum(cap),1); s=0;
for k=0:K0
    j=(0:cap(k+1)-1)';
    ring(s+1:s+numel(j))=k; th0(s+1:s+numel(j))=2*pi*j/cap(k+1); s=s+numel(j);
end
N=numel(ring); p=cfg.p_mm*1e-3; rr=p*(0:K0)';
r0=p*ring; x=r0.*cos(th0); y=r0.*sin(th0);

% ---- trajectory ----
switch cfg.scheme
    case 'straight'
        K=1; X=x; Y=y; S=ones(N,1); VX=zeros(N,1); VY=zeros(N,1);
    case 'twist'
        K=1; om=2*pi/(cfg.pitch_mm*1e-3);
        X=x; Y=y; VX=-om*y; VY=om*x; S=sqrt(1+(om*r0).^2);
    case {'braid','counter_braid'}
        K=max(4,round(cfg.K)); Lam=cfg.Lambda_mm*1e-3;
        tanA=tan(cfg.alpha_deg*pi/180);
        beta=tanA/max(rr);                              % rigid angular rate
        Fv=cumsum(cap)/N;                               % layer CDF
        phi=((1:N)'-0.5)/N;                             % uniform per-strand phase
        if strcmp(cfg.scheme,'counter_braid'), fam=2*mod(ring,2)-1; else, fam=ones(N,1); end
        zm=((0:K-1)+0.5)*(Lam/K);
        X=zeros(N,K); Y=zeros(N,K); VX=zeros(N,K); VY=zeros(N,K); S=zeros(N,K);
        drdzmag=(rr(end)-rr(1))/(Lam/2);                % triangle slope magnitude
        for m=1:K
            u=mod(phi+zm(m)/Lam,1);
            tri=2*min(u,1-u);                           % rising then falling, one cycle/Lam
            lay=layerq(tri,Fv);                         % occupancy-conserving layer
            % assert occupancy equals ring capacities
            cnt=accumarray(lay+1,1,[K0+1,1])';
            if any(cnt~=cap)
                error('q3_solver:occupancy','layer counts [%s] != caps at section %d', num2str(cnt), m);
            end
            r=p*lay;
            th=th0+fam*beta*zm(m);
            for k=0:K0
                idx=find(lay==k);
                [~,o]=sort(th(idx));
                th(idx(o))=2*pi*(0:numel(idx)-1)'/cap(k+1);   % music-chairs slots
            end
            drdz=drdzmag*(1-2*(u>=0.5));                % signed triangle slope
            dthdz=fam*beta;
            X(:,m)=r.*cos(th); Y(:,m)=r.*sin(th);
            VX(:,m)=drdz.*cos(th)-r.*sin(th).*dthdz;
            VY(:,m)=drdz.*sin(th)+r.*cos(th).*dthdz;
            S(:,m)=sqrt(1+drdz.^2+(r.*dthdz).^2);
        end
    otherwise
        error('q3_solver:scheme','unknown scheme %s',cfg.scheme);
end

% ---- mutual matrix, per-section averaging (no distance clamp: assert instead) ----
Macc=zeros(N); sacc=zeros(N,1); mind2=Inf;
for kk=1:K
    if K==1
        xk=X; yk=Y; vx=VX; vy=VY; sk=S;
    else
        xk=X(:,kk); yk=Y(:,kk); vx=VX(:,kk); vy=VY(:,kk); sk=S(:,kk);
    end
    dx=xk-xk'; dy=yk-yk'; d2=dx.^2+dy.^2;
    dtmp=d2; dtmp(1:N+1:end)=Inf; mind2=min(mind2,min(dtmp(:)));
    if mind2 < (2*a)^2
        error('q3_solver:spacing', ...
            'strand spacing %.1f um < 2a=%.1f um at section %d (scheme %s)', ...
            sqrt(mind2)*1e6, 2*a*1e6, kk, cfg.scheme);
    end
    d2(1:N+1:end)=a^2;
    Mseg=-0.5*log(d2);
    tdot=(1+vx*vx'+vy*vy')./(sk*sk');
    Macc=Macc+tdot.*Mseg; sacc=sacc+sk;
end
sbar=sacc/K; Mavg=Macc/K;

% ---- parallel common-voltage solve (Q2 convention) ----
Zmat=diag(Zint*sbar)+1i*(w*mu/(2*pi))*Mavg;
wi=1./sbar;
zsol=Zmat\wi;
Vdrive=Itot/(wi'*zsol);
Iwire=Vdrive*zsol;
Iz=Iwire./sbar;

Ibar=Itot/N;
eta_c=max(abs(Iz-Ibar))/Ibar*100;
eta_m=max(abs(abs(Iz)-Ibar))/Ibar*100;

P_peec=sum(real(Zint)*sbar.*abs(Iwire).^2);
Rac_peec=P_peec/Itot^2;
Rdc_bundle=(1/(sigma*As))/sum(1./sbar);
J_peec=Rac_peec/Rdc_bundle;

% ---- optional analytic intra-strand proximity term (uniform-field, RMS) ----
P_intra=0;
if cfg.intra
    Cintra=sigma*pi*a^4*w^2/4;      % W/m per strand per T_rms^2
    B2=zeros(N,1);
    for kk=1:K
        if K==1, xk=X; yk=Y; else, xk=X(:,kk); yk=Y(:,kk); end
        dx=xk-xk'; dy=yk-yk'; d2=max(dx.^2+dy.^2,(2*a)^2); d2(1:N+1:end)=Inf;
        B=(mu/(2*pi))*(((-dy+1i*dx)./d2)*Iwire);
        B2=B2+abs(B).^2;
    end
    P_intra=Cintra*sum(B2/K);
end
Rac_full=(P_peec+P_intra)/Itot^2;
J_full=Rac_full/Rdc_bundle;

% ---- auxiliary metrics ----
isOuter=ring==K0;
P_ring=real(Zint)*sbar.*abs(Iwire).^2;
ic=find(ring==0,1); if isempty(ic), ic=1; end

res.scheme=cfg.scheme;
res.J_peec=J_peec; res.J_full=J_full;
res.Rac_peec_mohm=Rac_peec*1e3; res.Rac_full_mohm=Rac_full*1e3; res.Rdc_mohm=Rdc_bundle*1e3;
res.P_peec_W_m=P_peec; res.P_intra_W_m=P_intra; res.P_total_W_m=P_peec+P_intra;
res.eta_c_percent=eta_c; res.eta_m_percent=eta_m;
res.sbar_mean=mean(sbar); res.sbar_max=max(sbar);
res.center_current_mA=abs(Iz(ic))*1e3;
res.outer_loss_share_percent=100*sum(P_ring(isOuter))/sum(P_ring);
res.min_center_dist_um=sqrt(mind2)*1e6;
res.Zint_ratio=real(Zint)/Rdc_strand;
res.N=N; res.I_z=Iz; res.ring=ring; res.r0=r0;
% geometry (for visualization): positions and length factors per section
res.X=X; res.Y=Y; res.S=S; res.K=K;
if isfield(cfg,'alpha_deg'), res.alpha_deg=cfg.alpha_deg; end
if isfield(cfg,'Lambda_mm'), res.Lambda_mm=cfg.Lambda_mm; end
end

function lay = layerq(u, Fv)
%LAYERQ  Quantile layer index: lay = #{CDF boundaries < u}, boundaries Fv.
u=u(:); lay=zeros(numel(u),1);
for j=1:numel(Fv)-1
    lay=lay+(u>Fv(j));
end
end
