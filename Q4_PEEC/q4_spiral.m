function res = q4_spiral(N, R, rin_frac, W, P, K, A_Cu)
%Q4_SPIRAL  Spiral full-exchange: stations on an area-quantile Archimedean spiral,
% driven by an EXACT cyclic shift. Satisfies simultaneously:
%   (1) cyclic-shift orbit (=> eta_I -> 0),
%   (2) area-quantile radial dwell (=> D_O -> 0),
%   (3) checkable spacing / bend feasibility.
%   loop:  r(phi)= sqrt(rin^2 + tri(phi/N)(R^2-rin^2)),  theta(phi)=2*pi*W*phi/N
%   strand k at z:  phi_k = mod(k + N z/P, N),  position = C(phi_k)
sigma=5.8e7; f=200e3; Itot=20; mu=4*pi*1e-7;
a=0.5*sqrt(4*A_Cu/(pi*N)); rin=rin_frac*R;
z=linspace(0,P,K+1); z(end)=[]; dz=P/K; speed=N/P; k=(0:N-1)';
tri=@(u) 2*min(u,1-u);
X=zeros(N,K); Y=zeros(N,K);
for s=1:K
    phi=mod(k+speed*z(s),N); u=phi/N;
    r=sqrt(rin^2+tri(u)*(R^2-rin^2)); th=2*pi*W*u;
    X(:,s)=r.*cos(th); Y(:,s)=r.*sin(th);
end
p=q4_peec_core(X,Y,z,a,sigma,f,Itot);

% min spacing over a few generic sections
mind=Inf;
for d0=[0.15 0.5 0.85]
    ph=mod(k+d0,N); u=ph/N; r=sqrt(rin^2+tri(u)*(R^2-rin^2)); th=2*pi*W*u;
    x0=r.*cos(th); y0=r.*sin(th);
    dd=(x0-x0').^2+(y0-y0').^2; dd(1:N+1:end)=Inf; mind=min(mind,min(dd(:)));
end
mind=sqrt(mind);

% 3D bend radius
VX=(X(:,[2:K,1])-X(:,[K,1:K-1]))/(2*dz); VY=(Y(:,[2:K,1])-Y(:,[K,1:K-1]))/(2*dz);
AX=(VX(:,[2:K,1])-VX(:,[K,1:K-1]))/(2*dz); AY=(VY(:,[2:K,1])-VY(:,[K,1:K-1]))/(2*dz);
k3=sqrt(AY.^2+AX.^2+(VX.*AY-VY.*AX).^2)./(1+VX.^2+VY.^2).^1.5;
Rb=1/max(k3(:));

% D_O vs equal-area bands
M=8; edges=linspace(rin,R,M+1); edges(1)=edges(1)-1e-12; edges(end)=edges(end)+1e-12;
rr=hypot(X,Y); bidx=discretize(rr,edges);
target=(edges(2:end).^2-edges(1:end-1).^2)/(R^2-rin^2);
dwell=zeros(N,M); for m=1:M, dwell(:,m)=mean(bidx==m,2); end
D_O=max(sum(abs(dwell-target),2)); D_O_inter=max(max(dwell)-min(dwell));

res=struct('N',N,'R_mm',R*1e3,'rin_frac',rin_frac,'W',W,'d_um',2*a*1e6, ...
    'J',p.J,'Rac',p.Rac,'eta',p.eta,'sbar',p.sbar,'mind_um',mind*1e6, ...
    'min_bend_d',Rb/(2*a),'D_O',D_O,'D_O_inter',D_O_inter,'X',X,'Y',Y,'z',z);
end
