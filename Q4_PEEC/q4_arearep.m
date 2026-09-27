function res = q4_arearep(side, half_width, rin_frac, P, K, A_Cu)
%Q4_AREAREP  Area-reparameterised transposition on the teammate square grid.
%   Angular schedule theta_k(z) from the teammate Hamilton loop;
%   radial profile replaced by the AREA QUANTILE of an annulus [rin,R]:
%       r_k(z) = sqrt(rin^2 + tri^2 (R^2-rin^2)),  tri = triangle(mod(k+n z/P)/n)
%   => dwell per radial band = band area  => D_O -> 0 (area fairness).
%   rin_frac = rin/R controls hollow packing (0 = filled).
sigma=5.8e7; f=200e3; Itot=20;
n=side^2; a=0.5*sqrt(4*A_Cu/(pi*n));
ids=ham_cycle(side);
x=linspace(-half_width,half_width,side); [XX,YY]=meshgrid(x,x);
pos=[XX(:),YY(:)]; Yc=pos(ids,:); [Mc,~]=per_spline_coef(Yc);
R=max(hypot(Yc(:,1),Yc(:,2))); rin=rin_frac*R;

z=linspace(0,P,K+1); z(end)=[]; dz=P/K; speed=n/P; kidx=(0:n-1)';
X=zeros(n,K); Y=zeros(n,K);
for s=1:K
    phi=mod(kidx+speed*z(s),n);
    S=per_spline_eval(Yc,Mc,phi);
    th=atan2(S(:,2),S(:,1));
    u=phi/n; tri=2*min(u,1-u);
    r=sqrt(rin^2+tri*(R^2-rin^2));
    X(:,s)=r.*cos(th); Y(:,s)=r.*sin(th);
end

p=q4_peec_core(X,Y,z,a,sigma,f,Itot);

% feasibility
s0=round(K/3); x0=X(:,s0); y0=Y(:,s0);
ddx=x0-x0'; ddy=y0-y0'; dd=ddx.^2+ddy.^2; dd(1:n+1:end)=Inf; mind=sqrt(min(dd(:)));
VX=(X(:,[2:K,1])-X(:,[K,1:K-1]))/(2*dz); VY=(Y(:,[2:K,1])-Y(:,[K,1:K-1]))/(2*dz);
AX=(VX(:,[2:K,1])-VX(:,[K,1:K-1]))/(2*dz); AY=(VY(:,[2:K,1])-VY(:,[K,1:K-1]))/(2*dz);
k3=sqrt(AY.^2+AX.^2+(VX.*AY-VY.*AX).^2)./(1+VX.^2+VY.^2).^1.5;
Rb=1/max(k3(:));

% D_O vs equal-area radial bands
M=8; edges=linspace(rin,R,M+1); edges(1)=edges(1)-1e-12; edges(end)=edges(end)+1e-12;
rr=hypot(X,Y); bidx=discretize(rr,edges);
target=(edges(2:end).^2-edges(1:end-1).^2)/(R^2-rin^2);
dwell=zeros(n,M); for m=1:M, dwell(:,m)=mean(bidx==m,2); end
D_O=max(sum(abs(dwell-target),2));
D_O_inter=max(max(dwell)-min(dwell));

res=struct('side',side,'N',n,'d_um',2*a*1e6,'P_mm',P*1e3,'rin_frac',rin_frac, ...
    'J',p.J,'Rac',p.Rac,'eta',p.eta,'sbar',p.sbar,'mind_um',mind*1e6, ...
    'min_bend_d',Rb/(2*a),'D_O',D_O,'D_O_inter',D_O_inter,'X',X,'Y',Y,'z',z);
end

function ids = ham_cycle(side)
xy=[];
for y=0:side-1, xy(end+1,:)=[0,y]; end %#ok<AGROW>
for xx=1:side-1
    if mod(xx,2)==1, ys=side-1:-1:1; else, ys=1:side-1; end
    for y=ys, xy(end+1,:)=[xx,y]; end %#ok<AGROW>
end
for xx=side-1:-1:1, xy(end+1,:)=[xx,0]; end %#ok<AGROW>
ids=xy(:,2)*side+xy(:,1)+1;
end

function [M,rhs]=per_spline_coef(Y)
n=size(Y,1);
rhs=6*([Y(2:end,:);Y(1,:)]-2*Y+[Y(end,:);Y(1:end-1,:)]);
lam=4+2*cos(2*pi*(0:n-1)'/n);
M=real(ifft(fft(rhs)./lam));
end

function [S,dS,d2S]=per_spline_eval(Y,M,phi)
n=size(Y,1);
i0=mod(floor(phi),n)+1; u=phi-floor(phi); i1=mod(i0,n)+1;
Y0=Y(i0,:); Y1=Y(i1,:); M0=M(i0,:); M1=M(i1,:);
S=Y0.*(1-u)+Y1.*u + (((1-u).^3-(1-u)).*M0 + (u.^3-u).*M1)/6;
dS=(Y1-Y0) + ((-3*(1-u).^2+1).*M0 + (3*u.^2-1).*M1)/6;
d2S=(1-u).*M0 + u.*M1;
end
