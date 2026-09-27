function res = q4_topology_teammate(side, half_width, P, K, A_Cu, w)
%Q4_TOPOLOGY_TEAMMATE  Teammate's full_exchange topology (even square grid +
% periodic-spline Hamiltonian cycle + cyclic phase schedule) with certificates.
%   centerlines: pos_k(z) = spline((k + n z/P) mod n), n = side^2
% optional w>0 : periodic moving-average smoothing of the loop (corner rounding).
% certificates: D_O over equal-area radial bands, min spacing, bend radius, sbar.
% A_Cu is the total normal copper area (m^2); d = sqrt(4 A_Cu/(pi n)).
if nargin<6, w=0; end
sigma=5.8e7; f=200e3; Itot=20; mu=4*pi*1e-7;

n=side^2; a=0.5*sqrt(4*A_Cu/(pi*n));
ids=ham_cycle(side);
x=linspace(-half_width,half_width,side); [XX,YY]=meshgrid(x,x);
pos=[XX(:),YY(:)];
Yc=pos(ids,:);                         % closed loop samples (n x 2)
if w>0, Yc=circ_movavg(Yc,w); end      % corner rounding
posC=Yc;                               % station set after smoothing
[Mc,~]=per_spline_coef(Yc);            % second-derivative coefficients

z=linspace(0,P,K+1); z(end)=[]; dz=P/K; speed=n/P;
X=zeros(n,K); Y=zeros(n,K); VX=zeros(n,K); VY=zeros(n,K); AX=zeros(n,K); AY=zeros(n,K);
kidx=(0:n-1)';
for s=1:K
    phi=mod(kidx+speed*z(s),n);
    [S,dS,d2S]=per_spline_eval(Yc,Mc,phi);
    X(:,s)=S(:,1); Y(:,s)=S(:,2);
    VX(:,s)=dS(:,1)*speed; VY(:,s)=dS(:,2)*speed;
    AX(:,s)=d2S(:,1)*speed^2; AY(:,s)=d2S(:,2)*speed^2;
end

% ---- PEEC ----
p=q4_peec_core(X,Y,z,a,sigma,f,Itot);

% ---- feasibility: min spacing (one generic section) + 3D wire bend radius ----
phi0=mod(kidx+0.5,n); [S0,~,~]=per_spline_eval(Yc,Mc,phi0);
ddx=S0(:,1)-S0(:,1)'; ddy=S0(:,2)-S0(:,2)'; dd=ddx.^2+ddy.^2; dd(1:n+1:end)=Inf; mind=sqrt(min(dd(:)));
AX=(VX(:,[2:K,1])-VX(:,[K,1:K-1]))/(2*dz); AY=(VY(:,[2:K,1])-VY(:,[K,1:K-1]))/(2*dz);
k3=sqrt(AY.^2+AX.^2+(VX.*AY-VY.*AX).^2)./(1+VX.^2+VY.^2).^1.5;
Rb=1/max(k3(:));                       % min 3D bend radius (includes z stretch)

% ---- D_O certificate (continuous spline dwell vs equal-count station bands, §2.3) ----
r=hypot(X,Y); rr=hypot(posC(:,1),posC(:,2));
M=8; rs=sort(rr); ei=round(linspace(1,n,M+1)); e=rs(ei);
for i=2:numel(e), if e(i)<=e(i-1), e(i)=e(i-1)+1e-12; end, end
e(1)=e(1)-1e-12; e(end)=e(end)+1e-12;
bidx=discretize(r,e);
target=zeros(1,M);
for m=1:M, target(m)=sum(rr>=e(m)&rr<e(m+1))/n; end
dwell=zeros(n,M);
for m=1:M, dwell(:,m)=mean(bidx==m,2); end
D_O=max(sum(abs(dwell-target),2));       % continuous (uniform-speed) area defect
D_O_inter=max(max(dwell)-min(dwell));    % inter-strand spread (0 => eta_I=0)

res=struct('side',side,'N',n,'d_um',2*a*1e6,'P_mm',P*1e3,'w',w,'J',p.J,'Rac',p.Rac, ...
    'P_W',p.P,'eta',p.eta,'sbar',p.sbar,'mind_um',mind*1e6,'min_bend_d',Rb/(2*a), ...
    'max_angle_deg',atan(max(sqrt(p.S(:).^2-1)))*180/pi,'D_O',D_O,'D_O_inter',D_O_inter, ...
    'X',X,'Y',Y,'z',z);
end

function Z=circ_movavg(Y,w)
n=size(Y,1); k=ones(2*w+1,1)/(2*w+1);
Ye=[Y(end-w+1:end,:);Y;Y(1:w,:)];
Z=conv2(Ye,k,'same'); Z=Z(w+1:w+n,:);
end

% ---------------- helpers ----------------
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
%Periodic cubic spline second derivatives for unit-spaced closed data Y (n x 2).
n=size(Y,1);
rhs=6*([Y(2:end,:);Y(1,:)]-2*Y+[Y(end,:);Y(1:end-1,:)]);
lam=4+2*cos(2*pi*(0:n-1)'/n);           % eigenvalues of circulant [1,4,1]
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

function q=qntl(v,p)
v=sort(v(:)); n=numel(v); idx=1+(n-1)*p;
lo=floor(idx); hi=ceil(idx); w=idx-lo;
q=v(lo).*(1-w)+v(hi).*w;
end
