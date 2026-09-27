function xyz = q3_ml_geom(factors, pitches_mm, A_mm2, gap_um, spp)
%Q3_ML_GEOM  Replicate litz_multilevel RecursiveCable geometry (numerical tangent),
% returns xyz (n x nt x 3), axial param t = column index (m).
if nargin<5, spp=128; end
factors=factors(:).'; pitches=pitches_mm(:).';
n=prod(factors); a=sqrt(A_mm2*1e-6/(pi*n)); gap=gap_um*1e-6;
offsets=cell(1,numel(factors)); env=a;
for L=1:numel(factors)
    xy=ml_layout(factors(L))*(2*env+gap);
    offsets{L}=xy; env=env+max(vecnorm(xy,2,2));
end
pu=round(abs(pitches)*1000); period=ml_lcm(pu)*1e-6;      % m
nt=max(129, ceil(period/(min(pu)*1e-6)*spp)+1);
t=linspace(0,period,nt);
c=zeros(1,nt,3); c(1,:,3)=t;
e1=zeros(1,nt,3); e1(1,:,1)=1;
e2=zeros(1,nt,3); e2(1,:,2)=1;
for level=numel(factors):-1:1
    xy=offsets{level}; r=vecnorm(xy,2,2); ph0=atan2(xy(:,2),xy(:,1));
    K=size(xy,1); M=size(c,1);
    phi=2*pi*t/(pitches(level)*1e-3) + ph0;               % K x nt
    cc=reshape(c,M,1,nt,3); e1b=reshape(e1,M,1,nt,3); e2b=reshape(e2,M,1,nt,3);
    uc=reshape(cos(phi),1,K,nt); us=reshape(sin(phi),1,K,nt);
    u=uc.*e1b + us.*e2b;                                  % M x K x nt x 3
    child=cc + reshape(r,1,K,1,1).*u;
    c=reshape(child,M*K,nt,3);
    u=reshape(u,M*K,nt,3);
    % tangent via periodic central difference of v = c - t*ez, plus ez
    v=c; v(:,:,3)=v(:,:,3)-t;
    dv=(v(:,[2:nt,1],:)-v(:,[nt,1:nt-1],:))/(2*(t(2)-t(1)));
    tg=dv; tg(:,:,3)=tg(:,:,3)+1; tg=tg./vecnorm(tg,2,3);
    uu=u - sum(u.*tg,3).*tg; uu=uu./vecnorm(uu,2,3);
    e1=uu; e2=cross(tg,e1,3);
end
xyz=c;
end

function xy=ml_layout(n)
if n==1, xy=[0 0]; return; end
if n==7
    a=(0:5)'*pi/3; xy=[[0 0];[cos(a) sin(a)]]; return;
end
if n>=2 && n<=6
    a=(0:n-1)'*2*pi/n; xy=[cos(a) sin(a)]/(2*sin(pi/n)); return;
end
k=ceil(sqrt(n)); P=[];
for q=-k:k
    for r=-k:k
        P(end+1,:)=[q+r/2, sqrt(3)/2*r]; %#ok<AGROW>
    end
end
r2=round(sum(P.^2,2),10); ang=atan2(P(:,2),P(:,1));
[~,idx]=sortrows([r2 ang],[1 2]);
xy=P(idx(1:n),:); xy=xy-mean(xy,1);
end

function L=ml_lcm(v)
L=1; for x=v(:)', L=L/gcd(L,x)*x; end
end
