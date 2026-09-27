function res = q4_peec_core(X,Y,z,a,sigma,f,Itot)
%Q4_PEEC_CORE  PEEC common-voltage solve for arbitrary periodic centerlines.
% X,Y,z : nz x N arrays (m), periodic trajectories sampled at uniform z.
% Uses the same formulas as q3_solver.m (Bessel internal Z, ln mutual, s_bar,
% transport constraint sum I/s = Itot).
mu=4*pi*1e-7; w=2*pi*f; nz=size(X,1); N=size(X,2);
dz=z(2)-z(1);
VX=(X([2:nz,1],:)-X([nz,1:nz-1],:))/(2*dz);
VY=(Y([2:nz,1],:)-Y([nz,1:nz-1],:))/(2*dz);
S=sqrt(1+VX.^2+VY.^2);
delta=sqrt(2/(w*mu*sigma)); kc=(1-1i)/delta;
Zint=kc/(2*pi*a*sigma)*besselj(0,kc*a)/besselj(1,kc*a);
Macc=zeros(N); sacc=zeros(N,1); mind2=Inf;
for k=1:nz
    xk=X(k,:).'; yk=Y(k,:).'; vx=VX(k,:).'; vy=VY(k,:).'; sk=S(k,:).';
    dx=xk-xk'; dy=yk-yk'; d2=dx.^2+dy.^2;
    dt=d2; dt(1:N+1:end)=Inf; mind2=min(mind2,min(dt(:)));
    d2(1:N+1:end)=a^2;
    Macc=Macc+((1+vx*vx'+vy*vy')./(sk*sk')).*(-0.5*log(d2));
    sacc=sacc+sk;
end
sbar=sacc/nz; Mavg=Macc/nz;
Zmat=diag(Zint*sbar)+1i*(w*mu/(2*pi))*Mavg;
wi=1./sbar; zsol=Zmat\wi; V=Itot/(wi'*zsol); Iwire=V*zsol; Iz=Iwire./sbar;
P=sum(real(Zint)*sbar.*abs(Iwire).^2); Rac=P/Itot^2;
As=pi*a^2; Rdc=(1/(sigma*As))/sum(1./sbar);
Ibar=Itot/N;
res=struct('J',Rac/Rdc,'Rac',Rac,'Rdc',Rdc,'P',P,'sbar',mean(sbar), ...
    'eta',max(abs(Iz-Ibar))/Ibar*100,'mind',sqrt(mind2),'Iz',Iz,'S',S);
end
