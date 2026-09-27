function q3_ingest_litz()
%Q3_INGEST_LITZ  Run the PEEC core on COMSOL (litz_q3) centerlines and compare.
% Reads litz_<case>_centerlines.csv (from litz_npz_to_csv.py), builds the mutual
% matrix from the actual 3D geometry, solves the parallel common-voltage system
% with the same formulas as q3_solver.m, and compares K=Rac/Rdc_geom to COMSOL.
% PEEC omits intra-strand proximity, so well-transposed cases are expected to
% sit below COMSOL K.
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end

sigma=5.8e7; mu=4*pi*1e-7; f=200e3; w=2*pi*f; Itot=20; A_Cu=6.0e-6;
delta=sqrt(2/(w*mu*sigma));

cases={ ...
 'rigid_twist',    4.7604979021993294, 13.696118621484092e-3; ...
 'shell_exchange', 4.663413488529197,  13.420472865211806e-3; ...
 'full_exchange',  1.5961987792728285, 4.918691463531148e-3; ...
 'final_484',      1.4971216482038827, 4.397906322047884e-3};

fprintf('========== PEEC on COMSOL (litz_q3) geometry ==========\n');
rows={};
for c=1:size(cases,1)
    tag=cases{c,1}; Kcomsol=cases{c,2}; RacComsol=cases{c,3};
    r=peec_on_case(tag,A_Cu,sigma,mu,w,delta,Itot);
    fprintf('%-15s N=%d K_peec=%.5f K_comsol=%.5f | Rac=%.4f vs %.4f mOhm | diff=%.2f%% eta=%.1f%% sbar=%.5f min=%.1fum\n', ...
        tag,r.N,r.J,Kcomsol,r.Rac*1e3,RacComsol*1e3,100*(r.J-Kcomsol)/Kcomsol,r.eta,r.sbar,r.mind*1e6);
    rows(end+1,:)={tag,r.N,r.J,Kcomsol,r.Rac*1e3,RacComsol*1e3,100*(r.J-Kcomsol)/Kcomsol,r.eta,r.sbar,r.mind*1e6}; %#ok<AGROW>
end
TT=cell2table(rows,'VariableNames',{'case','N','K_peec','K_comsol','Rac_peec_mohm','Rac_comsol_mohm','K_diff_pct','eta_pct','sbar','min_dist_um'});
writetable(TT,'litz_peec_vs_comsol.csv');
fprintf('\nSaved litz_peec_vs_comsol.csv\n'); disp(TT);

% per-strand current check for final_484 (COMSOL verified, P=0.5, N=484)
r=peec_on_case('final_484',A_Cu,sigma,mu,w,delta,Itot);
C=readmatrix('litz_final_484_currents.csv'); Icomsol=C(:,2)+1i*C(:,3);
Ibar=Itot/r.N;
g=(r.Iz'*Icomsol)/abs(r.Iz'*Icomsol); Iz_al=g*r.Iz;
fprintf('\n[final_484] per-strand COMSOL vs PEEC: |I| max dev=%.4f%%, complex max dev=%.4f%% (of Ibar), corr=%.8f\n', ...
    max(abs(abs(r.Iz)-abs(Icomsol)))/Ibar*100, max(abs(Iz_al-Icomsol))/Ibar*100, ...
    abs(r.Iz'*Icomsol)/(norm(r.Iz)*norm(Icomsol)));
fprintf('[final_484] COMSOL |I| min/mean/max = %.5f/%.5f/%.5f A (20/N=%.5f)\n', ...
    min(abs(Icomsol)),mean(abs(Icomsol)),max(abs(Icomsol)),Ibar);
end

function r=peec_on_case(tag,A_Cu,sigma,mu,w,delta,Itot)
T=readmatrix(sprintf('litz_%s_centerlines.csv',tag));
st=T(:,1); sd=T(:,2); x=T(:,3); y=T(:,4); z=T(:,5);
nz=max(st)+1; N=max(sd)+1;
X=reshape(x,N,nz).'; Y=reshape(y,N,nz).'; Z=reshape(z,N,nz).';   % nz x N
dz=Z(2,1)-Z(1,1);
a=sqrt(A_Cu/(pi*N)); As=pi*a^2;
kc=(1-1i)/delta; Zint=kc/(2*pi*a*sigma)*besselj(0,kc*a)/besselj(1,kc*a);

VX=zeros(nz,N); VY=zeros(nz,N);
VX(2:nz-1,:)=(X(3:nz,:)-X(1:nz-2,:))/(2*dz);
VY(2:nz-1,:)=(Y(3:nz,:)-Y(1:nz-2,:))/(2*dz);
VX(1,:)=VX(2,:); VX(nz,:)=VX(nz-1,:);
VY(1,:)=VY(2,:); VY(nz,:)=VY(nz-1,:);
S=sqrt(1+VX.^2+VY.^2);

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
Rdc=(1/(sigma*As))/sum(1./sbar);
Ibar=Itot/N;
r=struct('N',N,'J',Rac/Rdc,'Rac',Rac,'Rdc',Rdc,'P',P,'sbar',mean(sbar), ...
    'eta',max(abs(Iz-Ibar))/Ibar*100,'mind',sqrt(mind2),'Iz',Iz);
end
