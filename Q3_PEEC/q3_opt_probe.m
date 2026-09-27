function q3_opt_probe()
%Q3_OPT_PROBE  Is there head-room left? Extend Lambda, check alpha floor,
% and compare against the model's own theoretical floor (isolated-strand skin).
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end
base=struct('ringmax',10,'p_mm',0.15,'d_mm',0.14,'sigma',5.8e7,'f',200e3,'I',20,'K',256);

% model floor: perfect equalization + isolated-strand skin factor
a=0.5*base.d_mm*1e-3; w=2*pi*base.f;
sigma=base.sigma; mu=4*pi*1e-7;
delta=sqrt(2/(w*mu*sigma)); kc=(1-1i)/delta;
Zint=kc/(2*pi*a*sigma)*besselj(0,kc*a)/besselj(1,kc*a);
floorJ=real(Zint)/(1/(sigma*pi*a^2));
fprintf('model floor J (perfect equalization + isolated skin) = %.6f\n', floorJ);

alphas=[2 5 10 15];
Lams=[240 400 600 1000];
fprintf('\nbraid: J vs Lambda (K=256)\n%6s','alpha');
fprintf('%12d',Lams); fprintf('\n');
for al=alphas
    fprintf('%6d',al);
    for L=Lams
        c=base; c.scheme='braid'; c.alpha_deg=al; c.Lambda_mm=L;
        r=q3_solver(c);
        fprintf('%12.5f',r.J_peec);
    end
    fprintf('\n');
end

% physical-range candidate: alpha=10, longer Lambda
fprintf('\nphysical-range candidate (alpha=10 deg):\n');
fprintf('%8s %10s %10s %8s %8s\n','Lambda','Rac(mOhm)','J','sbar','eta%');
for L=[40 80 120 240]
    c=base; c.scheme='braid'; c.alpha_deg=10; c.Lambda_mm=L;
    r=q3_solver(c);
    fprintf('%8d %10.4f %10.5f %8.5f %8.1f\n',L,r.Rac_peec_mohm,r.J_peec,r.sbar_mean,r.eta_c_percent);
end
end
