function q3_ml_check()
%Q3_ML_CHECK  Cross-check litz_multilevel (N=64 multi-level stranding) with PEEC.
% Replicates their recursive geometry and runs our PEEC core; compares Rac/K to
% their published COMSOL values.
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end
A=6.0; gap=40; sigma=5.8e7; f=200e3; Itot=20;
a=0.5*sqrt(4*A*1e-6/(pi*64)); d_um=2*a*1e6;
fprintf('N=64, d=%.4f um, beta=d/delta=%.3f\n', d_um, 2*a/sqrt(2/(2*pi*f*4*pi*1e-7*sigma)));

cases={
  'baseline single-64 P=96',   [64],   [96],     4.967561407099009, 14.314942597278368e-3;
  'recommended 4x16 -16/+128', [4 16], [-16 128], 4.547793,          13.155524e-3};

fprintf('%-28s %8s %9s %9s %9s %9s %8s\n','case','PEEC K','PEEC Rac','COMSOL K','COMSOL Rac','K diff%','mind_um');
for c=1:size(cases,1)
    xyz=q3_ml_geom(cases{c,2},cases{c,3},A,gap,128);
    [n,nt,~]=size(xyz);
    pu=round(abs(cases{c,3})*1000); L=1; for x=pu, L=L/gcd(L,x)*x; end
    period=L*1e-6; z=linspace(0,period,nt);
    X=squeeze(xyz(:,:,1)).'; Y=squeeze(xyz(:,:,2)).';
    r=q4_peec_core(X,Y,z,a,sigma,f,Itot);
    Kc=cases{c,4}; Racc=cases{c,5};
    fprintf('%-28s %8.4f %9.4f %9.4f %9.4f %8.2f %8.1f\n', cases{c,1}, r.J, r.Rac*1e3, Kc, Racc*1e3, 100*(r.J-Kc)/Kc, r.mind*1e6);
end
end
