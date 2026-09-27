%% q3_solver_selftest.m — lean solver self-test (must pass before any sweep)
% S1 straight : reproduce Q2 PEEC circular anchor (Rac=14.357230 mOhm, J=4.242996)
% S2 twist    : reproduce Q2 PEEC twist P=40 anchor (Rac=14.035939, J=4.086757)
% S3 braid    : valid geometry (spacing >= 2a), conservation, transposition effective
% S4 counter  : valid geometry, conservation
clear; clc;
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end

base=struct('ringmax',10,'p_mm',0.15,'d_mm',0.14,'sigma',5.8e7,'f',200e3,'I',20, ...
    'K',16,'Lambda_mm',40,'alpha_deg',25,'pitch_mm',40);
nF=0;
fprintf('========== q3_solver self-test ==========\n');

% S1 straight
c=base; c.scheme='straight'; tic; rS=q3_solver(c); tS=toc;
nF=chk(nF,abs(rS.Rac_peec_mohm-14.357230)<5e-4, sprintf('S1 straight Rac=%.6f (Q2 14.357230)',rS.Rac_peec_mohm));
nF=chk(nF,abs(rS.J_peec-4.242996)<5e-4, sprintf('S1 straight J=%.6f (Q2 4.242996)',rS.J_peec));
nF=chk(nF,abs(rS.Zint_ratio-1.001048147)<1e-6, sprintf('S1 isolated skin factor=%.9f',rS.Zint_ratio));
nF=chk(nF,abs(sum(rS.I_z)-20)<1e-6,'S1 current conservation');

% S2 twist
c=base; c.scheme='twist'; tic; rT=q3_solver(c); tT=toc;
nF=chk(nF,abs(rT.Rac_peec_mohm-14.035939)<5e-4, sprintf('S2 twist Rac=%.6f (Q2 14.035939)',rT.Rac_peec_mohm));
nF=chk(nF,abs(rT.J_peec-4.086757)<5e-4, sprintf('S2 twist J=%.6f (Q2 4.086757)',rT.J_peec));
nF=chk(nF,abs(rT.sbar_mean-1.015070)<1e-4, sprintf('S2 twist sbar=%.6f (Q2 1.015070)',rT.sbar_mean));

% S3 braid (spacing assertion inside solver must not fire)
c=base; c.scheme='braid'; c.alpha_deg=25; tic; rB=q3_solver(c); tB=toc;
nF=chk(nF,rB.min_center_dist_um>=139.99, sprintf('S3 braid min spacing=%.1f um (>=2a=140)',rB.min_center_dist_um));
nF=chk(nF,abs(sum(rB.I_z)-20)<1e-6,'S3 braid current conservation');
nF=chk(nF,all(isfinite([rB.J_peec rB.J_full rB.P_total_W_m rB.sbar_mean])),'S3 braid metrics finite');
nF=chk(nF,rB.center_current_mA>5, sprintf('S3 braid center current=%.2f mA (transposition works)',rB.center_current_mA));
nF=chk(nF,rB.J_peec<rS.J_peec-0.3, sprintf('S3 braid J=%.3f << straight %.3f',rB.J_peec,rS.J_peec));
nF=chk(nF,rB.sbar_mean>1.0&&rB.sbar_mean<1.6, sprintf('S3 braid sbar=%.4f in range',rB.sbar_mean));

% S4 counter_braid
c=base; c.scheme='counter_braid'; c.alpha_deg=25; tic; rC=q3_solver(c); tC=toc;
nF=chk(nF,rC.min_center_dist_um>=139.99, sprintf('S4 counter min spacing=%.1f um',rC.min_center_dist_um));
nF=chk(nF,abs(sum(rC.I_z)-20)<1e-6,'S4 counter current conservation');
nF=chk(nF,all(isfinite([rC.J_peec rC.J_full rC.P_total_W_m])),'S4 counter metrics finite');

fprintf('\n%-16s %10s %8s %8s %9s %8s %8s %7s %6s\n', ...
    'scheme','Rac(mO/m)','J_peec','J_full','P(W/m)','sbar','eta%','ctr mA','minUm');
prt('straight',rS); prt('twist P=40',rT); prt('braid a=25',rB); prt('counter a=25',rC);
fprintf('\ntime: S1 %.2fs S2 %.2fs S3 %.2fs S4 %.2fs\n',tS,tT,tB,tC);

if nF==0, fprintf('\n>>> ALL PASS: q3_solver ready for sweep.\n');
else, fprintf('\n>>> %d FAILED\n',nF); error('selftest failed'); end

function nF=chk(nF,cond,msg)
if cond, fprintf('  [PASS] %s\n',msg); else, fprintf('  [FAIL] %s\n',msg); nF=nF+1; end
end
function prt(name,r)
fprintf('%-16s %10.4f %8.4f %8.4f %9.4f %8.5f %8.1f %7.2f %6.1f\n', ...
    name,r.Rac_peec_mohm,r.J_peec,r.J_full,r.P_total_W_m,r.sbar_mean, ...
    r.eta_c_percent,r.center_current_mA,r.min_center_dist_um);
end
