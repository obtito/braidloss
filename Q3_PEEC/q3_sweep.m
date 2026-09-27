function q3_sweep()
%Q3_SWEEP  Q3 tuning loop: scheme x alpha scan, then Lambda sensitivity.
% Locks K=256 (G1: braid 0.043%, counter 0.135% adjacent change < 0.5%).
% Writes q3_scan_alpha.csv, q3_sensitivity.csv, q3_optimal.json, q3_tuning_log.txt.
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end
logf=fopen('q3_tuning_log.txt','a');
fin=onCleanup(@() fclose(logf));
logline(logf,'===== Q3 sweep start =====');

base=struct('ringmax',10,'p_mm',0.15,'d_mm',0.14,'sigma',5.8e7,'f',200e3,'I',20, ...
    'K',256,'Lambda_mm',40,'alpha_deg',25,'pitch_mm',40);

% ---- anchors ----
rS=q3_solver(setfield(base,'scheme','straight'));
rT=q3_solver(setfield(base,'scheme','twist'));
logline(logf,'anchors: straight J=%.4f Rac=%.4f | twist J=%.4f Rac=%.4f', ...
    rS.J_peec,rS.Rac_peec_mohm,rT.J_peec,rT.Rac_peec_mohm);

% ---- Round 1: alpha scan ----
alphas=[2 4 6 8 10 15 20 25 30 35 40 45];
schemes={'braid','counter_braid'};
rows={};
rows=addrow(rows,'straight',NaN,rS);
rows=addrow(rows,'twist',NaN,rT);
for s=1:numel(schemes)
    for a=1:numel(alphas)
        c=base; c.scheme=schemes{s}; c.alpha_deg=alphas(a);
        r=q3_solver(c);
        logline(logf,'%-14s alpha=%5.1f  J=%.5f  Rac=%.4f  P=%.4f  sbar=%.5f  eta=%.1f%%  ctr=%.2fmA  min=%.1fum', ...
            schemes{s},alphas(a),r.J_peec,r.Rac_peec_mohm,r.P_total_W_m,r.sbar_mean,r.eta_c_percent,r.center_current_mA,r.min_center_dist_um);
        rows=addrow(rows,schemes{s},alphas(a),r);
    end
end
T=cell2table(rows,'VariableNames',{'scheme','alpha_deg','J','Rac_mohm','P_W_m','sbar','eta_pct','center_mA','outer_loss_pct','min_um'});
writetable(T,'q3_scan_alpha.csv');
logline(logf,'Round1 done: %d rows -> q3_scan_alpha.csv',height(T));

% ---- pick best alpha per active scheme ----
best=struct('scheme',{},'alpha',{},'J',{});
for s=1:numel(schemes)
    m=strcmp(T.scheme,schemes{s});
    [jmin,idx]=min(T.J(m));
    aa=T.alpha_deg(m); aa=aa(idx);
    best(s).scheme=schemes{s}; best(s).alpha=aa; best(s).J=jmin;
    logline(logf,'best %s: alpha*=%.1f J*=%.5f',schemes{s},aa,jmin);
end

% ---- Round 2: Lambda sensitivity at each best alpha ----
Lams=[20 30 40 60 80 120 160 240];
srows={};
for s=1:numel(best)
    for L=1:numel(Lams)
        c=base; c.scheme=best(s).scheme; c.alpha_deg=best(s).alpha; c.Lambda_mm=Lams(L);
        r=q3_solver(c);
        logline(logf,'%-14s alpha=%.1f Lambda=%2d  J=%.5f  Rac=%.4f  sbar=%.5f  eta=%.1f%%', ...
            best(s).scheme,best(s).alpha,Lams(L),r.J_peec,r.Rac_peec_mohm,r.sbar_mean,r.eta_c_percent);
        srows(end+1,:)={best(s).scheme,best(s).alpha,Lams(L),r.J_peec,r.Rac_peec_mohm,r.P_total_W_m,r.sbar_mean,r.eta_c_percent,r.center_current_mA,r.min_center_dist_um}; %#ok<AGROW>
    end
end
Ts=cell2table(srows,'VariableNames',{'scheme','alpha_deg','Lambda_mm','J','Rac_mohm','P_W_m','sbar','eta_pct','center_mA','min_um'});
writetable(Ts,'q3_sensitivity.csv');
logline(logf,'Round2 done: %d rows -> q3_sensitivity.csv',height(Ts));

% ---- optimal ----
[jopt,iopt]=min(Ts.J);
opt=struct('scheme',Ts.scheme{iopt},'alpha_deg',Ts.alpha_deg(iopt),'Lambda_mm',Ts.Lambda_mm(iopt), ...
    'J',jopt,'Rac_mohm',Ts.Rac_mohm(iopt),'P_W_m',Ts.P_W_m(iopt),'sbar',Ts.sbar(iopt), ...
    'eta_pct',Ts.eta_pct(iopt),'race_note','monotone to floor; see q3_sweep_outcomes.md for caveats');
[~,big]=min([rS.J_peec rT.J_peec]);
opt.baseline_straight_J=rS.J_peec; opt.baseline_twist_J=rT.J_peec;
opt.J_reduction_vs_straight_pct=100*(rS.J_peec-jopt)/rS.J_peec;
fid=fopen('q3_optimal.json','w'); fprintf(fid,'%s',jsonencode(opt)); fclose(fid);
logline(logf,'OPTIMAL: %s alpha=%.1f Lambda=%.1f J=%.5f (vs straight %.4f, -%.1f%%)', ...
    Ts.scheme{iopt},Ts.alpha_deg(iopt),Ts.Lambda_mm(iopt),jopt,rS.J_peec,opt.J_reduction_vs_straight_pct);
logline(logf,'===== Q3 sweep end =====');

fprintf('\n===== scan summary =====\n');
disp(T);
fprintf('===== Lambda sensitivity (best alpha) =====\n');
disp(Ts);
fprintf('OPTIMAL: %s alpha=%.1f Lambda=%.1f J=%.5f  (straight %.4f, -%.1f%%)\n', ...
    Ts.scheme{iopt},Ts.alpha_deg(iopt),Ts.Lambda_mm(iopt),jopt,rS.J_peec,opt.J_reduction_vs_straight_pct);
end

function rows=addrow(rows,scheme,alpha,r)
rows(end+1,:)={scheme,alpha,r.J_peec,r.Rac_peec_mohm,r.P_total_W_m,r.sbar_mean, ...
    r.eta_c_percent,r.center_current_mA,r.outer_loss_share_percent,r.min_center_dist_um};
end

function logline(fid,fmt,varargin)
msg=sprintf(fmt,varargin{:});
fprintf('%s\n',msg);
fprintf(fid,'%s  %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS'),msg);
end
