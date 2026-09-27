%% q3_convergence.m — segment-count (K) convergence gate G1
% For braid and counter_braid at alpha=25deg, Lambda=40mm, check that the
% adjacent K levels change J/Rac by < 0.5% and lock K*.
clear; clc;
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end

base=struct('ringmax',10,'p_mm',0.15,'d_mm',0.14,'sigma',5.8e7,'f',200e3,'I',20, ...
    'Lambda_mm',40,'alpha_deg',25);
Ks=[16 32 64 128 256];
schemes={'braid','counter_braid'};

fprintf('========== G1: K convergence (alpha=25, Lambda=40) ==========\n');
rows={};
for s=1:numel(schemes)
    fprintf('\n--- %s ---\n', schemes{s});
    fprintf('%4s %12s %10s %8s %8s %9s\n','K','Rac(mOhm/m)','J','sbar','eta%','minUm');
    prevJ=NaN;
    for kk=1:numel(Ks)
        c=base; c.scheme=schemes{s}; c.K=Ks(kk);
        r=q3_solver(c);
        dJ=100*abs(r.J_peec-prevJ)/prevJ;
        fprintf('%4d %12.6f %10.5f %8.5f %8.1f %9.1f', ...
            Ks(kk), r.Rac_peec_mohm, r.J_peec, r.sbar_mean, r.eta_c_percent, r.min_center_dist_um);
        if isnan(prevJ), fprintf('\n'); else, fprintf('   dJ=%.4f%%\n', dJ); end
        rows(end+1,:)={schemes{s}, Ks(kk), r.Rac_peec_mohm, r.J_peec, r.sbar_mean, r.eta_c_percent, r.min_center_dist_um, dJ}; %#ok<AGROW>
        prevJ=r.J_peec;
    end
end

T=cell2table(rows,'VariableNames',{'scheme','K','Rac_mohm','J','sbar','eta_pct','min_um','adjacent_dJ_pct'});
writetable(T,'q3_convergence.csv');
fprintf('\nSaved q3_convergence.csv\n');
fprintf('Gate G1 rule: adjacent Rac change < 0.5%% for K>=16 (informational).\n');
