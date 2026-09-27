function q3_tune_at_comsol()
%Q3_TUNE_AT_COMSOL  Re-run the PEEC tuning at the COMSOL design point and plot.
% COMSOL (litz_q3) design point: A_Cu = 6.0 mm^2, N ~ 400, d linked by area,
% f = 200 kHz, I = 20 A rms. Hex layout with ringmax=11 gives N=397 (~400).
% Panels: J(alpha), J(Lambda), scheme bars at the design point, and the direct
% PEEC-vs-COMSOL K comparison on the COMSOL geometry.
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end

A_Cu=6.0e-6; Ndes=397;                       % hex ringmax=11
d=sqrt(4*A_Cu/(pi*Ndes));                    % linked diameter
base=struct('ringmax',11,'p_mm',0.16,'d_mm',d*1e3,'sigma',5.8e7,'f',200e3,'I',20, ...
    'K',256,'Lambda_mm',120,'alpha_deg',15,'pitch_mm',250);
fprintf('COMSOL design point: A_Cu=6.0 mm^2, N=%d, d=%.4f mm, d/delta=%.4f\n', ...
    Ndes, d*1e3, d/sqrt(2/(2*pi*base.f*4*pi*1e-7*base.sigma)));

% anchors at design point
rS=q3_solver(setfield(base,'scheme','straight'));
rT=q3_solver(setfield(base,'scheme','twist'));
fprintf('anchors: straight J=%.4f Rac=%.4f | twist(P=250) J=%.4f Rac=%.4f\n', ...
    rS.J_peec,rS.Rac_peec_mohm,rT.J_peec,rT.Rac_peec_mohm);

alphas=[5 10 15 20 25 30];
schemes={'braid','counter_braid'};
Jb=zeros(size(alphas)); Jc=zeros(size(alphas)); rows={};
for s=1:2
    for a=1:numel(alphas)
        c=base; c.scheme=schemes{s}; c.alpha_deg=alphas(a);
        r=q3_solver(c);
        if s==1, Jb(a)=r.J_peec; else, Jc(a)=r.J_peec; end
        rows(end+1,:)={schemes{s},alphas(a),120,r.J_peec,r.Rac_peec_mohm,r.P_total_W_m,r.sbar_mean,r.eta_c_percent}; %#ok<AGROW>
    end
end
% Lambda sweep at best alpha (braid)
[~,ia]=min(Jb); aBest=alphas(ia);
Lams=[40 80 120 250]; JL=zeros(size(Lams));
for L=1:numel(Lams)
    c=base; c.scheme='braid'; c.alpha_deg=aBest; c.Lambda_mm=Lams(L);
    r=q3_solver(c); JL(L)=r.J_peec;
    rows(end+1,:)={'braid',aBest,Lams(L),r.J_peec,r.Rac_peec_mohm,r.P_total_W_m,r.sbar_mean,r.eta_c_percent}; %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'scheme','alpha_deg','Lambda_mm','J','Rac_mohm','P_W_m','sbar','eta_pct'});
writetable(T,'q3_tune_at_comsol.csv');
[~,il]=min(JL);
cb=base; cb.scheme='braid'; cb.alpha_deg=aBest; cb.Lambda_mm=Lams(il);
rBest=q3_solver(cb);

% ---- figure ----
fh=figure('Color','w','Position',[60 60 1340 860]);

% (1) J vs alpha
subplot(2,2,1); hold on; box on; grid on;
plot(alphas,Jb,'-o','LineWidth',1.8,'MarkerFaceColor','b');
plot(alphas,Jc,'-s','LineWidth',1.8,'MarkerFaceColor','r');
plot(alphas,rS.J_peec*ones(size(alphas)),'k--','LineWidth',1.2);
plot(alphas,rT.J_peec*ones(size(alphas)),'--','Color',[.4 .4 .4],'LineWidth',1.2);
xlabel('\alpha (deg)'); ylabel('J = R_{ac}/R_{dc}');
title(sprintf('编织角敏感性（A_{Cu}=6.0mm^2, N=%d, K=256, \\Lambda=120mm）',Ndes));
legend({'braid','counter\_braid','straight','twist P=250'},'Location','northwest');

% (2) J vs Lambda
subplot(2,2,2); hold on; box on; grid on;
plot(Lams,JL,'-^','LineWidth',1.8,'MarkerFaceColor','g');
xlabel('\Lambda (mm)'); ylabel('J = R_{ac}/R_{dc}');
title(sprintf('换位周期敏感性（braid, \\alpha=%d°）',aBest));
ylim([min(JL)*0.995 max(JL)*1.02]);

% (3) scheme bars at design point
subplot(2,2,3); hold on; box on; grid on;
vals=[rS.J_peec rT.J_peec min(Jb) min(Jc)];
b=bar(vals,'FaceColor','flat');
b.CData=[0.5 0.5 0.5;0.7 0.7 0.7;0.2 0.5 0.9;0.9 0.4 0.2];
set(gca,'XTickLabel',{'straight','twist','braid','counter'});
ylabel('J = R_{ac}/R_{dc}');
title('COMSOL 设计点各方案 J（PEEC，K=256）');
for k=1:4, text(k,vals(k)+0.05,sprintf('%.3f',vals(k)),'HorizontalAlignment','center'); end

% (4) PEEC vs COMSOL on COMSOL geometry
subplot(2,2,4); hold on; box on; grid on;
if isfile('litz_peec_vs_comsol.csv')
    L=readtable('litz_peec_vs_comsol.csv','VariableNamingRule','preserve');
    lbl=L{:,1}; kpeec=L{:,3}; kcom=L{:,4};
    x=1:numel(kpeec);
    bar(x-0.18,kpeec,0.36,'FaceColor',[0.25 0.6 0.85]);
    bar(x+0.18,kcom,0.36,'FaceColor',[0.85 0.45 0.25]);
    set(gca,'XTick',x,'XTickLabel',strrep(lbl,'_','\_'));
    ylabel('K = R_{ac}/R_{dc}');
    legend({'PEEC','COMSOL'},'Location','northeast');
    title('同一 COMSOL 几何上 PEEC vs COMSOL');
end

sgtitle(sprintf('问题三 MatLab/PEEC 调优 @ COMSOL 设计点：A_{Cu}=6.0mm^2, N=%d, d=%.1f\\mum, f=200kHz, I=20A', ...
    Ndes, d*1e6),'FontWeight','bold');
exportgraphics(fh,'q3_tune_at_comsol.png','Resolution',150);
fprintf('best braid at design: alpha=%d Lambda=%d J=%.5f Rac=%.4f P=%.4f sbar=%.5f\n', ...
    aBest,Lams(il),rBest.J_peec,rBest.Rac_peec_mohm,rBest.P_total_W_m,rBest.sbar_mean);
fprintf('saved q3_tune_at_comsol.png and q3_tune_at_comsol.csv\n');
end
