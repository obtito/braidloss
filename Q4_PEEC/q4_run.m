function q4_run()
%Q4_RUN  Approximation iteration on the teammate full_exchange topology.
% Sweeps grid side (N=400,484) and period P, certifies fairness/feasibility,
% evaluates J by PEEC, ranks, and writes CSV + figure.
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end
A_Cu=6.0e-6; K=512;
sides=[20 22]; hw=[0.00228 0.002284441841017299];
Ps=[0.1 0.1666666667 0.25 0.3333333333 0.5];

% teammate COMSOL full_exchange references (topology_control / final_results)
comsol=[20 0.25 1.5961987792728285; 22 0.5 1.4971216482038827];

fprintf('===== Q4 iteration on teammate full_exchange =====\n');
fprintf('%4s %5s %8s %8s %8s %8s %9s %8s %8s %8s\n','side','N','P(mm)','d(um)','J_peec','eta%','D_O','sbar','mind_um','bend_d');
rows={};
for si=1:numel(sides)
    for pi=1:numel(Ps)
        r=q4_topology_teammate(sides(si),hw(si),Ps(pi),K,A_Cu);
        fprintf('%4d %5d %8.1f %8.2f %8.5f %8.2f %9.2e %8.5f %8.1f %8.1f\n', ...
            r.side,r.N,r.P_mm,r.d_um,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d);
        rows(end+1,:)={r.side,r.N,r.P_mm,r.d_um,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d,r.max_angle_deg}; %#ok<AGROW>
    end
end
T=cell2table(rows,'VariableNames',{'side','N','P_mm','d_um','J_peec','eta_pct','D_O','sbar','mind_um','min_bend_d','max_angle_deg'});
writetable(T,'q4_teammate_sweep.csv');
fprintf('\nSaved q4_teammate_sweep.csv\n');

% ---- figure ----
fh=figure('Color','w','Position',[60 60 1200 760]);
subplot(1,2,1); hold on; box on; grid on;
for si=1:numel(sides)
    m=T.side==sides(si);
    plot(T.P_mm(m),T.J_peec(m),'-o','LineWidth',1.8,'MarkerFaceColor','b');
end
xlabel('period P (mm)'); ylabel('J_{peec} = R_{ac}/R_{dc}');
title('队友 full\_exchange 拓扑：PEEC 公平性口径 J(P)');
legend(arrayfun(@(s)sprintf('N=%d',s^2),sides,'UniformOutput',false),'Location','northeast');

subplot(1,2,2); hold on; box on; grid on;
x=1:size(comsol,1);
bar(x-0.18,[1.5962;1.4971],0.36,'FaceColor',[0.85 0.45 0.25]);
bar(x+0.18,[T.J_peec(T.side==20 & T.P_mm==250/1000); T.J_peec(T.side==22 & abs(T.P_mm-500)<1e-6)],0.36,'FaceColor',[0.25 0.6 0.85]);
set(gca,'XTick',x,'XTickLabel',{'N=400 P=250','N=484 P=500'});
ylabel('R_{ac}/R_{dc}');
legend({'COMSOL K (含股内邻近)','PEEC J (股间+趋肤)'},'Location','northwest');
title('同几何：COMSOL vs PEEC（差额=股内邻近 \Phi 目标）');

sgtitle('问题四 逼近迭代：队友 full\_exchange 拓扑复算与证书','FontWeight','bold');
exportgraphics(fh,'q4_teammate_iteration.png','Resolution',150);
disp(T);
end
