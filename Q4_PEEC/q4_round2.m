function q4_round2()
%Q4_ROUND2  Corner smoothing (bend feasibility) vs area-fairness trade-off,
% plus larger-N (smaller-d) sweep, on the teammate full_exchange topology.
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end
A_Cu=6.0e-6; K=512; P=0.5;

% ---- (1) smoothing width sweep, side=22 (N=484) ----
side=22; hw=0.002284441841017299;
ws=[0 1 2 3 5 8 12 16];
R1={};
fprintf('===== Round2a: corner smoothing (side=%d, P=%.0fmm) =====\n',side,P*1e3);
fprintf('%3s %10s %10s %9s %9s %8s %8s %8s\n','w','bend_d','D_O','mind_um','J_peec','eta%','sbar','angle');
for w=ws
    r=q4_topology_teammate(side,hw,P,K,A_Cu,w);
    fprintf('%3d %10.2f %10.4f %9.1f %9.5f %8.2f %8.5f %8.1f\n', ...
        w,r.min_bend_d,r.D_O,r.mind_um,r.J,r.eta,r.sbar,r.max_angle_deg);
    R1(end+1,:)={w,r.min_bend_d,r.D_O,r.mind_um,r.J,r.eta,r.sbar,r.max_angle_deg}; %#ok<AGROW>
end
T1=cell2table(R1,'VariableNames',{'w','bend_d','D_O','mind_um','J_peec','eta_pct','sbar','angle_deg'});
writetable(T1,'q4_round2_smoothing.csv');

% pick smallest w with bend>=5d and mind>=2a
ok=(T1.bend_d>=5)&(T1.mind_um>=1e6*sqrt(4*A_Cu/(pi*side^2)));   % 2a = d
if any(ok), wbest=min(T1.w(ok)); else, wbest=NaN; end
fprintf('best feasible smoothing w = %g\n',wbest);

% ---- (2) larger N (smaller d) sweep, w = best ----
sides=[20 22 24 26]; hws=[0.00228 0.002284441841017299 0.001555 0.001580];
wuse=wbest; if isnan(wuse), wuse=0; end
R2={};
fprintf('\n===== Round2b: N sweep (w=%g, P=%.0fmm) =====\n',wuse,P*1e3);
fprintf('%4s %5s %8s %9s %9s %8s %8s %8s\n','side','N','d_um','J_peec','eta%','D_O','sbar','bend_d');
for si=1:numel(sides)
    r=q4_topology_teammate(sides(si),hws(si),P,K,A_Cu,wuse);
    fprintf('%4d %5d %8.2f %9.5f %8.2f %8.4f %8.5f %8.2f\n', ...
        r.side,r.N,r.d_um,r.J,r.eta,r.D_O,r.sbar,r.min_bend_d);
    R2(end+1,:)={r.side,r.N,r.d_um,r.J,r.eta,r.D_O,r.sbar,r.min_bend_d,r.mind_um}; %#ok<AGROW>
end
T2=cell2table(R2,'VariableNames',{'side','N','d_um','J_peec','eta_pct','D_O','sbar','bend_d','mind_um'});
writetable(T2,'q4_round2_Nsweep.csv');

% ---- figure ----
fh=figure('Color','w','Position',[60 60 1280 800]);
subplot(2,2,1); hold on; box on; grid on;
plot(T1.D_O,T1.bend_d,'-o','LineWidth',1.6,'MarkerFaceColor','b');
yline(5,'r--','bend\\geq5d');
xlabel('D_O (area-fairness defect)'); ylabel('min bend radius / d');
title('拐角平滑权衡：弯曲可行性 vs 面积公平');
for i=1:height(T1), text(T1.D_O(i),T1.bend_d(i),sprintf(' w=%d',T1.w(i)),'FontSize',8); end

subplot(2,2,2); hold on; box on; grid on;
plot(T1.w,T1.mind_um,'-o','MarkerFaceColor','g');
yline(1e6*sqrt(4*A_Cu/(pi*side^2)),'r--','2a');
xlabel('smoothing width w'); ylabel('min center spacing (\mum)');
title('平滑宽度对最小间距的影响');

subplot(2,2,3); hold on; box on; grid on;
yyaxis left; plot(T2.N,T2.J_peec,'-o','LineWidth',1.8,'MarkerFaceColor','b'); ylabel('J_{peec}');
yyaxis right; plot(T2.N,T2.d_um,'-s','LineWidth',1.5,'MarkerFaceColor','r'); ylabel('d (\mum)');
xlabel('N = side^2'); title('更大 N（更细丝径）的公平性口径 J');

subplot(2,2,4); hold on; box on; grid on;
bar(1:height(T2),T2.eta_pct,'FaceColor',[0.3 0.6 0.9]);
set(gca,'XTick',1:height(T2),'XTickLabel',arrayfun(@(n)sprintf('N=%d',n),T2.N,'UniformOutput',false));
ylabel('\eta_I (%)'); title('各 N 的均流偏差（应趋 0）');

sgtitle('问题四 Round2：拐角平滑 × 更大 N','FontWeight','bold');
exportgraphics(fh,'q4_round2.png','Resolution',150);
fprintf('\nSaved q4_round2_smoothing.csv, q4_round2_Nsweep.csv, q4_round2.png\n');
end
