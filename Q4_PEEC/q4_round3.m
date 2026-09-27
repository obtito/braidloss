function q4_round3()
%Q4_ROUND3  Area reparameterisation (D_O->0) + hollow packing + smaller d,
% then export the best candidate centerlines (T1) and certificates.
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end
A_Cu=6.0e-6; K=512; P=0.5; side=22; hw=0.002284441841017299;

% (a) baseline (teammate loop) vs area-reparameterised
rb=q4_topology_teammate(side,hw,P,K,A_Cu,0);
ra=q4_arearep(side,hw,0,P,K,A_Cu);
fprintf('=== Round3a: teammate-loop vs area-reparam (N=%d) ===\n',side^2);
fprintf('%-12s %8s %8s %8s %8s %9s %8s\n','variant','J_peec','eta%','D_O','sbar','mind_um','bend_d');
fprintf('%-12s %8.5f %8.3f %8.4f %8.5f %9.1f %8.2f\n','loop',rb.J,rb.eta,rb.D_O,rb.sbar,rb.mind_um,rb.min_bend_d);
fprintf('%-12s %8.5f %8.3f %8.4f %8.5f %9.1f %8.2f\n','area-rep',ra.J,ra.eta,ra.D_O,ra.sbar,ra.mind_um,ra.min_bend_d);

% (b) hollow sweep (area-rep)
rins=[0 0.2 0.4 0.6 0.75];
R2={};
fprintf('\n=== Round3b: hollow sweep (area-rep, N=%d) ===\n',side^2);
fprintf('%6s %8s %8s %8s %8s %9s %8s\n','rin','J_peec','eta%','D_O','sbar','mind_um','bend_d');
for rf=rins
    r=q4_arearep(side,hw,rf,P,K,A_Cu);
    fprintf('%6.2f %8.5f %8.3f %8.4f %8.5f %9.1f %8.2f\n',rf,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d);
    R2(end+1,:)={rf,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d}; %#ok<AGROW>
end
T2=cell2table(R2,'VariableNames',{'rin_frac','J_peec','eta_pct','D_O','sbar','mind_um','bend_d'});
writetable(T2,'q4_round3_hollow.csv');

% (c) N sweep at rin=0.4
sides=[20 22 24 26]; hws=[0.00228 hw 0.001555 0.001580];
R3={};
fprintf('\n=== Round3c: N sweep (area-rep, rin=0.4) ===\n');
fprintf('%4s %5s %8s %8s %8s %8s %9s %8s\n','side','N','d_um','J_peec','eta%','D_O','mind_um','bend_d');
for si=1:numel(sides)
    r=q4_arearep(sides(si),hws(si),0.4,P,K,A_Cu);
    fprintf('%4d %5d %8.2f %8.5f %8.3f %8.4f %9.1f %8.2f\n', ...
        r.side,r.N,r.d_um,r.J,r.eta,r.D_O,r.mind_um,r.min_bend_d);
    R3(end+1,:)={r.side,r.N,r.d_um,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d}; %#ok<AGROW>
end
T3=cell2table(R3,'VariableNames',{'side','N','d_um','J_peec','eta_pct','D_O','sbar','mind_um','bend_d'});
writetable(T3,'q4_round3_Nsweep.csv');

% (d) T1 centerlines + certificates for best (area-rep rin=0.4, N=484)
best=q4_arearep(side,hw,0.4,P,K,A_Cu);
n=best.N; X=best.X; Y=best.Y; z=best.z;
r=hypot(X,Y); th=atan2(Y,X);
S=zeros(size(X));
dx=X(:,[2:end,1])-X(:,[end,1:end-1]); dy=Y(:,[2:end,1])-Y(:,[end,1:end-1]);
dz=P/K; Sv=sqrt(1+(dx/(2*dz)).^2+(dy/(2*dz)).^2); S=Sv;
[Zg,Sg]=ndgrid(z,1:n);
fid=fopen('q4_best_centerlines.csv','w');
fprintf(fid,'strand,z,r,theta,s\n');
for k=1:n
    for j=1:K
        fprintf(fid,'%d,%.6e,%.6e,%.6e,%.6f\n',k-1,z(j),r(k,j),th(k,j),S(k,j));
    end
end
fclose(fid);
writetable(table(best.N,best.d_um,best.P_mm,best.rin_frac,best.J,best.Rac,best.eta,best.D_O,best.sbar,best.mind_um,best.min_bend_d, ...
    'VariableNames',{'N','d_um','P_mm','rin_frac','J_peec','Rac_ohm','eta_pct','D_O','sbar','mind_um','bend_d'}),'q4_certificates.csv');

% figure
fh=figure('Color','w','Position',[60 60 1280 780]);
subplot(2,2,1); hold on; box on; grid on;
b=bar([rb.D_O ra.D_O; rb.D_O_inter ra.D_O_inter]','grouped');
set(gca,'XTickLabel',{'teammate loop','area-rep'}); ylabel('D_O'); legend({'area defect','inter-strand spread'},'Location','best');
title('面积重参数化把 D_O(面积) 压到 ~0');
subplot(2,2,2); hold on; box on; grid on;
plot(T2.rin_frac,T2.J_peec,'-o','LineWidth',1.8,'MarkerFaceColor','b'); xlabel('rin/R (hollow)'); ylabel('J_{peec}');
title('空心度对 J 的影响（面积重参数化）');
subplot(2,2,3); hold on; box on; grid on;
yyaxis left; plot(T3.N,T3.J_peec,'-o','LineWidth',1.8,'MarkerFaceColor','b'); ylabel('J_{peec}');
yyaxis right; plot(T3.N,T3.d_um,'-s','MarkerFaceColor','r'); ylabel('d (\mum)');
xlabel('N'); title('更大 N（更细丝径）');
subplot(2,2,4); hold on; box on; grid on;
bar(1:height(T2),T2.mind_um,'FaceColor',[0.3 0.7 0.4]);
yline(1e6*sqrt(4*A_Cu/(pi*side^2)),'r--','2a');
set(gca,'XTick',1:height(T2),'XTickLabel',arrayfun(@(x)sprintf('%.2f',x),T2.rin_frac,'UniformOutput',false));
xlabel('rin/R'); ylabel('min spacing (\mum)'); title('空心度对最小间距的影响');
sgtitle('问题四 Round3：面积重参数化 × 空心化 × 减 d','FontWeight','bold');
exportgraphics(fh,'q4_round3.png','Resolution',150);
fprintf('\nSaved q4_round3_hollow.csv, q4_round3_Nsweep.csv, q4_best_centerlines.csv, q4_certificates.csv, q4_round3.png\n');
end
