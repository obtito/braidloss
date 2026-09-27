function q4_round4()
%Q4_ROUND4  Corrected perfect-transposition family: spiral stations + exact
% cyclic shift. Scans W (winding), rin_frac (hollow), N; certifies D_O, eta_I,
% spacing, bend; exports certificates and best centerlines (T1).
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end
A_Cu=6.0e-6; P=0.25; K=256; N=484; R=2.3e-3;

fprintf('===== Round4: spiral full-exchange (N=%d, R=%.1fmm, P=%.0fmm) =====\n',N,R*1e3,P*1e3);
fprintf('%3s %6s %8s %8s %8s %8s %9s %8s %8s\n','W','rin','J_peec','eta%','D_O','sbar','mind_um','bend_d','2a_um');
R1={};
for rin=[0.30 0.45]
    for W=[8 12 16 20 28 40]
        r=q4_spiral(N,R,rin,W,P,K,A_Cu);
        fprintf('%3d %6.2f %8.5f %8.3f %8.4f %8.5f %9.1f %8.2f %8.1f\n', ...
            W,rin,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d,r.d_um);
        R1(end+1,:)={W,rin,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d,r.d_um}; %#ok<AGROW>
    end
end
T1=cell2table(R1,'VariableNames',{'W','rin_frac','J_peec','eta_pct','D_O','sbar','mind_um','bend_d','d_um'});
writetable(T1,'q4_round4_spiral.csv');

% feasible = spacing>=2a & bend>=5d & D_O<1e-3
ok=(T1.mind_um>=T1.d_um)&(T1.bend_d>=5)&(T1.D_O<1e-3);
fprintf('\nfeasible+perfect points: %d\n',sum(ok));
if any(ok)
    Tf=T1(ok,:); [~,i]=min(Tf.J_peec);
    Wb=Tf.W(i); rb=Tf.rin_frac(i);
    fprintf('best feasible: W=%d rin=%.2f  J=%.5f eta=%.3f%% D_O=%.2e sbar=%.5f mind=%.1fum bend=%.1fd\n', ...
        Wb,rb,Tf.J_peec(i),Tf.eta_pct(i),Tf.D_O(i),Tf.sbar(i),Tf.mind_um(i),Tf.bend_d(i));
else
    Wb=20; rb=0.45; fprintf('no feasible point; using W=%d rin=%.2f for T1 export\n',Wb,rb);
end

% N sweep at best
Ns=[331 400 484 576]; R2={};
fprintf('\n--- N sweep (rin=%.2f, W scaled with N) ---\n',rb);
for n=Ns
    Wn=max(6,round(Wb*484/n));   % keep tangential spacing ~constant
    r=q4_spiral(n,R,rb,Wn,P,K,A_Cu);
    fprintf('N=%3d W=%2d d=%.1f J=%.5f eta=%.3f%% D_O=%.2e sbar=%.5f mind=%.1f bend=%.1f\n', ...
        n,Wn,r.d_um,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d);
    R2(end+1,:)={n,Wn,r.d_um,r.J,r.eta,r.D_O,r.sbar,r.mind_um,r.min_bend_d}; %#ok<AGROW>
end
T2=cell2table(R2,'VariableNames',{'N','W','d_um','J_peec','eta_pct','D_O','sbar','mind_um','bend_d'});
writetable(T2,'q4_round4_Nsweep.csv');

% T1 centerlines + certificates for best feasible
best=q4_spiral(N,R,rb,Wb,P,max(K,512),A_Cu);
X=best.X; Y=best.Y; z=best.z; Kb=numel(z); n=best.N;
r=hypot(X,Y); th=atan2(Y,X);
dx=X(:,[2:Kb,1])-X(:,[Kb,1:Kb-1]); dy=Y(:,[2:Kb,1])-Y(:,[Kb,1:Kb-1]); dzb=z(2)-z(1);
S=sqrt(1+(dx/(2*dzb)).^2+(dy/(2*dzb)).^2);
fid=fopen('q4_best_centerlines.csv','w'); fprintf(fid,'strand,z,r,theta,s\n');
for kk=1:n, for j=1:Kb, fprintf(fid,'%d,%.6e,%.6e,%.6e,%.6f\n',kk-1,z(j),r(kk,j),th(kk,j),S(kk,j)); end, end
fclose(fid);
writetable(table(best.N,best.d_um,best.R_mm,best.rin_frac,best.W,best.J,best.Rac,best.eta,best.D_O,best.sbar,best.mind_um,best.min_bend_d, ...
    'VariableNames',{'N','d_um','R_mm','rin_frac','W','J_peec','Rac_ohm','eta_pct','D_O','sbar','mind_um','bend_d'}),'q4_certificates.csv');

% figure
fh=figure('Color','w','Position',[60 60 1280 780]);
subplot(2,2,1); hold on; box on; grid on;
for rf=unique(T1.rin_frac)'
    m=T1.rin_frac==rf; plot(T1.W(m),T1.mind_um(m),'-o','LineWidth',1.6);
end
yline(T1.d_um(1),'r--','2a'); xlabel('winding W'); ylabel('min spacing (\mum)');
title('螺旋绕数对最小间距的影响'); legend(arrayfun(@(x)sprintf('rin=%.2f',x),unique(T1.rin_frac),'UniformOutput',false));
subplot(2,2,2); hold on; box on; grid on;
for rf=unique(T1.rin_frac)'
    m=T1.rin_frac==rf; plot(T1.W(m),T1.J_peec(m),'-o','LineWidth',1.6);
end
xlabel('winding W'); ylabel('J_{peec}'); title('绕数对 J 的影响');
subplot(2,2,3); hold on; box on; grid on;
scatter(T1.D_O,T1.eta_pct,30,T1.W,'filled'); colorbar; xlabel('D_O'); ylabel('\eta_I (%)');
title('三证书：D_O 小且 \eta_I 小可同时满足');
subplot(2,2,4); hold on; box on; grid on;
yyaxis left; plot(T2.N,T2.J_peec,'-o','LineWidth',1.8,'MarkerFaceColor','b'); ylabel('J_{peec}');
yyaxis right; plot(T2.N,T2.d_um,'-s','MarkerFaceColor','r'); ylabel('d (\mum)');
xlabel('N'); title('更大 N（更细丝径）');
sgtitle('问题四 Round4：螺旋站位 + 精确循环移位（修正版）','FontWeight','bold');
exportgraphics(fh,'q4_round4.png','Resolution',150);
fprintf('\nSaved q4_round4_spiral.csv, q4_round4_Nsweep.csv, q4_best_centerlines.csv, q4_certificates.csv, q4_round4.png\n');
end
