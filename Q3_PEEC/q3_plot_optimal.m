function q3_plot_optimal()
%Q3_PLOT_OPTIMAL  Visualise the current optimal scheme with its parameters.
% Reads q3_optimal.json (fallback braid alpha=2 Lambda=240), re-runs q3_solver,
% and draws a 2x2 figure: 3D trajectories (one transposition period), parameter
% + metric panel, radial layer path r(z), cross-section snapshots.
% Saves q3_optimal_scheme.png
thisDir=fileparts(mfilename('fullpath')); if ~isempty(thisDir), cd(thisDir); end

if isfile('q3_optimal.json')
    o=jsondecode(fileread('q3_optimal.json'));
    scheme=o.scheme; alpha=o.alpha_deg; Lam=o.Lambda_mm;
else
    scheme='braid'; alpha=2; Lam=240;
end
cfg=struct('ringmax',10,'p_mm',0.15,'d_mm',0.14,'sigma',5.8e7,'f',200e3,'I',20, ...
    'K',256,'Lambda_mm',Lam,'alpha_deg',alpha,'pitch_mm',40,'scheme',scheme);
r=q3_solver(cfg);
rS=q3_solver(setfield(cfg,'scheme','straight'));
rT=q3_solver(setfield(cfg,'scheme','twist'));

X=r.X; Y=r.Y; N=r.N; K=r.K; ring=r.ring;
z1=((0:K-1)+0.5)*(Lam/K);          % mm, one period (Lam is in mm)
rpath=hypot(X,Y)*1e3;              % mm

fh=figure('Color','w','Position',[60 60 1320 860],'Name','Q3 optimal scheme');

% (1) 3D trajectories over one period
ax1=subplot(2,2,1); hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on');
for i=1:N
    plot3(ax1,X(i,:)*1e3, Y(i,:)*1e3, z1, 'Color',[0.82 0.85 0.9], 'LineWidth',0.2);
end
sel=[find(ring==0,1), find(ring==1,1), find(ring==3,1), find(ring==5,1), ...
     find(ring==7,1), find(ring==9,1), find(ring==10,1)];
sel=sel(~isnan(sel));
cols=lines(numel(sel));
for j=1:numel(sel)
    plot3(ax1,X(sel(j),:)*1e3, Y(sel(j),:)*1e3, z1, 'Color',cols(j,:), 'LineWidth',1.6);
end
xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
title(sprintf('最优方案三维轨迹：%s，\\alpha=%g°，\\Lambda=%d mm（1 个换位周期）', scheme, alpha, Lam));
view(40,20); axis tight; pbaspect([1 1 1]); grid on;
legend(arrayfun(@(k)sprintf('起始环 %d',k), ring(sel), 'UniformOutput',false), ...
    'Location','northeast','FontSize',7);

% (2) parameter + metric panel
ax2=subplot(2,2,2); axis(ax2,'off');
txt={ sprintf('scheme : %s', scheme), ...
      sprintf('alpha  : %g deg', alpha), ...
      sprintf('Lambda : %d mm', Lam), ...
      sprintf('K      : %d segments/period', K), ...
      sprintf('N      : %d strands', N), ...
      sprintf('d      : %.2f mm,  p = 0.15 mm', cfg.d_mm), ...
      sprintf('A_Cu   : 5.095 mm^2'), ...
      sprintf('f      : 200 kHz,  I = 20 A rms'), ...
      '----------------------------------------', ...
      sprintf('Rac      = %.4f mOhm/m', r.Rac_peec_mohm), ...
      sprintf('Rac/Rdc  = %.5f', r.J_peec), ...
      sprintf('P        = %.4f W/m', r.P_total_W_m), ...
      sprintf('eta_I    = %.1f %%', r.eta_c_percent), ...
      sprintf('sbar     = %.5f', r.sbar_mean), ...
      sprintf('center I = %.2f mA (ideal 60.42)', r.center_current_mA), ...
      '----------------------------------------', ...
      sprintf('vs straight  %.5f  (%.1f%%)', rS.J_peec, 100*(r.J_peec-rS.J_peec)/rS.J_peec), ...
      sprintf('vs twist     %.5f  (%.1f%%)', rT.J_peec, 100*(r.J_peec-rT.J_peec)/rT.J_peec) };
text(ax2,0.02,0.98,txt,'VerticalAlignment','top','Interpreter','none', ...
    'FontName','Consolas','FontSize',9.5,'BackgroundColor',[1 1 0.93],'EdgeColor',[0.4 0.4 0.4]);
title(ax2,'方案参数与指标');

% (3) radial layer path
ax3=subplot(2,2,3); hold(ax3,'on'); box(ax3,'on'); grid(ax3,'on');
for j=1:numel(sel)
    plot(ax3,[z1 z1+Lam],[rpath(sel(j),:) rpath(sel(j),:)],'Color',cols(j,:),'LineWidth',1.5);
end
xlabel('z (mm)'); ylabel('r (mm)'); ylim([0 1.6]);
title('径向层路径 r(z)：股线周期遍历内外层');

% (4) cross-section snapshots
ax4=subplot(2,2,4); hold(ax4,'on'); box(ax4,'on'); grid(ax4,'on');
snaps=[1 round(K/4) round(K/2) round(3*K/4)];
scol=[0.20 0.35 0.75; 0.85 0.33 0.10; 0.20 0.65 0.30; 0.55 0.25 0.70];
leg={};
for s=1:numel(snaps)
    plot(ax4,X(:,snaps(s))*1e3, Y(:,snaps(s))*1e3, '.', 'Color',scol(s,:), 'MarkerSize',8);
    leg{s}=sprintf('z = %.0f mm', z1(snaps(s)));
end
axis(ax4,'equal'); xlabel('x (mm)'); ylabel('y (mm)');
title('截面快照：每环股数守恒，间距 \geq 2a');
legend(leg,'Location','northeastoutside','FontSize',8);

exportgraphics(fh,'q3_optimal_scheme.png','Resolution',150);
fprintf('optimal scheme: %s alpha=%g Lambda=%d K=%d\n', scheme, alpha, Lam, K);
fprintf('Rac=%.4f mOhm/m  J=%.5f  P=%.4f W/m  eta=%.1f%%  sbar=%.5f\n', ...
    r.Rac_peec_mohm, r.J_peec, r.P_total_W_m, r.eta_c_percent, r.sbar_mean);
fprintf('saved q3_optimal_scheme.png\n');
end
