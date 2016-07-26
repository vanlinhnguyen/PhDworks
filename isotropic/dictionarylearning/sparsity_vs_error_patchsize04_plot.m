clear all; clc; close all;


%%
% fsize=28;
% fname='CMU Serif';
% 
% h=figure;
% set(h, 'Position', [100 100 800 600]);
% set(h,'color','w')
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% %
% load sparsity_vs_error_testing.mat
% ax(1)=subplot(1,2,1,'position',[0.15 0.15 0.75 0.75]); % top left
% % plot here
% hold on
% h1=plot(thres_pca, NRMSE_PCA_all,'r-s','LineWidth',3,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);  
% h2=plot(NON_ZERO_COEFFS_0DL_all./256, NRMSE_0DL_all,'k-^','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);  
% h3=plot(NON_ZERO_COEFFS_KSVD_all./256, NRMSE_KSVD_all,'b-o','LineWidth',3,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
% h4=plot(NON_ZERO_COEFFS_WL_all./128^2, NRMSE_WL_all,'m-s','LineWidth',3,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6);  
% hold off
% 
% xlabel('$Sparsity$','Interpreter','latex')
% ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
% ylim([0 0.3]);
% xlim([0 1]);
% 
% set(gca, 'YTickMode','manual');
% set(gca, 'YTick', 0:0.05:0.3);
% set (gca, 'YTickLabel', {'0.0','','0.1','','0.2','','0.3'},'FontSize',fsize)
% 
% set(gca, 'XTickMode','manual');
% set(gca, 'XTick', 0:0.2:1);
% set (gca, 'XTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)
% 
% leg =legend ([h1,h2,h3,h4],{'PCA','ODL','KSVD','Wavelet(Meyer)'},'interpreter', 'latex','location','northeast','FontSize',fsize); 
% 
% set(leg,'FontSize',fsize-2);
% legend boxoff
% box on
% grid on
% 
% 
% export_fig('./figures/Sparsity_vs_NRMSE_PCA_ODL_KSVD_WL_Meyer_patchsize04','-eps');
% close();

%%
% fsize=28;
% fname='CMU Serif';
% 
% h=figure;
% set(h, 'Position', [100 100 800 600]);
% set(h,'color','w')
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% %
% load sparsity_vs_error_testing.mat
% ax(1)=subplot(1,2,1,'position',[0.15 0.15 0.75 0.75]); % top left
% h1=loglog(thres_pca, NRMSE_PCA_all,'r-s','LineWidth',3,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);  
% hold on
% h2=loglog(NON_ZERO_COEFFS_0DL_all./256, NRMSE_0DL_all,'k-^','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);  
% h3=loglog(NON_ZERO_COEFFS_KSVD_all./256, NRMSE_KSVD_all,'b-o','LineWidth',3,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
% h4=loglog(NON_ZERO_COEFFS_WL_all./128^2, NRMSE_WL_all,'m-s','LineWidth',3,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6);  
% hold off
% 
% xlabel('$Sparsity$','Interpreter','latex')
% ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
% ylim([3*1e-5 2*1e-1]);
% xlim([1e-1 1e0]);
% % 
% % set(gca, 'YTickMode','manual');
% % set(gca, 'YTick', 0:0.05:0.3);
% % set (gca, 'YTickLabel', {'0.0','','0.1','','0.2','','0.3'},'FontSize',fsize)
% % 
% % set(gca, 'XTickMode','manual');
% % set(gca, 'XTick', 0:0.2:1);
% % set (gca, 'XTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)
% 
% leg =legend ([h1,h2,h3,h4],{'PCA','ODL','KSVD','Wavelet(Meyer)'},'interpreter', 'latex','location','southwest','FontSize',fsize); 
% 
% set(leg,'FontSize',fsize-2);
% legend boxoff
% box on
% grid on
% 
% 
% export_fig('./figures/Sparsity_vs_NRMSE_PCA_ODL_KSVD_WL_Meyer_patchsize04_zoom','-eps');
% close();

%%
fsize=26;
fname='CMU Serif';

load sparsity_vs_error_testing.mat

h=figure;
set(h, 'Position', [100 100 1400 600]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% 
ax(1)=subplot(1,2,1,'position',[0.075 0.15 0.4 0.7]); % top left

hold on
h1=plot(1-thres_pca, NRMSE_PCA_all,'r-s','LineWidth',3,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);  
h2=plot(1-NON_ZERO_COEFFS_0DL_all./256, NRMSE_0DL_all,'k-^','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);  
h3=plot(1-NON_ZERO_COEFFS_KSVD_all./256, NRMSE_KSVD_all,'b-o','LineWidth',3,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
h4=plot(1-NON_ZERO_COEFFS_WL_all./128^2, NRMSE_WL_all,'m-s','LineWidth',3,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6);  
hold off
set ( gca, 'xdir', 'reverse' )
xlabel('$Sparsity$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
ylim([0 0.3]);
xlim([0 1]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.3);
set (gca, 'YTickLabel', {'0.0','','0.1','','0.2','','0.3'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.2:1);
set (gca, 'XTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)


leg =legend ([h1,h2,h3,h4],{'PCA','ODL','KSVD','Wavelet(Meyer)'},'interpreter', 'latex','location','northeast','FontSize',fsize); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on
grid on

ax(2)=axes('position',[0.55 0.15 0.4 0.7]); % top right
set ( gca, 'xdir', 'reverse' )

% plot here
h1=semilogy(1-thres_pca, NRMSE_PCA_all,'r-s','LineWidth',3,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);  
hold on
h2=semilogy(1-NON_ZERO_COEFFS_0DL_all./256, NRMSE_0DL_all,'k-^','LineWidth',3,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',6);  
h3=semilogy(1-NON_ZERO_COEFFS_KSVD_all./256, NRMSE_KSVD_all,'b-o','LineWidth',3,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
h4=semilogy(1-NON_ZERO_COEFFS_WL_all./128^2, NRMSE_WL_all,'m-s','LineWidth',3,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6);  
hold off
set ( gca, 'xdir', 'reverse' )


xlabel('$Sparsity$','Interpreter','latex')
% ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
ylim([1*1e-4 5*1e-2]);
xlim([0 0.5]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-4,1e-3,1e-2,1e-1]);
% set (gca, 'YTickLabel', {num2str(-4, '10^{%d}'),'',num2str(-3, '10^{%d}'),...
%     num2str(-2, '10^{%d}'),num2str(-1, '10^{%d}')},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.1:0.5);
set (gca, 'XTickLabel', {'0.0','0.1','0.2','0.3','0.4','0.5'},'FontSize',fsize)

leg =legend ([h1,h2,h3,h4],{'PCA','ODL','KSVD','Wavelet(Meyer)'},'interpreter', 'latex','location','southwest','FontSize',fsize); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on
grid on

% export_fig('./figures/Sparsity_vs_NRMSE_PCA_ODL_KSVD_WL_Meyer_patchsize04','-eps');
% close();