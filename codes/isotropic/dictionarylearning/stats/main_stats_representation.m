% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

close all; clc; 
addpath('./funcs/')

%% INITIAL PARAMS
space_spacing = 4;
size_l = 4;
patchsize_h = size_l*space_spacing;
num_patch = 8*8;

filename_dicts=strcat('/data/PhDworks/isotropic/dictionarylearning/stats/dictstats_PCA_ODL_KSVD_patchsize',num2str(size_l,'%.1d'),'.mat');
filename_sparsity_error = strcat('/data/PhDworks/isotropic/dictionarylearning/stats/sparsity_vs_error_patchsize',num2str(size_l,'%.1d'),'.mat');

%% Learn dictionaries
% Extract patches
patches_HR_all = extract_patchesHR(num_patch, size_l);
dim_h = size(patches_HR_all,1);

% PCA 
[D_PCA, A_PCA,lambda_PCA] = learnPCA(patches_HR_all);

% ODL
params_ODL.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params_ODL.K=2*dim_h; 
params_ODL.lambda = 0.05;
params_ODL.lambda2 = 0; 
params_ODL.numThreads = 4; % number of threads
params_ODL.iter = 1000;  % max number of iterations.
[D_ODL,A_ODL] = learnODL(patches_HR_all,params_ODL);

% KSVD
params_KSVD.dictsize = 2*dim_h;
params_KSVD.Tdata = round(nnz(A_ODL)/size(A_ODL,2)); % maximal sparsity: TBD
params_KSVD.iternum = 30;
params_KSVD.memusage = 'high';
[D_KSVD, A_KSVD] = learnKSVD(patches_HR_all, params_KSVD);

% save to file
save(filename_dicts,'D_PCA','A_PCA','lambda_PCA','D_ODL','A_ODL','D_KSVD','A_KSVD');

%% Estimate Sparsity vs Error
num_snaps = 5;
[sparsity_PCA, sparsity_WL, sparsity_KSVD, sparsity_ODL, NRMSE_PCA_all, ...
    NRMSE_WL_all, NRMSE_0DL_all, NRMSE_KSVD_all] = sparsity_vs_error(filename_dicts, size_l, num_snaps);

save(filename_sparsity_error,'sparsity_PCA','sparsity_WL','sparsity_KSVD','sparsity_ODL',...
    'NRMSE_PCA_all','NRMSE_WL_all','NRMSE_0DL_all','NRMSE_KSVD_all');

%% Plot dictionaries
load(filename_dicts);

% PCA 
filename_fig = strcat('./figures/dictionary_PCA_patchsize',num2str(size_l,'%.1d'));
fig_PCA=figure(); set(fig_PCA, 'Color', 'w'); ImD_ODL=displayPatches(D_PCA); axis off;
export_fig(filename_fig,'-eps'); close();

% ODL
filename_fig = strcat('./figures/dictionary_ODL_patchsize',num2str(size_l,'%.1d'));
[~,ids]=sort(sum(A_ODL.^2,2),'descend');
fig_ODL=figure(); set(fig_ODL, 'Color', 'w'); ImD_ODL=displayPatches(D_ODL(:,ids(1:2:end))); axis off;
export_fig(filename_fig,'-eps'); close();

filename_fig = strcat('./figures/dictionary_ODL_patchsize',num2str(size_l,'%.1d'),'_full');
[~,ids]=sort(sum(A_ODL.^2,2),'descend');
fig_ODL=figure(); set(fig_ODL, 'Color', 'w'); ImD_ODL=displayPatches(D_ODL(:,ids(1:1:end))); axis off;
export_fig(filename_fig,'-eps'); close();

% KSVD
[~,ids]=sort(sum(A_KSVD.^2,2),'descend');

filename_fig = strcat('./figures/dictionary_KSVD_patchsize',num2str(size_l,'%.1d'));
fig_KSVD=figure(); set(fig_KSVD, 'Color', 'w'); ImD_KSVD=displayPatches(D_KSVD(:,ids(1:2:end))); axis off;
export_fig(filename_fig,'-eps'); close();

filename_fig = strcat('./figures/dictionary_KSVD_patchsize',num2str(size_l,'%.1d'),'_full');
fig_KSVD=figure(); set(fig_KSVD, 'Color', 'w'); ImD_KSVD=displayPatches(D_KSVD(:,ids(1:1:end))); axis off;
export_fig(filename_fig,'-eps'); close();


%% Accumulative energy of PCA
acc_ener=(cumsum(lambda_PCA)./sum(lambda_PCA));

fsize=20;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 800 600]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.15 0.2 0.725 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

% plot here
plot((1:patchsize_h^2)./patchsize_h^2,acc_ener,'r-','LineWidth',3);

xlabel('Sparsity')
ylabel('$\displaystyle \frac{\sum_{i=1}^{M} \left(\lambda_i\right)}{\sum_{i=1}^{D_h} \left(\lambda_i\right)}$','Interpreter','latex')

% ylabel('Accumulative variance')
ylim([0 1]);
% xlim([0 600]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.2:1);
set (gca, 'YTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.2:1);
set (gca, 'XTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)

box on
grid on

filename_fig = strcat('./figures/POD_accumulative_energy_patchsize',num2str(size_l,'%.1d'));
export_fig(filename_fig,'-eps');
close();


%% Plot sparsity vs error
load(filename_sparsity_error);

fsize=26;
fname='CMU Serif';

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
ax(1)=subplot(1,2,1,'position',[0.1 0.15 0.35 0.7]); % top left

hold on
h1=plot(sparsity_PCA, NRMSE_PCA_all,'r-s','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);  
h4=plot(sparsity_WL, NRMSE_WL_all,'m-s','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);  
h3=plot(sparsity_KSVD, NRMSE_KSVD_all,'b-o','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  
h2=plot(sparsity_ODL, NRMSE_0DL_all,'k-^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);  
hold off
set(ax(1),'XDir','reverse')

xlabel('$Sparsity$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
ylim([0 0.25]);
xlim([0 1]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.25);
set (gca, 'YTickLabel', {'0.0','','0.1','','0.2',''},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.2:1);
set (gca, 'XTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)

leg =legend ([h1,h2,h3,h4],{'PCA','ODL','KSVD','Wavelet'},'interpreter', 'latex','location','northeast','FontSize',fsize); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on
grid on

ax(2)=axes('position',[0.55 0.15 0.35 0.7]); % top right
% plot here
h1=semilogy(sparsity_PCA, NRMSE_PCA_all,'r-s','LineWidth',2,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',5);  
hold on
h4=semilogy(sparsity_WL, NRMSE_WL_all,'m-s','LineWidth',2,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',5);  
h3=semilogy(sparsity_KSVD, NRMSE_KSVD_all,'b-o','LineWidth',2,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',5);  
h2=semilogy(sparsity_ODL, NRMSE_0DL_all,'k-^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',5);  
hold off
set(ax(2),'XDir','reverse')

xlabel('$Sparsity$','Interpreter','latex')
ylim([1e-5 3*1e-2]);
xlim([0 0.5]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1]);
% set (gca, 'YTickLabel', {num2str(-4, '10^{%d}'),'',num2str(-3, '10^{%d}'),...
%     num2str(-2, '10^{%d}'),num2str(-1, '10^{%d}')},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.05:0.5);
set (gca, 'XTickLabel', {'0.0','','0.1','','0.2','','0.3','','0.4','','0.5'},'FontSize',fsize)

leg =legend ([h1,h2,h3,h4],{'PCA','ODL','KSVD','Wavelet'},'interpreter', 'latex','location','southwest','FontSize',fsize); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on
grid on

filename_fig = strcat('./figures/Sparsity_vs_NRMSE_PCA_ODL_KSVD_WL_Dau_patchsize',num2str(size_l,'%.1d'));
export_fig(filename_fig,'-eps');
close();
