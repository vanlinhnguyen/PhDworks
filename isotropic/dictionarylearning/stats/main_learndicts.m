% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

%% INITIAL PARAMS
space_spacing=4;
size_l=4; 
patchsize_h = size_l*space_spacing;
dim_h=patchsize_h^2;

%% Learn dictionaries
num_patch = 12*12;
patches_HR_all = extract_patches_HR_patchsize04(num_patch);

[D_PCA, A_PCA,lambda_PCA] = learnPCA_patchesHR_patchsize04(patches_HR_all);
[D_ODL,A_ODL] = learnODL_patchesHR_patchsize04(patches_HR_all);

%% Plot dictionaries
% PCA
fig_PCA=figure(); set(fig_PCA, 'Color', 'w'); ImD_ODL=displayPatches(D_PCA); axis off;
export_fig('./figures/dictionary_PCA_patchsize04','-eps'); close();

% ODL
fig_ODL=figure(); set(fig_ODL, 'Color', 'w'); ImD_ODL=displayPatches(D_ODL); axis off;
export_fig('./figures/dictionary_ODL_patchsize04','-eps'); close();

%%
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
plot((1:16^2)./16^2,acc_ener,'r-','LineWidth',3);

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
export_fig('./figures/POD_patchsize04_accumulative_energy','-eps');
close();