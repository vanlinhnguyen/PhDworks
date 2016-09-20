clear all; close all; clc;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';

%% Global constant
space_spacing=3; % subsampling ration in space
time_spacing=4; % subsampling ration in time (from 40Hz to 4Hz)
LTHS_idt=1:time_spacing:96;
t_midplanes=3:time_spacing:LTHS_idt(end);
%% Load data
fprintf('Compute errors ... \n ');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_LG=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_NLM_greedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_NLM_nongreedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_RR=strcat('/data/ISOTROPIC/regression/RR_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 
filename_KRR_rbf=strcat('/data/ISOTROPIC/regression/KRR_rbf_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_BF,'r');
nc5=netcdf(filename_LG,'r');
nc6=netcdf(filename_NLM_greedy,'r');
nc7=netcdf(filename_NLM_nongreedy,'r');
nc8=netcdf(filename_RR,'r');
nc9=netcdf(filename_KRR_rbf,'r');

x_ref_all=nc1{'velocity_x'}(:,:,:,t_midplanes);
x_interp_spatial_all=nc2{'Uinterp'}(:,:,:,t_midplanes);
x_interp_temporal_all=nc3{'Uinterp'}(:,:,:,t_midplanes);
x_fusion_BF_all=nc4{'Zhat_all'}(:,:,:,t_midplanes);
x_fusion_LG_all=nc5{'Zhat_all'}(:,:,:,t_midplanes);
x_NLM_greedy_all=permute(nc6{'x_HR_NLM_smallscales_all'}(:,t_midplanes,:,:),[1,3,4,2])+x_interp_spatial_all;
x_NLM_nongreedy_all=permute(nc7{'x_rec_all'}(:,t_midplanes,:,:),[1,3,4,2])+x_interp_spatial_all;
x_RR_all=nc8{'Urec'}(:,:,:,t_midplanes);
x_KRR_all=nc9{'Urec'}(:,:,:,t_midplanes);
close(nc1);  close(nc2);  close(nc3); close(nc4);
close(nc5);  close(nc6);  close(nc7); close(nc8); close(nc9);

%% Compute errors
NRMSE_interp_spatial=zeros(space_spacing+1,space_spacing+1); 
NRMSE_interp_temporal=zeros(space_spacing+1,space_spacing+1); 
NRMSE_fusion_BF=zeros(space_spacing+1,space_spacing+1); 
NRMSE_fusion_LG=zeros(space_spacing+1,space_spacing+1);
NRMSE_NLM_greedy=zeros(space_spacing+1,space_spacing+1);
NRMSE_NLM_nongreedy=zeros(space_spacing+1,space_spacing+1);
NRMSE_RR=zeros(space_spacing+1,space_spacing+1);
NRMSE_KRR=zeros(space_spacing+1,space_spacing+1);

for i=1:space_spacing+1
    for j=1:space_spacing+1
        x_ref = x_ref_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_interp_spatial = x_interp_spatial_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_interp_temporal = x_interp_temporal_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_fusion_BF = x_fusion_BF_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_fusion_LG = x_fusion_LG_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_NLM_greedy = x_NLM_greedy_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_NLM_nongreedy = x_NLM_nongreedy_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_RR = x_RR_all(:,i:space_spacing:end,j:space_spacing:end,:);
        x_KRR = x_KRR_all(:,i:space_spacing:end,j:space_spacing:end,:);
        
        NRMSE_interp_spatial(i,j) = sqrt(sum((x_ref(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_interp_temporal(i,j) = sqrt(sum((x_ref(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_fusion_BF(i,j) = sqrt(sum((x_ref(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_fusion_LG(i,j) = sqrt(sum((x_ref(:)-x_fusion_LG(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_NLM_greedy(i,j) = sqrt(sum((x_ref(:)-x_NLM_greedy(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_NLM_nongreedy(i,j) = sqrt(sum((x_ref(:)-x_NLM_nongreedy(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_RR(i,j) = sqrt(sum((x_ref(:)-x_RR(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_KRR(i,j) = sqrt(sum((x_ref(:)-x_KRR(:)).^2))/sqrt(sum((x_ref(:)).^2));
    end
end
NRMSE_interp_temporal = mean(mean(NRMSE_interp_temporal))*ones(size(NRMSE_interp_temporal));

%% Plot
[gridy,gridz]=meshgrid(linspace(0,1,space_spacing+1),linspace(0,1,space_spacing+1));

fsize=18;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 100 800 950]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% Change default colormap
colormap(h,parula(50)) % default one, but choose number of colors

% 
ax(1)=subplot(3,2,1,'position',[0.14 0.7 0.31 0.235]); % top left
hold on
uimagesc(gridy(1,:),gridz(:,1),NRMSE_interp_spatial);
ylabel('$\alpha/\Delta y$', 'interpreter', 'latex');
ylim([0,1]);
set(gca, 'YTick', 0:0.5:1);
set (gca, 'YTickLabel', {'0.0', '0.5','1.0'})

xlim([0,1]);
set(gca, 'XTick', []);

% [C,hfigc] = contour(gridz, gridy, NRMSE_interp_spatial,0.15:0.05:0.2);
% clabel(C,hfigc,'FontSize',fsize);
% set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
hold off;box on;

title1 = text('String', '$\mathbf{I}_s \mathbf{y}$', ...
        'interpreter', 'latex',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.5 1.03 0]);


ax(2)=axes('position',[0.51 0.7 0.31 0.235]); % top right
hold on
uimagesc(gridy(1,:),gridz(:,1),NRMSE_interp_temporal);
ylim([0,1]);
set(gca, 'YTick', []);

xlim([0,1]);
set(gca, 'XTick', []);
hold off;box on;

title2 = text('String', '$\mathbf{I}_t \mathbf{x}$', ...
        'interpreter', 'latex',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.5 1.03 0]);

    
ax(3)=axes('position',[0.14 0.41 0.31 0.235]); % bottom left
hold on
uimagesc(gridy(1,:),gridz(:,1),NRMSE_RR);
ylabel('$\alpha/\Delta y$', 'interpreter', 'latex');
ylim([0,1]);
set(gca, 'YTick', 0:0.5:1);
set (gca, 'YTickLabel', {'0.0', '0.5','1.0'})
xlim([0,1]);
set(gca, 'XTick', []);

set (gca, 'XTickLabel', {'0.0', '0.5','1.0'})
% [C,hfigc] = contour(gridz, gridy, NRMSE_interp_spatial,0.15:0.05:0.2);
% clabel(C,hfigc,'FontSize',fsize);
% set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
hold off;box on;
title3 = text('String', 'RR', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.5 1.03 0]);


ax(4)=axes('position',[0.51 0.41 0.31 0.235]);
hold on
uimagesc(gridy(1,:),gridz(:,1),NRMSE_KRR);

ylim([0,1]);
set(gca, 'YTick', []);
xlim([0,1]);
set(gca, 'XTick', []);

% [C,hfigc] = contour(gridz, gridy, NRMSE_interp_spatial,0.15:0.05:0.2);
% clabel(C,hfigc,'FontSize',fsize);
% set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
hold off;box on;

title4 = text('String', 'KRR', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.5 1.03 0]);

    
    
ax(5)=axes('position',[0.14 0.12 0.31 0.235]); % bottom left
hold on
uimagesc(gridy(1,:),gridz(:,1),NRMSE_NLM_nongreedy);
xlabel('$\beta/\Delta z$', 'interpreter', 'latex');
ylabel('$\alpha/\Delta y$', 'interpreter', 'latex');
ylim([0,1]);
set(gca, 'YTick', 0:0.5:1);
set (gca, 'YTickLabel', {'0.0', '0.5','1.0'})

xlim([0,1]);
set(gca, 'XTick', 0:0.5:1);
set (gca, 'XTickLabel', {'0.0', '0.5','1.0'})
% [C,hfigc] = contour(gridz, gridy, NRMSE_interp_spatial,0.15:0.05:0.2);
% clabel(C,hfigc,'FontSize',fsize);
% set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
hold off;box on;
title5 = text('String', 'Nongreedy propag', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.5 1.03 0]);


ax(6)=axes('position',[0.51 0.12 0.31 0.235]);
hold on
uimagesc(gridy(1,:),gridz(:,1),NRMSE_fusion_BF);
xlabel('$\beta/\Delta z$', 'interpreter', 'latex');
ylim([0,1]);
set(gca, 'YTick', []);

xlim([0,1]);
set(gca, 'XTick', 0:0.5:1);
set (gca, 'XTickLabel', {'0.0', '0.5','1.0'})
% [C,hfigc] = contour(gridz, gridy, NRMSE_interp_spatial,0.15:0.05:0.2);
% clabel(C,hfigc,'FontSize',fsize);
% set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
hold off;box on;

title6 = text('String', 'Bayesian fusion', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.5 1.03 0]);

clim=[0 0.2]; % Data range..
caxis(clim);
set(ax,'CLim',clim);
cb=colorbar;
set(cb,'position',[0.87 0.12 0.03 0.815]) 
set(cb, 'YTick', 0:0.1:0.2);
set (cb, 'YTickLabel', {'0.0','0.1','0.2'})

% export_fig('./figures/compare_boxin4HWs_sspacing3_tspacing4', '-a1','-q101','-eps');
% close();
