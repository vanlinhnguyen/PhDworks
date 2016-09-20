clear all; close all; clc;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';

%% Global constant
space_spacing=6; % subsampling ration in space
time_spacing=8; % subsampling ration in time (from 40Hz to 4Hz)
LTHS_idt=1:time_spacing:96;
t_midplanes=5:time_spacing:LTHS_idt(end);
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

x_ref_all=nc1{'velocity_x'}(:,4:space_spacing:end,4:space_spacing:end,1:1:LTHS_idt(end));
x_interp_spatial_all=nc2{'Uinterp'}(:,4:space_spacing:end,4:space_spacing:end,1:1:LTHS_idt(end));
x_interp_temporal_all=nc3{'Uinterp'}(:,4:space_spacing:end,4:space_spacing:end,1:1:LTHS_idt(end));
x_fusion_BF_all=nc4{'Zhat_all'}(:,4:space_spacing:end,4:space_spacing:end,1:1:LTHS_idt(end));
x_fusion_LG_all=nc5{'Zhat_all'}(:,4:space_spacing:end,4:space_spacing:end,1:1:LTHS_idt(end));
x_NLM_greedy_all=permute(nc6{'x_HR_NLM_smallscales_all'}(:,1:1:LTHS_idt(end),4:space_spacing:end,4:space_spacing:end),[1,3,4,2])+x_interp_spatial_all;
x_NLM_nongreedy_all=permute(nc7{'x_rec_all'}(:,1:1:LTHS_idt(end),4:space_spacing:end,4:space_spacing:end),[1,3,4,2])+x_interp_spatial_all;
x_RR_all=nc8{'Urec'}(:,4:space_spacing:end,4:space_spacing:end,1:1:LTHS_idt(end));
x_KRR_all=nc9{'Urec'}(:,4:space_spacing:end,4:space_spacing:end,1:1:LTHS_idt(end));
close(nc1);  close(nc2);  close(nc3); close(nc4);
close(nc5);  close(nc6);  close(nc7); close(nc8); close(nc9);

%% Compute errors
NRMSE_interp_spatial=zeros(time_spacing+1,1); 
NRMSE_interp_temporal=zeros(time_spacing+1,1); 
NRMSE_fusion_BF=zeros(time_spacing+1,1); 
NRMSE_fusion_LG=zeros(time_spacing+1,1); 
NRMSE_NLM_greedy=zeros(time_spacing+1,1); 
NRMSE_NLM_nongreedy=zeros(time_spacing+1,1); 
NRMSE_RR=zeros(time_spacing+1,1); 
NRMSE_KRR=zeros(time_spacing+1,1); 

for t=1:time_spacing
        x_ref = x_ref_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_interp_spatial = x_interp_spatial_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_interp_temporal = x_interp_temporal_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_fusion_BF = x_fusion_BF_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_fusion_LG = x_fusion_LG_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_NLM_greedy = x_NLM_greedy_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_NLM_nongreedy = x_NLM_nongreedy_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_RR = x_RR_all(:,:,:,t:time_spacing:LTHS_idt(end));
        x_KRR = x_KRR_all(:,:,:,t:time_spacing:LTHS_idt(end));

        NRMSE_interp_spatial(t,1) = sqrt(sum((x_ref(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_interp_temporal(t,1) = sqrt(sum((x_ref(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_fusion_BF(t,1) = sqrt(sum((x_ref(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_fusion_LG(t,1) = sqrt(sum((x_ref(:)-x_fusion_LG(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_NLM_greedy(t,1) = sqrt(sum((x_ref(:)-x_NLM_greedy(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_NLM_nongreedy(t,1) = sqrt(sum((x_ref(:)-x_NLM_nongreedy(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_RR(t,1) = sqrt(sum((x_ref(:)-x_RR(:)).^2))/sqrt(sum((x_ref(:)).^2));
        NRMSE_KRR(t,1) = sqrt(sum((x_ref(:)-x_KRR(:)).^2))/sqrt(sum((x_ref(:)).^2));
end

NRMSE_interp_spatial(:,1)=mean(NRMSE_interp_spatial(1:end-1,1));
NRMSE_RR(end,1)=NRMSE_RR(1,1);
NRMSE_KRR(end,1)=NRMSE_KRR(1,1);
%% Plot
fsize=20;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 800 500]);
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
cc=[96/255 96/255 96/255];
hold on
h1=plot(NRMSE_interp_spatial,'b-s','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
h2=plot(NRMSE_interp_temporal,'g-o','LineWidth',1,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);  
h3=plot(NRMSE_RR,'m-*','LineWidth',1,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6);  
h4=plot(NRMSE_KRR,'m-o','LineWidth',1,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6);  
h5=plot(NRMSE_NLM_greedy,'-*','color',cc,'LineWidth',1,'MarkerEdgeColor',cc,'MarkerFaceColor',cc,'MarkerSize',6);  
h6=plot(NRMSE_NLM_nongreedy,'-s','color',cc,'LineWidth',1,'MarkerEdgeColor',cc,'MarkerFaceColor',cc,'MarkerSize',6);  
h7=plot(NRMSE_fusion_BF,'r-d','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
hold off

xlabel('$\tau/\delta t$','Interpreter','latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
ylim([0 0.6]);
xlim([1 9]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.6);
set (gca, 'YTickLabel', {'0.0','','','','0.2','','','','0.4','','','','0.6'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 1:4:9);
set (gca, 'XTickLabel', {'0','N/2M','N/M'},'FontSize',fsize)

leg =legend ([h1,h2,h3,h4,h5,h6,h7],{'$\mathbf{I}_s \mathbf{y}$','$ \mathbf{I}_t \mathbf{x}$','RR','KRR','Greedy propag','Nongreedy propag','Bayesian fusion',},'interpreter', 'latex','location','south'); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on

% print
export_fig('./figures/compare_NRMSE_time_middlepoints_sspacing6_tspacing8','-eps');
close();


