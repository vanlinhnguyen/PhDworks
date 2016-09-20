% observation model: Z=I_sY+hn
% h=I_sS_sX-X
% Ny=257;Nz=288;Nt=5000;
% position of HWs in all domain:
% PIVids_y=1:Ny; PIVids_z=1:Nz;
% HWids_y=4:10:Ny; HWids_z=4:10:Nz;

clear all; close all; clc;

%% Global constant
Ny=257;Nz=288;Nt=10000;
space_spacing=10; % subsample from original grid
time_spacing=10; % from 40Hz (subsampled from 200Hz) to 4Hz

%% Extract interested zone
P=11; % number of HWs snapshots
Q=2;  % number of PIV snapshots
PIVids_t=1:time_spacing:P; % time indices of PIV snapshots

PIVids_t_all=5:time_spacing:Nt; % PIV snapshots (knots in time) 
HWsids_t_all=5:1:PIVids_t_all(end); % HW snapshots

Nt_Z=numel(HWsids_t_all); % total number of HW snapshots
Nt_X=numel(PIVids_t_all); % total number of PIV snapshots
num_blocks=numel(1:(P-1):Nt_Z)-1; % total number of blocks (P snapshots each)

%% Compute error
% point_id_y_fulldomain=129; % HWs (in y) are at 2 lines at yid=124 and 134
% point_id_z_fulldomain=9:10:279; % HWs (in z) are at 4:10:284
% 
% fprintf('Compute errors ... \n ');
% 
% nc1=netcdf('/data/DNSDATA/github/Bayesianfusion/Ufluc_40Hz.nc','r');
% nc2 = netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_spatialspacing_10.nc','r');
% nc3=netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_timespacing_10.nc','r');
% nc4 = netcdf('/data/DNSDATA/regression/RESULTS_PRED_RR_U_sspacing_10_tspacing_10.nc','r');
% nc5=netcdf('/data/DNSDATA/github/Bayesianfusion/FusedData_40Hz_alldomain_diagCov_timespacing_10_spacespacing_10.nc','r');
% 
% x_nonoise=[]; x_interp_spatial=[]; x_interp_temporal=[]; 
% x_RR=[]; x_MAP_diagCn=[];
% for t=1:num_blocks-1
%     t
%     i_HWs=(P-1)*(t-1)+5:(P-1)*t+5;
%     i_HWs_2save=i_HWs-4; 
%     
%     x_nonoise=[x_nonoise nc1{'Uall'}(i_HWs,point_id_y_fulldomain,point_id_z_fulldomain)];
%     x_interp_spatial=[x_interp_spatial nc2{'Uinterp'}(i_HWs,point_id_y_fulldomain,point_id_z_fulldomain)];  
%     x_interp_temporal=[x_interp_temporal nc3{'Uinterp'}(i_HWs,point_id_y_fulldomain,point_id_z_fulldomain)]; 
%     x_RR=[x_RR nc4{'U_pred'}(i_HWs,point_id_y_fulldomain,point_id_z_fulldomain)];     
%     x_MAP_diagCn= [x_MAP_diagCn nc5{'Zhat_all'}(i_HWs_2save,point_id_y_fulldomain,point_id_z_fulldomain)];
% end
% close(nc1); close(nc2); close(nc3); close(nc4); close(nc5); 
% 
% RMSE_interp_spatial=zeros(time_spacing+1,1);
% RMSE_interp_temporal=zeros(time_spacing+1,1);
% RMSE_RR=zeros(time_spacing+1,1);
% RMSE_MAP_diagCn=zeros(time_spacing+1,1);
% for tau=1:time_spacing+1
%     RMSE_interp_spatial(tau,1)=sqrt(sum((x_nonoise(tau,:)-x_interp_spatial(tau,:)).^2))/sqrt(sum((x_nonoise(tau,:)).^2));
%     RMSE_interp_temporal(tau,1)=sqrt(sum((x_nonoise(tau,:)-x_interp_temporal(tau,:)).^2))/sqrt(sum((x_nonoise(tau,:)).^2));
%     RMSE_RR(tau,1)=sqrt(sum((x_nonoise(tau,:)-x_RR(tau,:)).^2))/sqrt(sum((x_nonoise(tau,:)).^2));
%     RMSE_MAP_diagCn(tau,1)=sqrt(sum((x_nonoise(tau,:)-x_MAP_diagCn(tau,:)).^2))/sqrt(sum((x_nonoise(tau,:)).^2));
% end
% save('improper_diagCn_sspacing_10_tspacing_10_middlepoint_outer.mat','RMSE_interp_spatial','RMSE_interp_temporal','RMSE_RR','RMSE_MAP_diagCn')

%% Plot
load improper_diagCn_sspacing_10_tspacing_10_middlepoint_outer.mat
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
hold on
plot(RMSE_interp_spatial,'b-s','LineWidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',6);  
plot(RMSE_interp_temporal,'g-o','LineWidth',1,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6);  
plot(RMSE_RR,'m-*','LineWidth',1,'MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6);
plot(RMSE_MAP_diagCn,'r-d','LineWidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',6);
hold off


xlabel('$\tau/\delta t$','Interpreter','latex')
ylabel('$\epsilon(\tau)$', 'interpreter', 'latex')
ylim([0 0.6]);
xlim([1 11]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.6);
set (gca, 'YTickLabel', {'0.0','','','','0.2','','','','0.4','','','','0.6'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 1:5:11);
set (gca, 'XTickLabel', {'0','N/2M','N/M'},'FontSize',fsize)

leg =legend ({'$\mathbf{I}_s \mathbf{y}$','$ \mathbf{I}_t \mathbf{x}$','RR','Bayesian fusion',},'interpreter', 'latex','location','south'); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on

% print
export_fig('./figures/error_MAP_newOM1_middlepoint_outer','-eps');
close();


