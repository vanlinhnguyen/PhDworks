% observation model: Z=I_sY+hn
% h=I_sS_sX-X
% Ny=257;Nz=288;Nt=10000;
% position of HWs in all domain:
% PIVids_y=1:Ny; PIVids_z=1:Nz;
% HWids_y=4:space_spacing:Ny; HWids_z=4:space_spacing:Nz;
% t_all=1:Nt; t_knots=5:time_spacing:Nt;

clear all; close all; clc;

%% Global constant
Ny=257;Nz=288;Nt=10000;
space_spacing=10; % subsample from original grid
time_spacing=10; % from 40Hz (subsampled from 200Hz) to 4Hz
t_off=5; 

%% Extract interested zone
point_id_y_fulldomain=129; % HWs (in y) are at 2 lines at yid=124 and 134
point_id_z_fulldomain=149; % HWs (in z) are at 4:5:284

nc = netcdf('/data/DNSDATA/data/grid.nc','r');
gridy=nc{'gridy'}(point_id_y_fulldomain,point_id_z_fulldomain);
gridz=nc{'gridz'}(point_id_y_fulldomain,point_id_z_fulldomain);
close(nc)

%% ERRORS
% nc1=netcdf('/data/DNSDATA/github/Bayesianfusion/Ufluc_40Hz.nc','r');
% nc2 = netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_spatialspacing_10.nc','r');
% nc3=netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_timespacing_10.nc','r');
% nc4 = netcdf('/data/DNSDATA/regression/RESULTS_PRED_RR_U_sspacing_10_tspacing_10.nc','r');
% nc5=netcdf('/data/DNSDATA/github/Bayesianfusion/FusedData_40Hz_alldomain_diagCov_timespacing_10_spacespacing_10.nc','r');


Nt=9000;
nc1=netcdf('/data/DNSDATA/github/Bayesianfusion/Ufluc_40Hz.nc','r');
U_org_all=nc1{'Uall'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain); 
close(nc1);

nc2 = netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_spatialspacing_10.nc','r');
U_interp_spatial_all=nc2{'Uinterp'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain);  
close(nc2);

nc3=netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_timespacing_10.nc','r');
U_interp_temporal_all=nc3{'Uinterp'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain);  
close(nc3);

nc4=netcdf('/data/DNSDATA/github/Bayesianfusion/FusedData_40Hz_alldomain_diagCov_timespacing_10_spacespacing_10.nc','r');
U_MAP_all = nc4{'Zhat_all'}(1:Nt,point_id_y_fulldomain,point_id_z_fulldomain);
close(nc4); 

nc5 = netcdf('/data/DNSDATA/regression/RESULTS_PRED_RR_U_sspacing_10_tspacing_10.nc','r');
U_RR_all=nc5{'U_pred'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain); 
close(nc5);

err_interp_spatial = sqrt(sum((U_org_all(:)-U_interp_spatial_all(:)).^2))/sqrt(sum((U_org_all(:)).^2));
err_interp_temporal = sqrt(sum((U_org_all(:)-U_interp_temporal_all(:)).^2))/sqrt(sum((U_org_all(:)).^2));
err_RR = sqrt(sum((U_org_all(:)-U_RR_all(:)).^2))/sqrt(sum((U_org_all(:)).^2));
err_MAP = sqrt(sum((U_org_all(:)-U_MAP_all(:)).^2))/sqrt(sum((U_org_all(:)).^2));

CC_interp_spatial=sum((U_org_all-mean(U_org_all)).*(U_interp_spatial_all-mean(U_interp_spatial_all)))/sqrt(sum((U_org_all-mean(U_org_all)).^2)*sum((U_interp_spatial_all-mean(U_interp_spatial_all)).^2));
CC_interp_temporal=sum((U_org_all-mean(U_org_all)).*(U_interp_temporal_all-mean(U_interp_temporal_all)))/sqrt(sum((U_org_all-mean(U_org_all)).^2)*sum((U_interp_temporal_all-mean(U_interp_temporal_all)).^2));
CC_RR=sum((U_org_all-mean(U_org_all)).*(U_RR_all-mean(U_RR_all)))/sqrt(sum((U_org_all-mean(U_org_all)).^2)*sum((U_RR_all-mean(U_RR_all)).^2));
CC_MAP=sum((U_org_all-mean(U_org_all)).*(U_MAP_all-mean(U_MAP_all)))/sqrt(sum((U_org_all-mean(U_org_all)).^2)*sum((U_MAP_all-mean(U_MAP_all)).^2));

%%
Nt=900;

nc1=netcdf('/data/DNSDATA/github/Bayesianfusion/Ufluc_40Hz.nc','r');
U_org=nc1{'Uall'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain); 
close(nc1);

nc2 = netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_spatialspacing_10.nc','r');
U_interp_spatial=nc2{'Uinterp'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain);  
close(nc2);

nc3=netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_timespacing_10.nc','r');
U_interp_temporal=nc3{'Uinterp'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain);  
close(nc3);

nc4=netcdf('/data/DNSDATA/github/Bayesianfusion/FusedData_40Hz_alldomain_diagCov_timespacing_10_spacespacing_10.nc','r');
U_MAP = nc4{'Zhat_all'}(1:Nt,point_id_y_fulldomain,point_id_z_fulldomain);
close(nc4);

nc5 = netcdf('/data/DNSDATA/regression/RESULTS_PRED_RR_U_sspacing_10_tspacing_10.nc','r');
U_RR=nc5{'U_pred'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain); 
close(nc5);
%% Use export_fig to export large image
fsize=28;
fname='CMU Serif';

%% Join the two
ts=381:551; % zoom in region

h=figure();
% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [50 50 1400 700]);
set(gcf, 'Color', 'w');

ax(1)=subplot(3,1,1, 'Position', [0.1, 0.6, 0.85, 0.35]);
hold on
h1=plot(1:Nt,U_org,'k','LineWidth',2);
h2=plot(1:Nt,U_MAP,'r','LineWidth',1.5);
h3=plot(1:time_spacing:Nt,U_org(1:time_spacing:Nt),'ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
plot([ts(1) ts(1)],[-0.2 0.1],'k--','LineWidth',1.5);
plot([ts(end) ts(end)],[-0.2 0.1],'k--','LineWidth',1.5);

hold off
% xlabel('$t/\delta t$','interpreter', 'latex');
ylabel('Velocity');
ylim([-0.2,0.1]);
xlim([0, Nt]);

set(gca, 'XTick', 0:300:Nt);
set(gca, 'YTick', -0.2:0.1:0.1);
set (gca, 'YTickLabel', {'-0.2', '-0.1', '0.0','0.1'})
box on

leg =legend ({'Reference','Bayesian fusion','LTHS positions'},'interpreter', 'latex','location','southeast'); 

set(leg,'FontSize',fsize-2);
legend boxoff
box on
text(395,-0.15,' (to zoom)','FontWeight', 'Bold', 'FontName', fname,'FontSize', fsize);

ax(2)=axes('position',[0.1, 0.15, 0.4, 0.35]);
hold on
plot(ts,U_org(ts),'k','LineWidth',2);
plot(ts,U_MAP(ts),'r','LineWidth',3);
h4=plot(ts,U_interp_spatial(ts),'b','LineWidth',2);
h5=plot(ts,U_interp_temporal(ts),'g','LineWidth',2);
plot(ts(1:time_spacing:end),U_org(ts(1:time_spacing:end)),'ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
hold off
xlabel('$t/\delta t$','interpreter', 'latex');
ylim([-0.1,0.075]);
xlim([min(ts), max(ts)]);

set(gca, 'XTick', 400:100:500);
set(gca, 'YTick', -0.1:0.05:0.075);
set (gca, 'YTickLabel', {'-0.10', '-0.05', '0.00','0.05'})
box on

ax(3)=axes('position',[0.55, 0.15, 0.4, 0.35]);
hold on
plot(ts,U_org(ts),'k','LineWidth',2);
plot(ts,U_MAP(ts),'r','LineWidth',3);
h6=plot(ts,U_RR(ts),'m','LineWidth',2);
plot(ts(1:time_spacing:end),U_org(ts(1:time_spacing:end)),'ko','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',5);
hold off
xlabel('$t/\delta t$','interpreter', 'latex');
% ylabel('$z_s$','interpreter', 'latex');
ylim([-0.1,0.075]);
xlim([min(ts), max(ts)]);

set(gca, 'XTick', 400:100:500);
set(gca, 'YTick', -0.1:0.05:0.075);
set (gca, 'YTickLabel', {[]})
box on


leg={'Reference','Bayesian fusion','LTHS positions','$\mathbf{I}_s \mathbf{y}$','$\mathbf{I}_t \mathbf{x}$','RR'};
legendflex([h1,h2,h3,h4,h5,h6],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 5], ...
                       'nrow',1, ...
                       'fontsize',fsize-2,...
                       'Interpreter','Latex',...
                       'box','off');
filename=strcat('./figures/improper_point_spacespacing_10_timespacing_10_yid',num2str(point_id_y_fulldomain,'%.3d'),'_zid',num2str(point_id_z_fulldomain,'%.3d'));
export_fig(filename,'-eps','-q101','-a4');
close();
