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

idsy_all_Z=4:254; % indices of lager box in y-dir, middle at yid= 128
idsz_all_Z=4:284; % indices of lager box in z-dir


PIVids_t_all=5:time_spacing:Nt; % PIV snapshots (knots in time) 
HWsids_t_all=5:1:PIVids_t_all(end); % HW snapshots

Nt_Z=numel(HWsids_t_all); % total number of HW snapshots
Nt_X=numel(PIVids_t_all); % total number of PIV snapshots
num_blocks=numel(1:(P-1):Nt_Z)-1; % total number of blocks (P snapshots each)

%% Load data
point_id_y_fulldomain=124:134; % HWs (in y) are at 2 lines at yid=124 and 134
point_id_z_fulldomain=4:284; % HWs (in z) are at 14:10:274
point_id_y_MAP=point_id_y_fulldomain-3; % truncated region in MAP as idsy_all_Z=4:254
point_id_z_MAP=point_id_z_fulldomain-3; % truncated region in MAP as idsz_all_Z=14:274

fprintf('Compute errors ... \n ');

nc1=netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/data/Ufluc_40Hz.nc','r');
nc2 = netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/data/Uinterp_spatialspacing_10.nc','r');
nc3=netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/data/Uinterp_timespacing_050.nc','r');
nc4 = netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/LinearRegression/RESULTS_PRED_RR_U_sspacing_10_tspacing_10.nc','r');
nc5=netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/fused/variousspacings/FusedData_40Hz_alldomain_diagCn_timespacing_10_spacespacing_10_improper.nc','r');

t=1:num_blocks-1;
i_HWs=(P-1)*(t-1)+10; % if take only the worse time (between 2 PIVs)
i_HWs_2save=i_HWs-4;
x_nonoise_all=nc1{'Uall'}(i_HWs,point_id_y_fulldomain,idsz_all_Z);
x_interp_spatial_all=nc2{'Uinterp'}(i_HWs,point_id_y_fulldomain,idsz_all_Z);  
x_interp_temporal_all=nc3{'Uinterp'}(i_HWs,point_id_y_fulldomain,idsz_all_Z); 
x_RR_all=nc4{'U_pred'}(i_HWs,point_id_y_fulldomain,idsz_all_Z);    
x_MAP_all=nc5{'Zhat_all'}(i_HWs_2save,point_id_y_MAP,point_id_z_MAP);

x_nonoise=[]; x_interp_spatial=[]; x_interp_temporal=[];
x_RR=[]; x_MAP=[];
for k=point_id_z_MAP(1:space_spacing:end-space_spacing)
    point_id_z_onebox=k:k+space_spacing;
    x_nonoise=[x_nonoise; x_nonoise_all(:,:,point_id_z_onebox)];
    x_interp_spatial=[x_interp_spatial; x_interp_spatial_all(:,:,point_id_z_onebox)];
    x_interp_temporal=[x_interp_temporal; x_interp_temporal_all(:,:,point_id_z_onebox)];
    x_RR=[x_RR; x_RR_all(:,:,point_id_z_onebox)];        
    x_MAP=[x_MAP; x_MAP_all(:,:,point_id_z_onebox)]; 
end
close(nc1); close(nc2); close(nc3); close(nc4); close(nc5);

%% Compute errors
RMSE_interp_spatial=zeros(space_spacing+1,space_spacing+1); 
RMSE_interp_temporal=zeros(space_spacing+1,space_spacing+1); 
RMSE_RR=zeros(space_spacing+1,space_spacing+1); 
RMSE_MAP=zeros(space_spacing+1,space_spacing+1);

for i=1:space_spacing+1
    for j=1:space_spacing+1
        RMSE_interp_spatial(i,j)= sqrt(sum((x_nonoise(:,i,j)-x_interp_spatial(:,i,j)).^2))/sqrt(sum((x_nonoise(:,i,j)).^2));
        RMSE_interp_temporal(i,j)=sqrt(sum((x_nonoise(:,i,j)-x_interp_temporal(:,i,j)).^2))/sqrt(sum((x_nonoise(:,i,j)).^2));
        RMSE_RR(i,j)= sqrt(sum((x_nonoise(:,i,j)-x_RR(:,i,j)).^2))/sqrt(sum((x_nonoise(:,i,j)).^2));       
        RMSE_MAP(i,j)= sqrt(sum((x_nonoise(:,i,j)-x_MAP(:,i,j)).^2))/sqrt(sum((x_nonoise(:,i,j)).^2));
    end
end
    

%% Plot
[gridy,gridz]=meshgrid(linspace(0,1,11),linspace(0,1,11));

fsize=18;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 725 550]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% Change default colormap
colormap(h,parula(15)) % default one, but choose number of colors

% 
ax(1)=subplot(2,2,1,'position',[0.14 0.6 0.325 0.35]); % top left
hold on
pcolor(gridz,gridy,RMSE_interp_spatial);shading interp;
ylabel('$\alpha/\Delta y$', 'interpreter', 'latex');
ylim([0,1]);
set(gca, 'YTick', 0:0.5:1);
set (gca, 'YTickLabel', {'0.0', '0.5','1.0'})

xlim([0,1]);
set(gca, 'XTick', []);

[C,hfigc] = contour(gridz,gridy,RMSE_interp_spatial,0.3:0.05:0.5);
clabel(C,hfigc,'FontSize',fsize);
set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
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



ax(3)=axes('position',[0.52 0.6 0.325 0.35]); % top right
hold on
pcolor(gridz,gridy,RMSE_interp_temporal);shading interp;
ylim([0,1]);
set(gca, 'YTick', []);

xlim([0,1]);
set(gca, 'XTick', []);
hold off;box on;

title3 = text('String', '$\mathbf{I}_t \mathbf{x}$', ...
        'interpreter', 'latex',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.5 1.03 0]);

    
ax(2)=axes('position',[0.14 0.15 0.325 0.35]); % bottom left
hold on
pcolor(gridz,gridy,RMSE_RR);shading interp;
xlabel('$\beta/\Delta z$', 'interpreter', 'latex');
ylabel('$\alpha/\Delta y$', 'interpreter', 'latex');
ylim([0,1]);
set(gca, 'YTick', 0:0.5:1);
set (gca, 'YTickLabel', {'0.0', '0.5','1.0'})

xlim([0,1]);
set(gca, 'XTick', 0:0.5:1);
set (gca, 'XTickLabel', {'0.0', '0.5','1.0'})
[C,hfigc] = contour(gridz,gridy,RMSE_RR,0.3:0.05:0.5);
clabel(C,hfigc,'FontSize',fsize);
set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
hold off;box on;
title2 = text('String', 'Penalized LSE', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.6 1.03 0]);


ax(4)=axes('position',[0.52 0.15 0.325 0.35]);
hold on
pcolor(gridz,gridy,RMSE_MAP);shading interp;
xlabel('$\beta/\Delta z$', 'interpreter', 'latex');
ylim([0,1]);
set(gca, 'YTick', []);

xlim([0,1]);
set(gca, 'XTick', 0:0.5:1);
set (gca, 'XTickLabel', {'0.0', '0.5','1.0'})
[C,hfigc] = contour(gridz,gridy,RMSE_MAP,0.3:0.05:0.5);
clabel(C,hfigc,'FontSize',fsize);
set(hfigc,'LineWidth',0.25,'Color', [0.5 0.5 0.5]);
hold off;box on;

title4 = text('String', 'Fusion model', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+5, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.6 1.03 0]);


clim=[0 0.5]; % Data range..
caxis(clim);
set(ax,'CLim',clim);
cb=colorbar;
set(cb,'position',[0.885 0.15 0.03 0.80]) 

set(cb, 'YTick', 0:0.1:0.5);
set (cb, 'YTickLabel', {'0.0','0.1','0.2','0.3','0.4','0.5'})

export_fig('./Figures/error_MAP_newOM1_boxin4HWs_outer','-eps');
close();
