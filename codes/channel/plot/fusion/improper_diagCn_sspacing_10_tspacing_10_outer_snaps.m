% observation model: Z=I_sY+hn
% h=I_sS_sX-X
% Ny=257;Nz=288;Nt=10000;
% position of HWs in all domain:
% PIVids_y=1:Ny; PIVids_z=1:Nz;
% HWids_y=4:space_spacing:Ny; HWids_z=4:space_spacing:Nz;
% t_all=1:Nt; t_knots=5:time_spacing:Nt;

% clear all; close all; clc;

%% Global constant
Ny=257;Nz=288;Nt=10000;
space_spacing=10; % subsample from original grid
time_spacing=10; % from 40Hz (subsampled from 200Hz) to 4Hz
t_off=5; 

%% Extract interested zone
outer_id_y_fulldomain=60:200; % HWs (in y) are at 2 lines at yid=124 and 134
outer_id_z_fulldomain=4:284; % HWs (in z) are at 4:5:284
outer_id_y_MAP=outer_id_y_fulldomain-3; % truncated region in MAP as idsy_all_Z=4:254
outer_id_z_MAP=outer_id_z_fulldomain-3; % truncated region in MAP as idsz_all_Z=14:274

nc = netcdf('/data/DNSDATA/data/grid.nc','r');
gridy=nc{'gridy'}(outer_id_y_fulldomain,outer_id_z_fulldomain);
gridz=nc{'gridz'}(outer_id_y_fulldomain,outer_id_z_fulldomain);
close(nc)

%%
ttd=6;

nc1=netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/data/Ufluc_40Hz.nc','r');
U_org=nc1{'Uall'}(ttd+t_off-1,outer_id_y_fulldomain,outer_id_z_fulldomain); 
close(nc1);

nc2 = netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/data/variousspacings/Uinterp_spatialspacing_10.nc','r');
U_interp_spatial=nc2{'Uinterp'}(ttd+t_off-1,outer_id_y_fulldomain,outer_id_z_fulldomain);  
close(nc2);

nc3=netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/data/variousspacings/Uinterp_timespacing_10.nc','r');
U_interp_temporal=nc3{'Uinterp'}(ttd+t_off-1,outer_id_y_fulldomain,outer_id_z_fulldomain);  
close(nc3);

nc4=netcdf('/data/DNSDATA/Fusion/BayesFusion_112014/fused/variousspacings/FusedData_40Hz_alldomain_diagCn_timespacing_10_spacespacing_10_improper.nc','r');
U_MAP = nc4{'Zhat_all'}(ttd,outer_id_y_MAP,outer_id_z_MAP);
close(nc4);

%% Plot together
fsize=22;
fname='CMU Serif';

[z_resol,y_resol,U_org_resol]=increase_resol (gridz(1,:),gridy(:,1),U_org,5);
[~,~,U_interp_spatial_resol]=increase_resol (gridz(1,:),gridy(:,1),U_interp_spatial,5);
[~,~,U_interp_temporal_resol]=increase_resol (gridz(1,:),gridy(:,1),U_interp_temporal,5);
[~,~,U_MAP_resol]=increase_resol (gridz(1,:),gridy(:,1),U_MAP,5);

h=figure;
set(h, 'Position', [100 100 1800 800]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)


% 
ax(1)=subplot(2,2,1,'position',[0.14 0.6 0.35 0.35]); % top left
uimagesc(z_resol,y_resol,flipud(U_org_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylabel('$s_2/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})
box on


title1 = text('String', 'Reference', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize+2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.2 1.85 0]);


ax(3)=axes('position',[0.52 0.6 0.35 0.35]); % top right
uimagesc(z_resol,y_resol,flipud(U_MAP_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', [])

box on


title3 = text('String', 'Fusion model', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize+2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0.2 1.85 0]);

    
ax(2)=axes('position',[0.14 0.15 0.35 0.35]); % bottom left
uimagesc(z_resol,y_resol,flipud(U_interp_spatial_resol));
set(gca,'YDir','normal'); % to reverse y-axis

xlabel('$s_1/H$','interpreter', 'latex');
ylabel('$s_2/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})
box on

title2 = text('String', '$\mathbf{I}_s \mathbf{y}$', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize+2, ...
        'FontWeight', 'Bold', ...
        'interpreter', 'latex',...
        'FontName', fname, ...
        'Position', [0 1.85 0]);


ax(4)=axes('position',[0.52 0.15 0.35 0.35]);
uimagesc(z_resol,y_resol,flipud(U_interp_temporal_resol));
set(gca,'YDir','normal'); % to reverse y-axis

xlabel('$s_1/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', [])
box on

title4 = text('String', '$\mathbf{I}_t \mathbf{x}$', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize+2, ...
        'FontWeight', 'Bold', ...
        'interpreter', 'latex',...
        'FontName', fname, ...
        'Position', [0 1.85 0]);


clim=[-0.20 0.2]; % Data range..
caxis(clim);
set(ax,'CLim',clim);
cb=colorbar;
set(cb,'position',[0.89 0.15 0.02 0.80]) 

set(cb, 'YTick', -0.2:0.1:0.2);
set (cb, 'YTickLabel', {'-0.2','-0.1','0.0','0.1','0.2'},'FontSize',fsize,'FontName',fname)


filename=strcat('./Figures/improper_outer_spacespacing_10_timespacing_10_subplots_t',num2str(ttd,'%.3d'));
export_fig(filename,'-eps','-q101','-a4');
close();