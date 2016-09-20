clear all; close all; clc;

%% Global constant
Ny=257;Nz=288;Nt=10000;
space_spacing=10; % subsample from original grid
time_spacing=10; % from 40Hz (subsampled from 200Hz) to 4Hz
t_off=5; 

%% Extract interested zone
outer_id_y_fulldomain=60:200; % HWs (in y) are at 2 lines at yid=124 and 134
outer_id_z_fulldomain=24:264; % HWs (in z) are at 4:5:284
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


%% Plot together
[z_resol,y_resol,U_org_resol]=increase_resol (gridz(1,:),gridy(:,1),U_org,1);

%%
fsize=12;
fname='CMU Serif';
h=figure;
set(h, 'Position', [500 10 800 1000]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)


% 

uimagesc(z_resol,y_resol,flipud(U_org_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylabel('$y/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})

box on


title1 = text('String','$\mathbf{I}_s \mathbf{y}$', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize, ...
        'interpreter', 'latex',...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);

    
ax(2)=axes('position',[0.51 0.83 0.375 0.13]); % top right
uimagesc(z_resol,y_resol,flipud(U_interp_spatial_largescales_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', [])

box on


title2 = text('String','$\mathbf{I}_s \mathbf{y}$', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize, ...
        'interpreter', 'latex',...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);
    
    
ax(3)=axes('position',[0.1 0.66 0.375 0.13]); % top right
uimagesc(z_resol,y_resol,flipud(U_interp_temporal_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylabel('$y/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})

box on


title3 = text('String','$\mathbf{I}_t \mathbf{x}$', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize, ...
        'interpreter', 'latex',...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);

    
ax(4)=axes('position',[0.51 0.66 0.375 0.13]); % top right
uimagesc(z_resol,y_resol,flipud(U_interp_temporal_largescales_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', [])

box on


title4 = text('String','$\mathbf{I}_t \mathbf{x}$', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize, ...
        'interpreter', 'latex',...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);

ax(5)=axes('position',[0.1 0.49 0.375 0.13]); % top right
uimagesc(z_resol,y_resol,flipud(U_RR_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylabel('$y/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})

box on


title5 = text('String','Penalized LSE', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize-2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);

    
ax(6)=axes('position',[0.51 0.49 0.375 0.13]); % top right
uimagesc(z_resol,y_resol,flipud(U_RR_largescales_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', [])

box on


title6 = text('String','Penalized LSE', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize-2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);


    
ax(7)=axes('position',[0.1 0.32 0.375 0.13]); % top right
uimagesc(z_resol,y_resol,flipud(U_MAP_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylabel('$y/H$','interpreter', 'latex');

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})

box on


title7 = text('String','Fusion model', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize',fsize-2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);

    
ax(8)=axes('position',[0.51 0.32 0.375 0.13]); % top right
uimagesc(z_resol,y_resol,flipud(U_MAP_largescales_resol));
set(gca,'YDir','normal'); % to reverse y-axis

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', [])
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', [])

box on


title8 = text('String','Fusion model', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize',fsize-2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);
    
    
ax(9)=axes('position',[0.1 0.15 0.375 0.13]); % top left
uimagesc(z_resol,y_resol,flipud(U_org_resol));
set(gca,'YDir','normal'); % to reverse y-axis

xlabel('$z/H$','interpreter', 'latex');
ylabel('$y/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})
box on


title9 = text('String', 'Reference', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize-2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);


ax(10)=axes('position',[0.51 0.15 0.375 0.13]); % top left
uimagesc(z_resol,y_resol,flipud(U_org_largescales_resol));
set(gca,'YDir','normal'); % to reverse y-axis

xlabel('$z/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', {[]})
box on


title10 = text('String', 'Reference', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize-2, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);

clim=[-0.20 0.2]; % Data range..
caxis(clim);
set(ax,'CLim',clim);
cb=colorbar;
set(cb,'position',[0.92 0.15 0.025 0.81]) 

set(cb, 'YTick', -0.2:0.05:0.2);
set (cb, 'YTickLabel', {'-0.20','-0.15','-0.10','-0.05','0.00','0.05','0.10','0.15','0.20'},'FontSize',fsize,'FontName',fname)


% filename=strcat('./Figures/improper_outer_spacespacing_10_timespacing_10_subplots_all_t',num2str(ttd,'%.3d'));
% export_fig(filename,'-eps','-q101','-a4');
% close();