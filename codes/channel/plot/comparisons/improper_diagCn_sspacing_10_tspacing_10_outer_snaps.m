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

nc1=netcdf('/data/DNSDATA/github/Bayesianfusion/Ufluc_40Hz.nc','r');
U_org=nc1{'Uall'}(ttd+t_off-1,outer_id_y_fulldomain,outer_id_z_fulldomain); 
close(nc1);

nc2 = netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_spatialspacing_10.nc','r');
U_interp_spatial=nc2{'Uinterp'}(ttd+t_off-1,outer_id_y_fulldomain,outer_id_z_fulldomain);  
close(nc2);

nc3=netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_timespacing_10.nc','r');
U_interp_temporal=nc3{'Uinterp'}(ttd+t_off-1,outer_id_y_fulldomain,outer_id_z_fulldomain);  
close(nc3);

nc4=netcdf('/data/DNSDATA/github/Bayesianfusion/FusedData_40Hz_alldomain_diagCov_timespacing_10_spacespacing_10.nc','r');
U_MAP = nc4{'Zhat_all'}(ttd,outer_id_y_MAP,outer_id_z_MAP);
close(nc4);

%% Use export_fig to export large image
fsize=20;
fname='Arial';

%% ORG
% h=figure();
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% set(gcf, 'Position', [200 200 900 400]);
% set(gcf, 'Color', 'w');
% 
% % pcolor is not good for export_fig
% [z_resol,y_resol,U_org_resol]=increase_resol (gridz(1,:),gridy(:,1),U_org,5);
% uimagesc(z_resol,y_resol,flipud(U_org_resol))
% set(gca,'YDir','normal'); % to reverse y-axis
% 
% xlabel('$z/H$','interpreter', 'latex');
% ylabel('$y/H$','interpreter', 'latex');
% ylim([min(min(gridy)), max(max(gridy))]);
% xlim([min(min(gridz)), max(max(gridz))]);
% 
% set(gca, 'XTick', -1.5:0.5:1.5);
% set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
% set(gca, 'YTick', 0:0.5:2);
% set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})
% box on
% 
% caxis([-0.2, 0.2]); 
% colorbar('location','EastOutside');
% 
% cbh = findobj( gcf, 'tag', 'Colorbar');
% set(cbh,'YTick',-0.2:0.1:0.2);
% set( cbh,'YTickLabel',{'-0.2', '-0.1', '0.0','0.1','0.2'});
% 
% pos_fig=get(gca,'position');
% pos_cbh=get(cbh,'Position');
% pos_cbh(3)=0.03;
% set(cbh,'Position',pos_cbh)
% 
% % pos_fig(2)=pos_fig(2)+0.00;
% % pos_fig(4)=pos_fig(4)-0.00;
% set(gca,'position',pos_fig)
% 
% filename=strcat('./Figures/improper_outer_spacespacing_10_timespacing_10_Uorg_t',num2str(ttd,'%.3d'));
% export_fig(filename,'-eps','-q101','-a4');
% close();




%% INTERP_SPATIAL
% h=figure();
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% set(gcf, 'Position', [200 200 900 400]);
% set(gcf, 'Color', 'w');
% 
% % pcolor is not good for export_fig
% [z_resol,y_resol,U_interp_spatial_resol]=increase_resol (gridz(1,:),gridy(:,1),U_interp_spatial,5);
% uimagesc(z_resol,y_resol,flipud(U_interp_spatial_resol))
% set(gca,'YDir','normal'); % to reverse y-axis
% 
% xlabel('$z/H$','interpreter', 'latex');
% ylabel('$y/H$','interpreter', 'latex');
% ylim([min(min(gridy)), max(max(gridy))]);
% xlim([min(min(gridz)), max(max(gridz))]);
% 
% set(gca, 'XTick', -1.5:0.5:1.5);
% set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
% set(gca, 'YTick', 0:0.5:2);
% set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})
% box on
% 
% caxis([-0.2, 0.2]); 
% colorbar('location','EastOutside');
% 
% cbh = findobj( gcf, 'tag', 'Colorbar');
% set(cbh,'YTick',-0.2:0.1:0.2);
% set( cbh,'YTickLabel',{'-0.2', '-0.1', '0.0','0.1','0.2'});
% 
% pos_fig=get(gca,'position');
% pos_cbh=get(cbh,'Position');
% pos_cbh(3)=0.03;
% set(cbh,'Position',pos_cbh)
% 
% % pos_fig(2)=pos_fig(2)+0.00;
% % pos_fig(4)=pos_fig(4)-0.00;
% set(gca,'position',pos_fig)
% 
% filename=strcat('./Figures/improper_outer_spacespacing_10_timespacing_10_Uinterp_spatial_t',num2str(ttd,'%.3d'));
% export_fig(filename,'-eps','-q101','-a4');
% close();



%% INTERP_TEMPORAL
% h=figure();
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% set(gcf, 'Position', [200 200 900 400]);
% set(gcf, 'Color', 'w');
% 
% % pcolor is not good for export_fig
% [z_resol,y_resol,U_interp_temporal_resol]=increase_resol (gridz(1,:),gridy(:,1),U_interp_temporal,5);
% uimagesc(z_resol,y_resol,flipud(U_interp_temporal_resol))
% set(gca,'YDir','normal'); % to reverse y-axis
% 
% xlabel('$z/H$','interpreter', 'latex');
% ylabel('$y/H$','interpreter', 'latex');
% ylim([min(min(gridy)), max(max(gridy))]);
% xlim([min(min(gridz)), max(max(gridz))]);
% 
% set(gca, 'XTick', -1.5:0.5:1.5);
% set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
% set(gca, 'YTick', 0:0.5:2);
% set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})
% box on
% 
% caxis([-0.2, 0.2]); 
% colorbar('location','EastOutside');
% 
% cbh = findobj( gcf, 'tag', 'Colorbar');
% set(cbh,'YTick',-0.2:0.1:0.2);
% set( cbh,'YTickLabel',{'-0.2', '-0.1', '0.0','0.1','0.2'});
% 
% pos_fig=get(gca,'position');
% pos_cbh=get(cbh,'Position');
% pos_cbh(3)=0.03;
% set(cbh,'Position',pos_cbh)
% 
% % pos_fig(2)=pos_fig(2)+0.00;
% % pos_fig(4)=pos_fig(4)-0.00;
% set(gca,'position',pos_fig)
% 
% filename=strcat('./Figures/improper_outer_spacespacing_10_timespacing_10_Uinterp_temporal_t',num2str(ttd,'%.3d'));
% export_fig(filename,'-eps','-q101','-a4');
% close();






%% MAP
% h=figure();
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% set(gcf, 'Position', [200 200 900 400]);
% set(gcf, 'Color', 'w');
% 
% % pcolor is not good for export_fig
% [z_resol,y_resol,U_MAP_resol]=increase_resol (gridz(1,:),gridy(:,1),U_MAP,5);
% uimagesc(z_resol,y_resol,flipud(U_MAP_resol))
% set(gca,'YDir','normal'); % to reverse y-axis
% 
% xlabel('$z/H$','interpreter', 'latex');
% ylabel('$y/H$','interpreter', 'latex');
% ylim([min(min(gridy)), max(max(gridy))]);
% xlim([min(min(gridz)), max(max(gridz))]);
% 
% set(gca, 'XTick', -1.5:0.5:1.5);
% set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
% set(gca, 'YTick', 0:0.5:2);
% set (gca, 'YTickLabel', {'0.0', '0.5', '1.0','1.5','2.0'})
% box on
% 
% caxis([-0.2, 0.2]); 
% colorbar('location','EastOutside');
% 
% cbh = findobj( gcf, 'tag', 'Colorbar');
% set(cbh,'YTick',-0.2:0.1:0.2);
% set( cbh,'YTickLabel',{'-0.2', '-0.1', '0.0','0.1','0.2'});
% 
% pos_fig=get(gca,'position');
% pos_cbh=get(cbh,'Position');
% pos_cbh(3)=0.03;
% set(cbh,'Position',pos_cbh)
% 
% % pos_fig(2)=pos_fig(2)+0.00;
% % pos_fig(4)=pos_fig(4)-0.00;
% set(gca,'position',pos_fig)
% 
% filename=strcat('./Figures/improper_outer_spacespacing_10_timespacing_10_UMAP_t',num2str(ttd,'%.3d'));
% export_fig(filename,'-eps','-q101','-a4');
% close();


%% Plot together
z = gridz(1,:);
y = gridy(:,1);

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
uimagesc(z,y,flipud(U_org));
set(gca,'YDir','normal'); % to reverse y-axis

ylabel('$y/H$','interpreter', 'latex');
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
        'FontSize', fsize, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);


ax(3)=axes('position',[0.52 0.6 0.35 0.35]); % top right
uimagesc(z,y,flipud(U_MAP));
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
        'FontSize', fsize, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);

    
ax(2)=axes('position',[0.14 0.15 0.35 0.35]); % bottom left
uimagesc(z,y,flipud(U_interp_spatial));
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

title2 = text('String', 'Spatial interpolation', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);


ax(4)=axes('position',[0.52 0.15 0.35 0.35]);
uimagesc(z,y,flipud(U_interp_temporal));
set(gca,'YDir','normal'); % to reverse y-axis

xlabel('$z/H$','interpreter', 'latex');
ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:0.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '-1.0','-0.5', '0.0','0.5','1.0','1.5'})
set(gca, 'YTick', 0:0.5:2);
set (gca, 'YTickLabel', [])
box on

title4 = text('String', 'Temporal interpolation', ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontSize', fsize, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [0 1.85 0]);


clim=[-0.20 0.2]; % Data range..
caxis(clim);
set(ax,'CLim',clim);
cb=colorbar;
set(cb,'position',[0.89 0.15 0.02 0.80]) 

set(cb, 'YTick', -0.2:0.1:0.2);
set (cb, 'YTickLabel', {'-0.2','-0.1','0.0','0.1','0.2'},'FontSize',fsize,'FontName',fname)


filename=strcat('./figures/improper_outer_spacespacing_10_timespacing_10_subplots_t',num2str(ttd,'%.3d'));
export_fig(filename, '-a4','-q101','-eps','-painters');
close();