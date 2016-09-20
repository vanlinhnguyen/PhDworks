%% LOAD TO NETCDF
clear all; close all; clc;

%% PLANES TO LEARN THE DICTIONARIES
nc = netcdf('/data/DNSDATA/data/Ufluc_40Hz.nc','r');
u=nc{'Uall'}(1,:,:);
close(nc)

nc = netcdf('/data/DNSDATA/data/grid.nc','r');
gridz=nc{'gridz'}(:,:);
gridy=nc{'gridy'}(:,:);
close(nc)


%% PLOT 2D plane
fsize=20;
fname='CMU Serif';

u_HR_2D=squeeze(u(:,:,1));

h=figure;
set(h, 'Position', [200 200 1200 400]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.1 0.2 0.7 0.7]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

ax(1)=subplot(2,2,1,'position',[0.1 0.2 0.3 0.7]); % top left
% pcolor(gridz,gridy,u); shading flat; caxis([-0.25,0.25]); axis equal; 
uimagesc(gridz(1,:),gridy(:,1),flipud(u)); caxis([-0.25,0.25]); axis equal;

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:1.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '0.0','1.5'})
set(gca, 'YTick', 0:1:2);
set (gca, 'YTickLabel', {'0', '1','2'})
box on
xlabel('$z/H$','interpreter', 'latex');
ylabel('$y/H$','interpreter', 'latex');

ax(3)=axes('position',[0.43 0.2 0.3 0.7]); % top right
% pcolor(gridz(1:10:end,1:10:end),gridy(1:10:end,1:10:end),u(1:10:end,1:10:end)); shading flat; caxis([-0.25,0.25]); axis equal; 
uimagesc(gridz(1,1:10:end),gridy(1:10:end,1),flipud(u(1:10:end,1:10:end))); caxis([-0.25,0.25]); axis equal;

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:1.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '0.0','1.5'})
set(gca, 'YTick', 0:1:2);
set (gca, 'YTickLabel', {'', '',''})
box on
xlabel('$z/H$','interpreter', 'latex');

cb=colorbar;
set(cb,'position',[0.75 0.265 0.017 0.57]) 
set(cb, 'YTick', -0.25:0.25:0.25);
set(cb, 'YTickLabel', {'-0.25','0.00','0.25'})

export_fig('./figures/samplesnap_2D', '-a4','-q101','-eps','-painters');
close()