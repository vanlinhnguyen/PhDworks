%% LOAD TO NETCDF
clear all; close all; clc;

%% PLANES TO LEARN THE DICTIONARIES
nc = netcdf('/data/PhDworks/isotropic/refdata_downsampled4.nc','r');
velocity_x=nc{'velocity_x'}(1,:,:,:);
close(nc)

% nc = netcdf('/data/PhDworks/isotropic/FIELD-020.nc','r');
% velocity_x=nc{'velocity_x'}(:,:,:);
% close(nc)

Nh=size(velocity_x,1);
[gridz, gridy, gridx] = meshgrid(0:1/(Nh-1):1, 0:1/(Nh-1):1, 0:1/(Nh-1):1);

%% Plot 3D isosurface
fig1=figure();
axis equal;
box on;
set(fig1,'color','w')

fv = isosurface (gridz,gridy,gridx,velocity_x,0);
p = patch(fv);
set(p,'FaceColor','red','EdgeColor','none');
camlight;
lighting gouraud;
xlabel('x');
ylabel('y');
zlabel('z');

view(45,45)
camlight(-45,45)

xlim([0,1]);
ylim([0,1]);
zlim([0,1]);
axis off; box on;
print('./figures/isosurf_uref_HR_3D','-dpng')

%% PLOT 2D plane
fsize=28;
fname='CMU Serif';

u_HR_2D=squeeze(velocity_x(:,:,1));

h=figure;
set(h, 'Position', [200 200 1600 800]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.1 0.1 0.7 0.7]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

ax(1)=subplot(2,2,1,'position',[0.1 0.1 0.3 0.7]); % top left
imagesc(u_HR_2D); caxis([-3,3]); axis off; axis equal; 


ax(3)=axes('position',[0.415 0.1 0.3 0.7]); % top right
imagesc(u_HR_2D(1:3:end,1:3:end)); caxis([-3,3]); axis off; axis equal; 

cb=colorbar;
set(cb,'position',[0.73 0.151 0.015 0.5975]) 
set(cb, 'YTick', -3:3:3);
set(cb, 'YTickLabel', {'-3','0','3'})

export_fig('./figures/uref_HR_2D', '-a1','-q101','-eps','-painters');
close()