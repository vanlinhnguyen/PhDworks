%% LOAD TO NETCDF
clear all; close all; clc;
Ny=257; Nz=288; Nt=10000;
%% PLANES TO LEARN THE DICTIONARIES
nc = netcdf('/data/DNSDATA/data/Ufluc_40Hz.nc','r');
u=nc{'Uall'}(1,:,:);
close(nc)

nc = netcdf('/data/DNSDATA/data/grid.nc','r');
gridz=nc{'gridz'}(:,:);
gridy=nc{'gridy'}(:,:);
close(nc)

space_spacing=10; % subsample from original grid

yid=128; % midle of the domain and midle of 2 HWs "line"

%%
nc = netcdf('/data/DNSDATA/data/grid.nc','r');
y=nc{'gridy'}(yid,1);
z=nc{'gridz'}(1,1:Nz);
dz=z(2)-z(1);
close(nc)

%% cutoff frequency
k_org=(0:(Nz-1))*2*pi/(Nz*dz);
f_cutoff1=k_org(ceil(Nz/2))/5;
f_cutoff2=k_org(ceil(Nz/2))/10;
f_cutoff3=k_org(ceil(Nz/2))/20;
%% Compute fft
N=2*Nz-2; % number of points after periodization
k=(0:N-1)*2*pi/(2*Nz*dz); % wavenumber (rad/m), even positions are consistent with k_org

nc1=netcdf('/data/DNSDATA/data/Ufluc_40Hz.nc','r');

pdft_ref=zeros(1,N);

for t=1:Nt
    xref = make_periodic(nc1{'Uall'}(t,yid,1:Nz)',z);        

    xdft_ref = fft(xref);
    pdft_temp=abs(xdft_ref).^2/(N*dz);
    pdft_ref=pdft_ref+pdft_temp;
end
close(nc1);

pdft_ref=pdft_ref./Nt;

%%
fsize=30;
fname='CMU Serif';

%% PLOT E(k)
close all;
k1=-5/3;
xline1=[6*10^0,2*10^1];
yline1=1*10^2*[exp(k1*log(xline1(1))),exp(k1*log(xline1(2)))];

k2=-4;
xline2=[3*10^1,1*10^2];
yline2=2*10^5*[exp(k2*log(xline2(1))),exp(k2*log(xline2(2)))];

k3=-8;
xline3=[1.3*10^2,2.5*10^2];
yline3=2*10^13*[exp(k3*log(xline3(1))),exp(k3*log(xline3(2)))];

h=figure(); 
set(gcf, 'Position', [400 100 1000 700]);
set(h, 'Color', 'w');

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define paper setup to print
set(gcf,'Units','normal');
set(gca,'Position',[0.15 0.15 0.8 0.8]); % [x_leftlowcorner y__leftlowcorner width height]
set(gcf,'Units','pixels');

h1=loglog(k(3:2:N/2),pdft_ref(3:2:N/2),'r-','LineWidth',2);
hold on;
loglog(xline1,yline1,'k--','LineWidth',2)
loglog(xline2,yline2,'k--','LineWidth',2)
loglog(xline3,yline3,'k--','LineWidth',2)

plot([f_cutoff1,f_cutoff1] ,[10^-8,8*10^-3],'b--','LineWidth',2)
plot([f_cutoff2,f_cutoff2] ,[10^-8,10^-1],'b--','LineWidth',2)
plot([f_cutoff3,f_cutoff3] ,[10^-8,5.7*10^-1],'b--','LineWidth',2)

hold off

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2]);
set (gca, 'YTickLabel', {num2str(-8, '10^{%d}'),'',num2str(-6, '10^{%d}'),...
    '',num2str(-4, '10^{%d}'),'',num2str(-2, '10^{%d}'),...
    '',num2str(0, '10^{%d}'),'',num2str(2, '10^{%d}')},...
    'FontSize',fsize)

text(15,10,'-5/3','HorizontalAlignment','right','FontSize',fsize-4)
text(70,1e-1,'-4','HorizontalAlignment','right','FontSize',fsize-4)
text(2.3*1e2,5*1e-5,'-8','HorizontalAlignment','right','FontSize',fsize-4)

xlim([2*10^0 3*10^2]);
ylim([10^-8 10^2]);

xlabel('$k$','interpreter','latex'); ylabel('$E(k)$','interpreter','latex');

grid on
grid minor

export_fig('./figures/spectrum_spanwise_refDNS', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close();





%% PLOT E(k)
close all;

fsize=24;
fname='CMU Serif';

u_HR_2D=squeeze(u(:,:,1));

h=figure;
set(h, 'Position', [200 200 1200 800]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.1 0.2 0.3 0.4]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

ax(1)=subplot(2,2,1,'position',[0.1 0.55 0.3 0.3]); % top left
uimagesc(gridz(1,:),gridy(:,1),flipud(u)); caxis([-0.25,0.25]); %axis equal;

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:1.5:1.5);
set (gca, 'XTickLabel', {'', '',''})
set(gca, 'YTick', 0:1:2);
set (gca, 'YTickLabel', {'0', '1','2'})
box on
% xlabel('$z/H$','interpreter', 'latex');
ylabel('$y/H$','interpreter', 'latex');

ax(2)=axes('position',[0.1 0.2 0.3 0.3]); % top right
uimagesc(gridz(1,1:10:end),gridy(1:10:end,1),flipud(u(1:10:end,1:10:end))); caxis([-0.25,0.25]); %axis equal;

ylim([min(min(gridy)), max(max(gridy))]);
xlim([min(min(gridz)), max(max(gridz))]);

set(gca, 'XTick', -1.5:1.5:1.5);
set (gca, 'XTickLabel', {'-1.5', '0.0','1.5'})
set(gca, 'YTick', 0:1:2);
set (gca, 'YTickLabel', {'0', '1','2'})
box on
xlabel('$z/H$','interpreter', 'latex');
ylabel('$y/H$','interpreter', 'latex');

cb=colorbar('northoutside');
set(cb,'position',[0.1 0.895 0.3 0.017])
set(cb, 'YTick', -0.25:0.25:0.25);
set(cb, 'YTickLabel', {'-0.25','0.00','0.25'})


% ==================== spectrum =======================
k1=-5/3;
xline1=[6*10^0,2*10^1];
yline1=1*10^2*[exp(k1*log(xline1(1))),exp(k1*log(xline1(2)))];

k2=-4;
xline2=[3*10^1,1*10^2];
yline2=2*10^5*[exp(k2*log(xline2(1))),exp(k2*log(xline2(2)))];

k3=-8;
xline3=[1.3*10^2,2.5*10^2];
yline3=2*10^13*[exp(k3*log(xline3(1))),exp(k3*log(xline3(2)))];

ax(3)=axes('position',[0.43 0.2 0.45 0.65]); % top right

h1=loglog(k(3:2:N/2),pdft_ref(3:2:N/2),'r-','LineWidth',2);
hold on;
loglog(xline1,yline1,'k--','LineWidth',2)
loglog(xline2,yline2,'k--','LineWidth',2)
loglog(xline3,yline3,'k--','LineWidth',2)

plot([f_cutoff1,f_cutoff1] ,[10^-8,8*10^-3],'b--','LineWidth',2)
plot([f_cutoff2,f_cutoff2] ,[10^-8,10^-1],'b--','LineWidth',2)
plot([f_cutoff3,f_cutoff3] ,[10^-8,5.7*10^-1],'b--','LineWidth',2)

hold off

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0,1e1,1e2]);
set (gca, 'YTickLabel', {num2str(-8, '10^{%d}'),'',num2str(-6, '10^{%d}'),...
    '',num2str(-4, '10^{%d}'),'',num2str(-2, '10^{%d}'),...
    '',num2str(0, '10^{%d}'),'',num2str(2, '10^{%d}')},...
    'FontSize',fsize, 'YAxisLocation', 'right')

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1,1e2]);
set (gca, 'XTickLabel', {num2str(0, '10^{%d}'),num2str(1, '10^{%d}'),...
    num2str(2, '10^{%d}')},'FontSize',fsize)

text(15,10,'-5/3','HorizontalAlignment','right','FontSize',fsize-4)
text(70,1e-1,'-4','HorizontalAlignment','right','FontSize',fsize-4)
text(2.3*1e2,5*1e-5,'-8','HorizontalAlignment','right','FontSize',fsize-4)

xlim([2*10^0 3*10^2]);
ylim([10^-8 10^2]);

xlabel('$k$','interpreter','latex'); 
ylabel('$E(k)$','interpreter','latex','Rotation',-90, 'Units', 'Normalized', 'Position', [1.2, 0.5, 0]);

text(80,2*1e-8,'$k_{1}$','HorizontalAlignment','right','FontSize',fsize-2,'interpreter','latex')
text(40,2*1e-8,'$k_{2}$','HorizontalAlignment','right','FontSize',fsize-2,'interpreter','latex')
text(20,2*1e-8,'$k_{3}$','HorizontalAlignment','right','FontSize',fsize-2,'interpreter','latex')

grid on
grid minor

export_fig('./figures/samplevels_spectrum_spanwise_refDNS', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close();
