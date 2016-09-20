N_HR=384; % Number of grid-point in streamwise direction
k_max=N_HR/2;
k_max_real=2*k_max/3;

%% Load 3D field
nc=netcdf('/data/ISOTROPIC/data/FIELD-020.nc','r'); 
velocity_x_resolved=nc{'velocity_x'}(:,:,:);
close(nc)

%% 3D spectra of full resolution
% fftn to frequency space
F = fftn(velocity_x_resolved);

k_2D_HR=zeros(N_HR,1);
E_HR=zeros(N_HR,1);
k_1D_HR=[0:N_HR/2 -N_HR/2+1:1:-1];

E_ref=zeros(N_HR,1);
for m=1:N_HR
    E_ref = E_ref + 1/N_HR*estimate_spect_2D(velocity_x_resolved(:,:,m));     
end


%% Plot
close all;

k1=-5/3;
xline1=[5,10];
yline1=3.5*[exp(k1*log(xline1(1))),exp(k1*log(xline1(2)))];

k2=-11/3;
xline2=[15,30];
yline2=500*[exp(k2*log(xline2(1))),exp(k2*log(xline2(2)))];

k3=-8;
xline3=[60,100];
yline3=1e10*[exp(k3*log(xline3(1))),exp(k3*log(xline3(2)))];

%
fsize=30;
fname='CMU Serif';
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

nc=netcdf('/data/ISOTROPIC/data/STAT.nc','r');
energy_spectra = nc{'energy_spectra'}(:,:);
close(nc)

fig1=figure(); 
set(gcf, 'Position', [400 100 1000 700]);
set(fig1, 'Color', 'w');

loglog(2:k_max_real,E_ref(2:k_max_real),'r-','LineWidth',2)
hold on
plot([k_max_real/4,k_max_real/4] ,[1*10^-8,6*10^-4],'b--','LineWidth',2)
plot([k_max_real/(3*4),k_max_real/(3*4)] ,[1*10^-8,3*10^-2],'b--','LineWidth',2)
plot([k_max_real/(4*4),k_max_real/(4*4)] ,[1*10^-8,5*10^-2],'b--','LineWidth',2)
plot([k_max_real/(6*4),k_max_real/(6*4)] ,[1*10^-8,1.14*10^-1],'b--','LineWidth',2)
loglog(xline1,yline1,'k--','LineWidth',2)
loglog(xline2,yline2,'k--','LineWidth',2)
loglog(xline3,yline3,'k--','LineWidth',2)
hold off
xlim([2,128]);
ylim([2*1e-8,1e0]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-8,1e-6,1e-4,1e-2,1e0]);
% set (gca, 'YTickLabel', {'','','','',''},'FontSize',fsize)

text(9,0.29,'$-5/3$','HorizontalAlignment','right','FontSize',fsize-4,'interpreter','latex');
text(32,0.015,'$-11/3$','HorizontalAlignment','right','FontSize',fsize-4,'interpreter','latex');
text(90,2*1e-5,'$-8$','HorizontalAlignment','right','FontSize',fsize-4,'interpreter','latex');

text(48,3.6*1e-8,'$k_{max}$','HorizontalAlignment','right','FontSize',fsize-2,'interpreter','latex')
text(13,3.6*1e-8,'$k_{1}$','HorizontalAlignment','right','FontSize',fsize-2,'interpreter','latex')
text(9.7,3.6*1e-8,'$k_{2}$','HorizontalAlignment','right','FontSize',fsize-2,'interpreter','latex')
text(6.4,3.6*1e-8,'$k_{3}$','HorizontalAlignment','right','FontSize',fsize-2,'interpreter','latex')

xlabel('$k$','interpreter','latex'); ylabel('$E(k)$','interpreter','latex')

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1,1e2]);
set (gca, 'XTickLabel', {num2str(0, '10^{%d}'),num2str(1, '10^{%d}'),...
    num2str(2, '10^{%d}')},'FontSize',fsize)

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);
set (gca, 'YTickLabel', {num2str(-8, '10^{%d}'),'',num2str(-6, '10^{%d}'),...
    '',num2str(-4, '10^{%d}'),'',num2str(-2, '10^{%d}'),...
    '',num2str(0, '10^{%d}')},...
    'FontSize',fsize)

grid on
grid minor

export_fig('./figures/spectrum2d_DNS', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close();