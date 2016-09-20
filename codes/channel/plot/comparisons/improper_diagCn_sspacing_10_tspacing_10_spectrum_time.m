% observation model: Z=I_sY+hn
% h=I_sS_sX-X
% Ny=257;Nz=288;Nt=10000;
% position of HWs in all domain:
% PIVids_y=1:Ny; PIVids_z=1:Nz;
% HWids_y=4:space_spacing:Ny; HWids_z=4:space_spacing:Nz;
% t_all=1:Nt; t_knots=5:time_spacing:Nt;

clear all; close all; clc;

Ny=257;Nz=288;
space_spacing=10; % subsample from original grid
time_spacing=10; % from 40Hz (subsampled from 200Hz) to 4Hz
t_off=5; 

yid=128; % midle of the domain and midle of 2 HWs "line"
dt=1/40; % data is at 40Hz

%%
point_id_y_fulldomain=[119,129,139]; % HWs (in y) are at 2 lines at yid=124 and 134
point_id_z_fulldomain=9:space_spacing:279; % HWs (in z) are at 4:10:284
point_id_y_MAP=point_id_y_fulldomain-3; % truncated region in MAP as idsy_all_Z=4:254
point_id_z_MAP=point_id_z_fulldomain-3; % truncated region in MAP as idsz_all_Z=14:274

nc = netcdf('/data/DNSDATA/data/grid.nc','r');
y=nc{'gridy'}(point_id_y_fulldomain,point_id_z_fulldomain);
z=nc{'gridz'}(point_id_y_fulldomain,point_id_z_fulldomain);
close(nc)

Nt=9000;
nc1=netcdf('/data/DNSDATA/github/Bayesianfusion/Ufluc_40Hz.nc','r');
U_org=nc1{'Uall'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain); 
close(nc1);

nc2 = netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_spatialspacing_10.nc','r');
U_interp_space=nc2{'Uinterp'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain);  
close(nc2);

nc3=netcdf('/data/DNSDATA/github/Bayesianfusion/Uinterp_timespacing_10.nc','r');
U_interp_time=nc3{'Uinterp'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain);  
close(nc3);

nc4=netcdf('/data/DNSDATA/github/Bayesianfusion/FusedData_40Hz_alldomain_diagCov_timespacing_10_spacespacing_10.nc','r');
U_MAP = nc4{'Zhat_all'}(1:Nt,point_id_y_MAP,point_id_z_MAP);
close(nc4); 

nc5 = netcdf('/data/DNSDATA/regression/RESULTS_PRED_RR_U_sspacing_10_tspacing_10.nc','r');
U_RR=nc5{'U_pred'}(t_off:Nt+t_off-1,point_id_y_fulldomain,point_id_z_fulldomain); 
close(nc5); 

%% Compute fft
Fs=40; window=Nt/40;
dt=1/Fs;
f=[0:window-1]'/(window*dt);

f_org=(0:Nt-1)/(Nt*dt);
f_cutoff=f_org(ceil(Nt/2))/time_spacing;

Pxx_ref=0;
Pxx_MAP=0;
Pxx_RR=0;
Pxx_interp_space=0;
Pxx_interp_time=0;
for i=1:size(U_org,2)
    for j=1:size(U_org,3)
        xref = U_org(:,i,j);        
        x_MAP =U_MAP(:,i,j);
        x_RR =U_RR(:,i,j);
        x_interp_space=U_interp_space(:,i,j);
        x_interp_time=U_interp_time(:,i,j);

        Pxx_temp = pwelch(xref,window);
        Pxx_ref=Pxx_ref+(1/(size(U_org,2)*size(U_org,3)))*Pxx_temp;

        Pxx_temp = pwelch(x_MAP,window);
        Pxx_MAP=Pxx_MAP+(1/(size(U_org,2)*size(U_org,3)))*Pxx_temp;

        Pxx_temp = pwelch(x_RR,window);
        Pxx_RR=Pxx_RR+(1/(size(U_org,2)*size(U_org,3)))*Pxx_temp;

        Pxx_temp = pwelch(x_interp_space,window);
        Pxx_interp_space=Pxx_interp_space+(1/(size(U_org,2)*size(U_org,3)))*Pxx_temp;

        Pxx_temp = pwelch(x_interp_time,window);
        Pxx_interp_time=Pxx_interp_time+(1/(size(U_org,2)*size(U_org,3)))*Pxx_temp;
    end
end

%%
fsize=25;
fname='CMU Serif';

% -5/3 line
xline=[1,3];
yline=0.3*10^-1*[exp((-5/3)*log(1)),exp((-5/3)*log(3))];

%% PLOT E(k)
h=figure;

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [400 100 1000 800]);
set(gcf, 'Color', 'w');

h1=loglog(f(2:floor(window/2)),Pxx_ref(2:end-1),'k-','LineWidth',3);
hold on;
h2=loglog(f(2:floor(window/2)),Pxx_MAP(2:end-1),'r-','LineWidth',3);
h3=loglog(f(2:floor(window/2)),Pxx_RR(2:end-1),'m-','LineWidth',3);
h4=loglog(f(2:floor(window/2)),Pxx_interp_space(2:end-1),'b-','LineWidth',2);
h5=loglog(f(2:floor(window/2)),Pxx_interp_time(2:end-1),'g-','LineWidth',3);

plot([f_cutoff,f_cutoff] ,[10^-7,4.4*10^-3],'k--','LineWidth',3)
loglog(xline,yline,'k-','LineWidth',2)
hold off

text(2,0.03,'-5/3','HorizontalAlignment','right','FontSize',fsize)

xlim([1.5*10^-1 2*10^1]);
ylim([10^-7 9*10^-2]);

xlabel('$f$', 'interpreter', 'latex'); ylabel('$E(f)$', 'interpreter', 'latex');
% xlabel('Frequency'); ylabel('Energy');

leg=legend([h1 h2 h3 h4 h5],{'Reference','Bayesian fusion','RR','$\mathbf{I}_s \mathbf{y}$','$\mathbf{I}_t \mathbf{x}$'},'interpreter', 'latex','location','northeast');
set(leg,'FontSize',fsize-4);
legend boxoff

box on
grid on


filename=strcat('./figures/improper_sspacing_10_tspacing_10_spectrum_time');
export_fig(filename,'-eps','-q101','-a4','-nocrop');
close();

%% Zoom
h=figure;

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [400 100 1000 800]);
set(gcf, 'Color', 'w');

h1=loglog(f(3:12),Pxx_ref(3:12),'k-','LineWidth',3.5);
hold on;
h2=loglog(f(3:12),Pxx_MAP(3:12),'r-','LineWidth',3.5);
h3=loglog(f(3:12),Pxx_RR(3:12),'m-','LineWidth',3.5);
h4=loglog(f(3:12),Pxx_interp_space(3:12),'b-','LineWidth',2);
h5=loglog(f(3:12),Pxx_interp_time(3:12),'g-','LineWidth',3.5);

plot([f_cutoff,f_cutoff] ,[10^-7,4.4*10^-3],'k--','LineWidth',3.5)
% loglog(xline,yline,'k-','LineWidth',2)
hold off

% text(40,20,'-5/3','HorizontalAlignment','right','FontSize',fsize+2)

xlim([min(f(3:12))*0.9 max(f(3:12))*1.1]);
ylim([2*10^-3 6*10^-2]);
set (gca, 'YTickLabel', {[]})
set (gca, 'XTickLabel', {[]})
box on

filename=strcat('./figures/improper_sspacing_10_tspacing_10_spectrum_time_zoom');
export_fig(filename,'-eps','-q101','-a4','-nocrop');
close();