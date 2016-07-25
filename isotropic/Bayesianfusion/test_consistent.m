clear all; close all; clc;

%% Define file locations (to load or to save)
space_spacing=05; % subsampling ration in space
time_spacing=04; % subsampling ration in time (from 40Hz to 4Hz)

filename_grid='/data/DNSDATA/github/Bayesianfusion/grid.nc';
filename_ref='/data/DNSDATA/github/Bayesianfusion/Ufluc_40Hz.nc';
filename_interp_space=strcat('/data/DNSDATA/github/Bayesianfusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/DNSDATA/github/Bayesianfusion/Uinterp_timespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/DNSDATA/github/Bayesianfusion/FusedData_40Hz_alldomain_diagCov_timespacing_',num2str(space_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nz = nc('Nz').itsDimsize;
close(nc);

%% Define parameters
HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space

HTLS_idy= 4:space_spacing:254; % row indices of HTLS in space
HTLS_idz= 4:space_spacing:284; % column indices of HTLS in space

My=numel(HTLS_idy); % number of rows of HTLS in space
Mz=numel(HTLS_idz); % number of column of HTLS in space

t_off=2; % indicate the first LTHS snapshot 
LTHS_idt=t_off:time_spacing:Nt; % indices of all LTHS snapshots in time 
HTLS_idt=1:1:Nt; % indices of all HTHS snapshots in time

N=Ny*Nz; % total number of HTHS points in each snapshot
M=My*Mz; % total number of HTLS points in each snapshot

P=numel(HTLS_idt); % total number of HW snapshots
Q=numel(LTHS_idt); % total number of PIV snapshots


%% Fusion
t_fuse=1112;
pos_t=rem(t_fuse-t_off,time_spacing)+1;

nc1=netcdf('/data/DNSDATA/Fusion/BayesFusion_final/fused/variousspacings/FusedData_40Hz_alldomain_diagCn_timespacing_04_spacespacing_05_improper.nc','r');
Zhat_old=nc1{'Zhat_all'}(t_fuse-t_off+1,:,:);
close(nc1);

nc1=netcdf(filename_fusion,'r');
Zhat_new=nc1{'Zhat_all'}(t_fuse-t_off+1,4:254,4:284);
close(nc1);

nc3=netcdf(filename_ref,'r');
Zorg=nc3{'Uall'}(t_fuse,4:254,4:284);
close(nc3);

SNR_fused = 20*log10(norm(Zorg(:))/norm(Zorg(:)-Zhat_new(:)))
SNR_fused_old = 20*log10(norm(Zorg(:))/norm(Zorg(:)-Zhat_old(:)))