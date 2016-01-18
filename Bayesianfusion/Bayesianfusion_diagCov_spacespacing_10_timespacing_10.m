clear all; close all; clc;

%% Define file locations (to load or to save)
addpath('./funcs/');
space_spacing=10; % subsampling ration in space
time_spacing=10; % subsampling ration in time (from 40Hz to 4Hz)

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

t_off=5; % indicate the first LTHS snapshot 
LTHS_idt=t_off:time_spacing:Nt; % indices of all LTHS snapshots in time 
HTLS_idt=1:1:Nt; % indices of all HTHS snapshots in time

N=Ny*Nz; % total number of HTHS points in each snapshot
M=My*Mz; % total number of HTLS points in each snapshot

P=numel(HTLS_idt); % total number of HW snapshots
Q=numel(LTHS_idt); % total number of PIV snapshots

%% Interpolation step =====================================================
nc = netcdf(filename_grid,'r');
gridy_HS = nc{'gridy'}(:,:);
gridz_HS = nc{'gridz'}(:,:);
close(nc);

interp_space(gridy_HS, gridz_HS, HTHS_idy, HTHS_idz, HTLS_idy,HTLS_idz,filename_ref, filename_interp_space);
interp_time(LTHS_idt, filename_ref, filename_interp_time);

%% Estimating parameters for fusion step ==================================
fprintf('\n ********** START ESTIMATING FUSION PARAMS ************ \n');

%% Compute Cns 
start=tic();
nc1 = netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r'); 
hs=nc1{'Uall'}(LTHS_idt,HTHS_idy,HTHS_idz)-nc2{'Uinterp'}(LTHS_idt,HTHS_idy,HTHS_idz);
close(nc1); close(nc2);

S=squeeze(var(hs,1)); clearvars hs;
fprintf('Learning Cns: %.2f seconds.\n', toc(start));

% average in z-direction
S_ave=zeros(Ny,Nz);
offset_s=space_spacing; % avoid border effect by interpolation
for i=1:space_spacing
    S_ave_onecol=mean(S(:,i+offset_s:space_spacing:Nz-offset_s),2);
    S_ave(:,i:space_spacing:Nz)=repmat(S_ave_onecol,[1,numel(i:space_spacing:Nz)]);
end

%% Compute Cnt 
nc=netcdf(filename_grid,'r');
y_HTLS = nc{'gridy'}(HTLS_idy,1);
y_HTHS = nc{'gridy'}(HTHS_idy,1);
close(nc);

start=tic();
T_ave=zeros(time_spacing,Ny,Nz);

nc1 = netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_time,'r'); 

offset_t=time_spacing; % avoid border effect by interpolation
for t=1:time_spacing
    ht=nc1{'Uall'}((t+offset_t:time_spacing:P-offset_t)+t_off-1,HTLS_idy,HTLS_idz)-nc2{'Uinterp'}((t+offset_t:time_spacing:P-offset_t)+t_off-1,HTLS_idy,HTLS_idz);
    T=squeeze(mean(ht.^2,1)); clearvars ht;
    T=mean(T,2);
    T_interp=interp1(y_HTLS,T,y_HTHS);
    T_ave(t,:,:)=repmat(T_interp,[1,Nz]);
end
close(nc1); close(nc2);
T_ave(isnan(T_ave))=0;

fprintf('\n ********** FINISH ESTIMATING FUSION PARAMS ************ \n');

%% Fusion step ============================================================
fprintf('\n ********** START FUSION STEP ************ \n');
nc1 = netcdf(filename_interp_space,'r');
nc2 = netcdf(filename_interp_time,'r');
nc = netcdf(filename_fusion,'clobber');
nc('Nt')=0; 
nc('Ny')=Ny;
nc('Nz')=Nz;
nc{'Zhat_all'}=ncfloat('Nt','Ny','Nz');

% for t_fuse=t_off:LTHS_idt(end)
for t_fuse=t_off:1:100
    fprintf('Fusion %.4d-th snapshot.\n', t_fuse);
    
    pos_t=rem(t_fuse-t_off,time_spacing)+1;
    
    % estimate weights
    if(pos_t==1)
        Wt=ones(Ny,Nz);
        Ws=zeros(Ny,Nz);
    else
        Wt=S_ave./(S_ave+squeeze(T_ave(pos_t,:,:)));
        Ws=squeeze(T_ave(pos_t,:,:))./(S_ave+squeeze(T_ave(pos_t,:,:)));
    end
    
    % fusion
    IsY =  nc1{'Uinterp'}(t_fuse,HTHS_idy,HTHS_idz);
    ItX =  nc2{'Uinterp'}(t_fuse,HTHS_idy,HTHS_idz);

    Zhat = Wt.*ItX+Ws.*IsY;

    nc{'Zhat_all'}(t_fuse-t_off+1,:,:) = Zhat;
end 
close(nc1); close(nc2); close(nc);

fprintf('\n ********** FINISH FUSION STEP ************ \n');