clear all; close all; clc;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
Nt=37; Nx=96;
Ny=96; Nz=96;
%%  TIME: 04; SPACE: 03
space_spacing=03; % subsampling ration in space
time_spacing=04; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 04  SPACE: 03\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_fusion_linear=strcat('/data/ISOTROPIC/fusion/FusedData_linearGauss_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');
nc5=netcdf(filename_fusion_linear,'r');

x_nonoise=nc1{'velocity_x'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_interp_spatial=nc2{'Uinterp'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_interp_temporal=nc3{'Uinterp'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_MAP=nc4{'Zhat_all'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_MAP_linear=nc5{'Zhat_all'}(1:Nt,1:Nz,1:Ny,1:Nx);

close(nc1); close(nc2); close(nc3); close(nc4); close(nc5);

x_nonoise_lf = zeros(size(x_nonoise));
x_interp_spatial_lf = x_nonoise_lf;
x_interp_temporal_lf = x_nonoise_lf;
x_MAP_lf = x_nonoise_lf;
x_MAP_linear_lf=x_nonoise_lf;

for t=1:Nt
    for i=1:Nx
        x_nonoise_lf(t,:,:,i) = filter_2D(squeeze(x_nonoise(t,:,:,i)),space_spacing,space_spacing);
        x_interp_spatial_lf(t,:,:,i) = filter_2D(squeeze(x_interp_spatial(t,:,:,i)),space_spacing,space_spacing);
        x_interp_temporal_lf(t,:,:,i) = filter_2D(squeeze(x_interp_temporal(t,:,:,i)),space_spacing,space_spacing);
        x_MAP_lf(t,:,:,i) = filter_2D(squeeze(x_MAP(t,:,:,i)),space_spacing,space_spacing);
        x_MAP_linear_lf(t,:,:,i) = filter_2D(squeeze(x_MAP_linear(t,:,:,i)),space_spacing,space_spacing);
    end
end
x_nonoise_hf=x_nonoise-x_nonoise_lf;
x_interp_spatial_hf=x_interp_spatial-x_interp_spatial_lf;
x_interp_temporal_hf=x_interp_temporal-x_interp_temporal_lf;
x_MAP_hf=x_MAP-x_MAP_lf;
x_MAP_linear_hf=x_MAP_linear-x_MAP_linear_lf;

err_interp_spatial_time4_space3_lf=sqrt(sum((x_nonoise_lf(:)-x_interp_spatial_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_interp_temporal_time4_space3_lf=sqrt(sum((x_nonoise_lf(:)-x_interp_temporal_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_fusion_time4_space3_lf=sqrt(sum((x_nonoise_lf(:)-x_MAP_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_fusion_linear_time4_space3_lf=sqrt(sum((x_nonoise_lf(:)-x_MAP_linear_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));

err_interp_spatial_time4_space3_hf=sqrt(sum((x_nonoise_hf(:)-x_interp_spatial_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_interp_temporal_time4_space3_hf=sqrt(sum((x_nonoise_hf(:)-x_interp_temporal_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_fusion_time4_space3_hf=sqrt(sum((x_nonoise_hf(:)-x_MAP_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_fusion_linear_time4_space3_hf=sqrt(sum((x_nonoise_hf(:)-x_MAP_linear_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));

%%  TIME: 06; SPACE: 04
space_spacing=4; % subsampling ration in space
time_spacing=6; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 06  SPACE: 04\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_fusion_linear=strcat('/data/ISOTROPIC/fusion/FusedData_linearGauss_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');
nc5=netcdf(filename_fusion_linear,'r');

x_nonoise=nc1{'velocity_x'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_interp_spatial=nc2{'Uinterp'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_interp_temporal=nc3{'Uinterp'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_MAP=nc4{'Zhat_all'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_MAP_linear=nc5{'Zhat_all'}(1:Nt,1:Nz,1:Ny,1:Nx);

close(nc1); close(nc2); close(nc3); close(nc4); close(nc5);

x_nonoise_lf = zeros(size(x_nonoise));
x_interp_spatial_lf = x_nonoise_lf;
x_interp_temporal_lf = x_nonoise_lf;
x_MAP_lf = x_nonoise_lf;
x_MAP_linear_lf=x_nonoise_lf;

for t=1:Nt
    for i=1:Nx
        x_nonoise_lf(t,:,:,i) = filter_2D(squeeze(x_nonoise(t,:,:,i)),space_spacing,space_spacing);
        x_interp_spatial_lf(t,:,:,i) = filter_2D(squeeze(x_interp_spatial(t,:,:,i)),space_spacing,space_spacing);
        x_interp_temporal_lf(t,:,:,i) = filter_2D(squeeze(x_interp_temporal(t,:,:,i)),space_spacing,space_spacing);
        x_MAP_lf(t,:,:,i) = filter_2D(squeeze(x_MAP(t,:,:,i)),space_spacing,space_spacing);
        x_MAP_linear_lf(t,:,:,i) = filter_2D(squeeze(x_MAP_linear(t,:,:,i)),space_spacing,space_spacing);
    end
end
x_nonoise_hf=x_nonoise-x_nonoise_lf;
x_interp_spatial_hf=x_interp_spatial-x_interp_spatial_lf;
x_interp_temporal_hf=x_interp_temporal-x_interp_temporal_lf;
x_MAP_hf=x_MAP-x_MAP_lf;
x_MAP_linear_hf=x_MAP_linear-x_MAP_linear_lf;

err_interp_spatial_time6_space4_lf=sqrt(sum((x_nonoise_lf(:)-x_interp_spatial_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_interp_temporal_time6_space4_lf=sqrt(sum((x_nonoise_lf(:)-x_interp_temporal_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_fusion_time6_space4_lf=sqrt(sum((x_nonoise_lf(:)-x_MAP_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_fusion_linear_time6_space4_lf=sqrt(sum((x_nonoise_lf(:)-x_MAP_linear_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));

err_interp_spatial_time6_space4_hf=sqrt(sum((x_nonoise_hf(:)-x_interp_spatial_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_interp_temporal_time6_space4_hf=sqrt(sum((x_nonoise_hf(:)-x_interp_temporal_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_fusion_time6_space4_hf=sqrt(sum((x_nonoise_hf(:)-x_MAP_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_fusion_linear_time6_space4_hf=sqrt(sum((x_nonoise_hf(:)-x_MAP_linear_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));

%%  TIME: 08; SPACE: 06
space_spacing=6; % subsampling ration in space
time_spacing=8; % subsampling ration in time (from 40Hz to 4Hz)

fprintf('compute errors TIME: 08  SPACE: 06\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_fusion=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_fusion_linear=strcat('/data/ISOTROPIC/fusion/FusedData_linearGauss_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_fusion,'r');
nc5=netcdf(filename_fusion_linear,'r');

x_nonoise=nc1{'velocity_x'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_interp_spatial=nc2{'Uinterp'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_interp_temporal=nc3{'Uinterp'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_MAP=nc4{'Zhat_all'}(1:Nt,1:Nz,1:Ny,1:Nx);
x_MAP_linear=nc5{'Zhat_all'}(1:Nt,1:Nz,1:Ny,1:Nx);

close(nc1); close(nc2); close(nc3); close(nc4); close(nc5);

x_nonoise_lf = zeros(size(x_nonoise));
x_interp_spatial_lf = x_nonoise_lf;
x_interp_temporal_lf = x_nonoise_lf;
x_MAP_lf = x_nonoise_lf;
x_MAP_linear_lf=x_nonoise_lf;

for t=1:Nt
    for i=1:Nx
        x_nonoise_lf(t,:,:,i) = filter_2D(squeeze(x_nonoise(t,:,:,i)),space_spacing,space_spacing);
        x_interp_spatial_lf(t,:,:,i) = filter_2D(squeeze(x_interp_spatial(t,:,:,i)),space_spacing,space_spacing);
        x_interp_temporal_lf(t,:,:,i) = filter_2D(squeeze(x_interp_temporal(t,:,:,i)),space_spacing,space_spacing);
        x_MAP_lf(t,:,:,i) = filter_2D(squeeze(x_MAP(t,:,:,i)),space_spacing,space_spacing);
        x_MAP_linear_lf(t,:,:,i) = filter_2D(squeeze(x_MAP_linear(t,:,:,i)),space_spacing,space_spacing);
    end
end
x_nonoise_hf=x_nonoise-x_nonoise_lf;
x_interp_spatial_hf=x_interp_spatial-x_interp_spatial_lf;
x_interp_temporal_hf=x_interp_temporal-x_interp_temporal_lf;
x_MAP_hf=x_MAP-x_MAP_lf;
x_MAP_linear_hf=x_MAP_linear-x_MAP_linear_lf;

err_interp_spatial_time8_space6_lf=sqrt(sum((x_nonoise_lf(:)-x_interp_spatial_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_interp_temporal_time8_space6_lf=sqrt(sum((x_nonoise_lf(:)-x_interp_temporal_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_fusion_time8_space6_lf=sqrt(sum((x_nonoise_lf(:)-x_MAP_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));
err_fusion_linear_time8_space6_lf=sqrt(sum((x_nonoise_lf(:)-x_MAP_linear_lf(:)).^2))/sqrt(sum((x_nonoise_lf(:)).^2));

err_interp_spatial_time8_space6_hf=sqrt(sum((x_nonoise_hf(:)-x_interp_spatial_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_interp_temporal_time8_space6_hf=sqrt(sum((x_nonoise_hf(:)-x_interp_temporal_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_fusion_time8_space6_hf=sqrt(sum((x_nonoise_hf(:)-x_MAP_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));
err_fusion_linear_time8_space6_hf=sqrt(sum((x_nonoise_hf(:)-x_MAP_linear_hf(:)).^2))/sqrt(sum((x_nonoise_hf(:)).^2));