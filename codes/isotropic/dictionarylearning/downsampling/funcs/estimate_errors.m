%ESTIMATE_ERRORS  estimate errors of reconstructions by interpolation, SR1,
%SR2 and SR3 for all midplanes (far from LTHS planes)
%
%IN:
% err - anonymous function to estimate error (NRMSE, RMSE or MSE...)
%OUT:
% mean and std of error of reconstruction by interpolation, SR1 (HR-LR),
% SR2 (HR-LRinterp) and SR3 (features)

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [err_interp_mean,err_SR1_HRLR_mean,err_SR2_HRLRinterp_mean,err_SR3_feas_mean,...
    err_interp_std,err_SR1_HRLR_std,err_SR2_HRLRinterp_std,err_SR3_feas_std] = estimate_errors(err)

%% INITIAL PARAMS
space_spacing = 4;
time_spacing = 6;
Nh = 96;
Nt = 37;
filename_ref = '/data/PhDworks/isotropic/refdata_downsampled4.nc';
filename_interp = '/data/PhDworks/isotropic/interp_downsample_sspacing04.nc';
midplane_ids = 4:time_spacing:Nh;
n_midplanes = numel(midplane_ids);

% boder to expand the field
size_l=6; 
border = (size_l*space_spacing)/2;
left = border;
bottom=border;

%% ESTIMATE ERRORS
nc1=netcdf(filename_ref,'r');
x_ref_all = nc1{'velocity_x'}(1:Nt,:,:,midplane_ids);
close(nc1); 

nc2=netcdf(filename_interp,'r');
x_interp_all = nc2{'X_LRinterp_enlarged'}(1:Nt,bottom+1:bottom+Nh,left+1:left+Nh,midplane_ids);
close(nc2);
 
filename_SR_feas='/data/PhDworks/isotropic/dictionarylearning/downsampling/SR_couplefeatures_K1314_lambda020_allmidplanes.nc';
nc=netcdf(filename_SR_feas,'r');
x_SR_feas_all = nc{'Zhat_all'}(1:Nt,1:Nh,1:Nh,1:n_midplanes);
close(nc);

filename_SR_LRinterp='/data/PhDworks/isotropic/dictionarylearning/downsampling/SR_HR_LRinterp_K1152_lambda020_allmidplanes.nc';
nc=netcdf(filename_SR_LRinterp,'r');
x_SR_LRinterp_all = nc{'Zhat_all'}(1:Nt,1:Nh,1:Nh,1:n_midplanes);
close(nc);

filename_SR_LR='/data/PhDworks/isotropic/dictionarylearning/downsampling/SR_HR_LR_K1224_lambda010_allmidplanes.nc';
nc=netcdf(filename_SR_LR,'r');
x_SR_LR_all = nc{'Zhat_all'}(1:Nt,1:Nh,1:Nh,1:n_midplanes);
close(nc);

err_SR3_feas_all = zeros(Nt*n_midplanes,1);
err_SR2_HRLRinterp_all = zeros(Nt*n_midplanes,1);
err_SR1_HRLR_all = zeros(Nt*n_midplanes,1);
err_interp_all = zeros(Nt*n_midplanes,1);
count=0;
for t=1:Nt
    for i=1:n_midplanes
        count=count+1;
        
        x_ref=squeeze(x_ref_all(t,:,:,i));
        x_interp = squeeze(x_interp_all(t,:,:,i));
        
        x_SR_feas=squeeze(x_SR_feas_all(t,:,:,i));
        x_SR_LRinterp=squeeze(x_SR_LRinterp_all(t,:,:,i));
        x_SR_LR=squeeze(x_SR_LR_all(t,:,:,i));

        % Error 
        err_SR3_feas_all(count,1) = err(x_SR_feas,x_ref);
        err_SR2_HRLRinterp_all(count,1) = err(x_SR_LRinterp,x_ref);
        err_SR1_HRLR_all(count,1) = err(x_SR_LR,x_ref);
        err_interp_all(count,1) = err(x_interp,x_ref);
    end
end

err_SR3_feas_mean=mean(err_SR3_feas_all);
err_SR3_feas_std=std(err_SR3_feas_all);

err_SR2_HRLRinterp_mean=mean(err_SR2_HRLRinterp_all);
err_SR2_HRLRinterp_std=std(err_SR2_HRLRinterp_all);

err_SR1_HRLR_mean=mean(err_SR1_HRLR_all);
err_SR1_HRLR_std=std(err_SR1_HRLR_all);

err_interp_mean=mean(err_interp_all);
err_interp_std=std(err_interp_all);

