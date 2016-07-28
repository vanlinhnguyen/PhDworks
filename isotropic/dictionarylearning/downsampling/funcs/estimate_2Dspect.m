%ESTIMATE_2DSPECT  estimate 2D spectra of reconstructions by interpolation,
%SR1, SR2 and SR3 for all midplanes (far from LTHS planes)
%
%IN:
% err - anonymous function to estimate error (NRMSE, RMSE or MSE...)
%OUT:
% mean and std of error of reconstruction by interpolation, SR1 (HR-LR),
% SR2 (HR-LRinterp) and SR3 (features)

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016
function [E_ref,E_interp,E_SR1_HRLR,E_SR2_HRLRinterp, E_SR3_feas, ...
    E_err_interp, E_err_SR1_HRLR, E_err_SR2_HRLRinterp, E_err_SR3_feas] = estimate_2Dspect()

%% INITIAL PARAMS
space_spacing = 4;
time_spacing = 6;
Nh = 96;
Nt = 37;
midplane_ids = 4:time_spacing:Nh;
n_midplanes = numel(midplane_ids);

% boder to expand the field
size_l=6; 
border = (size_l*space_spacing)/2;
left = border;
bottom=border;

filename_ref = '/data/PhDworks/isotropic/refdata_downsampled4.nc';
filename_interp = '/data/PhDworks/isotropic/interp_downsample_sspacing04.nc';

%% FFT
% fftn to frequency space
nc1=netcdf(filename_ref,'r');
x_ref_all = nc1{'velocity_x'}(1:Nt,:,:,midplane_ids);
close(nc1);

nc2=netcdf(filename_interp,'r');
x_interp_all = nc2{'X_LRinterp_enlarged'}(1:Nt,bottom+1:bottom+Nh,left+1:left+Nh,midplane_ids);
close(nc2);

filename_SR_feas='/data/ISOTROPIC/dictionary_learning/restart/sspacing04/downsampling/SR_couplefeatures_K1314_lambda020_allmidplanes.nc';
nc=netcdf(filename_SR_feas,'r');
x_SR_feas_all = nc{'Zhat_all'}(1:Nt,1:Nh,1:Nh,1:n_midplanes);
close(nc);

filename_SR_LRinterp='/data/ISOTROPIC/dictionary_learning/restart/sspacing04/downsampling/SR_HR_LRinterp_K1152_lambda020_allmidplanes.nc';
nc=netcdf(filename_SR_LRinterp,'r');
x_SR_LRinterp_all = nc{'Zhat_all'}(1:Nt,1:Nh,1:Nh,1:n_midplanes);
close(nc);

filename_SR_LR='/data/ISOTROPIC/dictionary_learning/restart/sspacing04/downsampling/SR_HR_LR_K1224_lambda020_allmidplanes.nc';
nc=netcdf(filename_SR_LR,'r');
x_SR_LR_all = nc{'Zhat_all'}(1:Nt,1:Nh,1:Nh,1:n_midplanes);
close(nc);

E_ref=zeros(Nh,1);
E_interp=zeros(Nh,1);
E_SR3_feas=zeros(Nh,1);
E_SR2_HRLRinterp=zeros(Nh,1);
E_SR1_HRLR=zeros(Nh,1);

E_err_interp=zeros(Nh,1);
E_err_SR3_feas=zeros(Nh,1);
E_err_SR2_HRLRinterp=zeros(Nh,1);
E_err_SR1_HRLR=zeros(Nh,1);

norm_factor = 1/(Nt*n_midplanes);
for t=1:Nt
    for i=1:n_midplanes
        x_ref=squeeze(x_ref_all(t,:,:,i));
        x_interp = squeeze(x_interp_all(t,:,:,i));

        x_SR_feas=squeeze(x_SR_feas_all(t,:,:,i));
        x_SR_LRinterp=squeeze(x_SR_LRinterp_all(t,:,:,i));
        x_SR_LR=squeeze(x_SR_LR_all(t,:,:,i));
       
        % Spect
        E = estimate_spect_2D(x_ref); 
        E_ref = E_ref + norm_factor*E;

        E = estimate_spect_2D(x_interp);
        E_interp = E_interp + norm_factor*E; 
              
        E = estimate_spect_2D(x_SR_feas);
        E_SR3_feas = E_SR3_feas + norm_factor*E;

        E = estimate_spect_2D(x_SR_LRinterp);
        E_SR2_HRLRinterp = E_SR2_HRLRinterp + norm_factor*E;
        
        E = estimate_spect_2D(x_SR_LR);
        E_SR1_HRLR = E_SR1_HRLR + norm_factor*E;
        
        % Error
        E = estimate_spect_2D(x_ref - x_interp);
        E_err_interp = E_err_interp + norm_factor*E;
               
        E = estimate_spect_2D(x_ref - x_SR_feas);
        E_err_SR3_feas = E_err_SR3_feas + norm_factor*E;  

        E = estimate_spect_2D(x_ref - x_SR_LRinterp);
        E_err_SR2_HRLRinterp = E_err_SR2_HRLRinterp + norm_factor*E;  

        E = estimate_spect_2D(x_ref - x_SR_LR);
        E_err_SR1_HRLR = E_err_SR1_HRLR + norm_factor*E;  
    end
end