%ERROR_FUSION_LR estimate error of fusion results by Linear Gaussian model
%in space and time are 3x3 and 4 respectively
%
%IN:
%    space_spacing - subsampling in space
%    time_spacing - subsampling in time
%    err - anonymous function to estimate error (RMSE, MSE, NRMSE...)
%OUT:
%    err_mean - mean error
%    err_max - maximum error (from worse points)
%    err_LF - error of large scales
%    err_HF - error of small scales

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [err_mean,err_max, err_LF, err_HF] = error_fusion_LG(space_spacing, time_spacing,err)

filename_ref='/data/PhDworks/isotropic/refdata_downsampled4.nc';
filename_fusion=strcat('/data/PhDworks/isotropic/Bayesianfusion/FusedData_LG_tspacing',num2str(time_spacing,'%.1d'),'_sspacing',num2str(space_spacing,'%.1d'),'.nc');

nc1=netcdf(filename_ref,'r');
nc2=netcdf(filename_fusion,'r');

x_ref_all=nc1{'velocity_x'}(:,:,:,:);
x_fusion_all=nc2{'Zhat_all'}(:,:,:,:);

close(nc1); close(nc2);

% mean
err_mean = err(x_ref_all,x_fusion_all);

% max
if mod(space_spacing,2) == 0 % don't need to check time_spacing, always even in our tests
    id_firstmid = space_spacing/2+1;
    x_ref_mid = x_ref_all(:,id_firstmid:space_spacing:end,...
        id_firstmid:space_spacing:end,time_spacing/2+1:time_spacing:end);    
    x_fusion_mid = x_fusion_all(:,id_firstmid:space_spacing:end,...
        id_firstmid:space_spacing:end,time_spacing/2+1:time_spacing:end);        
else
    id_firstmid = round(space_spacing/2);
    x_ref_mid = [x_ref_all(:,id_firstmid:space_spacing:end,id_firstmid:space_spacing:end,time_spacing/2+1:time_spacing:end);...
        x_ref_all(:,id_firstmid+1:space_spacing:end,id_firstmid:space_spacing:end,time_spacing/2+1:time_spacing:end);...
        x_ref_all(:,id_firstmid:space_spacing:end,id_firstmid+1:space_spacing:end,time_spacing/2+1:time_spacing:end);...
        x_ref_all(:,id_firstmid+1:space_spacing:end,id_firstmid+1:space_spacing:end,time_spacing/2+1:time_spacing:end)];
    x_fusion_mid = [x_fusion_all(:,id_firstmid:space_spacing:end,id_firstmid:space_spacing:end,time_spacing/2+1:time_spacing:end);...
        x_fusion_all(:,id_firstmid+1:space_spacing:end,id_firstmid:space_spacing:end,time_spacing/2+1:time_spacing:end);...
        x_fusion_all(:,id_firstmid:space_spacing:end,id_firstmid+1:space_spacing:end,time_spacing/2+1:time_spacing:end);...
        x_fusion_all(:,id_firstmid+1:space_spacing:end,id_firstmid+1:space_spacing:end,time_spacing/2+1:time_spacing:end)];
end
err_max = err(x_ref_mid,x_fusion_mid);

% error by scales

x_ref_LF = zeros(size(x_ref_all));
x_fusion_LF = x_ref_LF;

for t=1:size(x_ref_all,1) % Nt
    for i=1:size(x_ref_all,4) % Nx
        x_ref_LF(t,:,:,i) = filter_2D(squeeze(x_ref_all(t,:,:,i)),space_spacing,space_spacing);
        x_fusion_LF(t,:,:,i) = filter_2D(squeeze(x_fusion_all(t,:,:,i)),space_spacing,space_spacing);
    end
end
x_ref_HF=x_ref_all - x_ref_LF;
x_fusion_HF=x_fusion_all - x_fusion_LF;

err_LF = err(x_ref_LF,x_fusion_LF);
err_HF = err(x_ref_HF,x_fusion_HF);