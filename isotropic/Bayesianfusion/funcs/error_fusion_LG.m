%BAYESIANFUSION_DIAGCOV_SSPACING3_TSPACING4 fusion with subsampling ratios
%in space and time are 3x3 and 4 respectively
%
%IN:
%    space_spacing - subsampling in space
%    time_spacing - subsampling in time

%OUT:

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [err_mean, err_max] = error_fusion_LG(space_spacing, time_spacing,err)

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