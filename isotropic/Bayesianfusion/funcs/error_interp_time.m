%ERROR_INTERP_TIME estimate error by time interpolation
%
%IN:
%    space_spacing - subsampling in space
%    err - anonymous function to estimate error (RMSE, MSE, NRMSE...)
%OUT:
%    err_mean - mean error
%    err_max - maximum error (from worse points)
%    err_LF - error of large scales
%    err_HF - error of small scales
% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [err_mean,err_max, err_LF, err_HF] = error_interp_time(time_spacing, err)

% Initial params
Nt = 37;
Nx = 96;
Ny = 96;
Nz = 96;
LTHS_idt=1:time_spacing:Nx;

% Load reference
filename_ref='/data/PhDworks/isotropic/refdata_downsampled4.nc';
nc1=netcdf(filename_ref,'r');
x_ref_all = nc1{'velocity_x'}(:,:,:,:);
close(nc1); 

% Interpolate
left=4*time_spacing; right = 4*time_spacing-(Nx-LTHS_idt(end));
HTHS_idt_enlarged = 1:Nx+left+right;
LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);

nc1=netcdf(filename_ref,'r');

x_interp_all = zeros(Nt,Nz,Ny,Nx);

for t=1:Nt
    PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
    PIV_sampled=cat(3, PIV_sampled(:,:,Nx-left+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:right));
    
    PIV_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
    x_interp_all(t,:,:,:)=permute(PIV_interp(left+1:left+Nx,:,:),[2,3,1]);
end
close(nc1);

% mean
err_mean = err(x_ref_all,x_interp_all);

% max
x_ref_mid = x_ref_all(:,:,:,time_spacing/2+1:time_spacing:end);    
x_interp_mid = x_interp_all(:,:,:,time_spacing/2+1:time_spacing:end);        
err_max = err(x_ref_mid,x_interp_mid);

% error by scales

x_ref_LF = zeros(size(x_ref_all));
x_interp_LF = x_ref_LF;

for t=1:size(x_ref_all,1) % Nt
    for i=1:size(x_ref_all,4) % Nx
        x_ref_LF(t,:,:,i) = filter_2D(squeeze(x_ref_all(t,:,:,i)),space_spacing,space_spacing);
        x_interp_LF(t,:,:,i) = filter_2D(squeeze(x_interp_all(t,:,:,i)),space_spacing,space_spacing);
    end
end
x_ref_HF=x_ref_all - x_ref_LF;
x_interp_HF=x_interp_all - x_interp_LF;

err_LF = err(x_ref_LF,x_interp_LF);
err_HF = err(x_ref_HF,x_interp_HF);