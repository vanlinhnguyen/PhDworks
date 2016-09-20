%ERROR_INTERP_SPACE estimate error by spatial interpolation
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

function [err_mean,err_max, err_LF, err_HF] = error_interp_space(space_spacing, err)

% Initial params
Nt = 37;
Nx = 96;
Ny = 96;
Nz = 96;

HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

% Load reference
filename_ref='/data/PhDworks/isotropic/refdata_downsampled4.nc';
nc1=netcdf(filename_ref,'r');
x_ref_all=nc1{'velocity_x'}(:,:,:,:);
close(nc1); 

% Interpolate
left=4*space_spacing; right = 4*space_spacing-(Ny-HTLS_idy(end));
bottom=4*space_spacing; top = 4*space_spacing-(Nz-HTLS_idz(end));
[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nz+left+right,1:Ny+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

nc1=netcdf(filename_ref,'r');
x_interp_all = zeros(Nt,Nz,Ny,Nx);
for t=1:Nt
    for i=1:Nx
        % interp
        asnap_LTHS=nc1{'velocity_x'}(t,:,:,i);
        asnap_LTHS_enlarged=enlarge_2D(asnap_LTHS,left, right, bottom, top);
        
        asnap_LTLS=asnap_LTHS_enlarged(1:space_spacing:end,1:space_spacing:end);
        asnap_LTHS_interp=interp2(gridz_LS_enlarged, gridy_LS_enlarged, asnap_LTLS, gridz_HS_enlarged, gridy_HS_enlarged,'spline');
               
        x_interp_all(t,:,:,i) = asnap_LTHS_interp(bottom+1:bottom+Ny,left+1:left+Nz);
    end
end 
close(nc1);

% mean
err_mean = err(x_ref_all,x_interp_all);

% max
if mod(space_spacing,2) == 0 % don't need to check time_spacing, always even in our tests
    id_firstmid = space_spacing/2+1;
    x_ref_mid = x_ref_all(:,id_firstmid:space_spacing:end,...
        id_firstmid:space_spacing:end,:);    
    x_interp_mid = x_interp_all(:,id_firstmid:space_spacing:end,...
        id_firstmid:space_spacing:end,:);        
else
    id_firstmid = round(space_spacing/2);
    x_ref_mid = [x_ref_all(:,id_firstmid:space_spacing:end,id_firstmid:space_spacing:end,:);...
        x_ref_all(:,id_firstmid+1:space_spacing:end,id_firstmid:space_spacing:end,:);...
        x_ref_all(:,id_firstmid:space_spacing:end,id_firstmid+1:space_spacing:end,:);...
        x_ref_all(:,id_firstmid+1:space_spacing:end,id_firstmid+1:space_spacing:end,:)];
    x_interp_mid = [x_interp_all(:,id_firstmid:space_spacing:end,id_firstmid:space_spacing:end,:);...
        x_interp_all(:,id_firstmid+1:space_spacing:end,id_firstmid:space_spacing:end,:);...
        x_interp_all(:,id_firstmid:space_spacing:end,id_firstmid+1:space_spacing:end,:);...
        x_interp_all(:,id_firstmid+1:space_spacing:end,id_firstmid+1:space_spacing:end,:)];
end
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