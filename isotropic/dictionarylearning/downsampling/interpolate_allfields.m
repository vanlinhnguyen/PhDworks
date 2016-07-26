% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

%% INITIAL PARAMS
space_spacing=4;
time_spacing=6;
Nt = 37;
Nh = 96;
Nl=Nh/space_spacing;

size_l=6; 
size_h=size_l*space_spacing;

% boder to expand the field
border = size_h/2;
left = border; right = border;
bottom=border; top = border;

%% Loop all blocks
fprintf('\n ****************** START INTERPOLATION ****************** \n');

filename_interp = '/data/PhDworks/isotropic/IsY_sspacing04.nc';
nc1 = netcdf(filename_interp,'clobber');

nc1('Nt')=0;
nc1('Nx')=Nh;
nc1('Nyh')=Nh+size_h;
nc1('Nzh')=Nh+size_h;
nc1{'X_LRinterp_enlarged'}=ncfloat('Nt','Nzh','Nyh','Nx');

nc1('Nyl')=(Nh+size_h)/space_spacing;
nc1('Nzl')=(Nh+size_h)/space_spacing;
nc1{'X_LR_enlarged'}=ncfloat('Nt','Nzl','Nyl','Nx');

filename_ref = '/data/PhDworks/isotropic/refdata_downsampled4.nc';
nc2 = netcdf(filename_ref,'r');

for t=1:Nt
    X_LRinterp_oneblock = zeros(Nh+size_h,Nh+size_h,Nh);
    X_LR_oneblock = zeros((Nh+size_h)/space_spacing,(Nh+size_h)/space_spacing,Nh);
    for i=1:Nh
        X_HR_org=nc2{'velocity_x'}(t,1:Nh,1:Nh,i);
        X_HR_org_enlarged=enlarge_2D(X_HR_org,left, right, bottom, top);
        
        X_LR_enlarged=resize_nguyen(X_HR_org_enlarged, 1/space_spacing,'bicubic');
        X_LR_interp_enlarged=resize_nguyen(X_LR_enlarged, space_spacing,'bicubic');
        
        X_LRinterp_oneblock(:,:,i) = X_LR_interp_enlarged;
        X_LR_oneblock(:,:,i) = X_LR_enlarged;
    end
    nc1{'X_LRinterp_enlarged'}(t,:,:,:) = X_LRinterp_oneblock;
    nc1{'X_LR_enlarged'}(t,:,:,:) = X_LR_oneblock;
end
close(nc1); close(nc2);
