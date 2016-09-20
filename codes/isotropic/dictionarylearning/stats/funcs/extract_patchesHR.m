%EXTRACT_PATCHESHR  extract coupled LR-HR patches (SR1 and
%SR2). Patch size at LR is size_l x size_l
%
%IN:
%   num_patch - number of patches extracted from each plane
%   size_l - size of LR patches (size_l x size_l)

%OUT:
%   filename - file name where data is stored

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function patches_HR_all = extract_patchesHR(num_patch, size_l)
fprintf('\n ***************** EXTRACT PATCHES LR-HR  **************** \n');

%% INITIAL PARAMS
space_spacing=4;
time_spacing=6;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';

nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;

LTHS_idt=1:time_spacing:Nh;
x_ref_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,LTHS_idt);
close(nc)

Nl=Nh/space_spacing;

patchsize_l = size_l; 
patchsize_h = size_l*space_spacing;
dim_h=patchsize_h^2;

%% EXTRACT PATCHES
num_planes=Nt*numel(LTHS_idt);

patches_HR_all = zeros(dim_h,num_patch*num_planes);

plane_id=0;
for t=1:Nt
    for i=1:numel(LTHS_idt)
        plane_id=plane_id+1;
        
        % random positions of patches
        [~,~,xrow_h,ycol_h] = random_indices(Nl,patchsize_l,space_spacing, num_patch);

        % prepare the field: load, enlarge, downsample and interpolate
        X_HR_ref = squeeze(x_ref_all(t,:,:,i));
        
        % Extract random patches: HR and LR features
        patches_HR = zeros(dim_h, num_patch);

        % loop all the patches
        for ii = 1:num_patch
            row_h = xrow_h(ii); col_h = ycol_h(ii);

            % HR
            patch_HR = X_HR_ref(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_HR(:,ii) = patch_HR(:);              
        end 
        patches_HR_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id)= patches_HR;
    end
end
