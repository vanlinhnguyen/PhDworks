%EXTRACT_PATCHES_LRHR_PATCHSIZE06  extract coupled LR-HR patches (SR1 and
%SR2). Patch size at LR is 6 x 6
%
%IN:
%   num_patch - number of patches extracted from each plane


%OUT:
%   filename - file name where data is stored

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function filename = extract_patches_LRHR_patchsize06(num_patch)
fprintf('\n ***************** EXTRACT PATCHES LR-HR  **************** \n');

%% INITIAL PARAMS
space_spacing=4;
time_spacing=6;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_interp = '/data/PhDworks/isotropic/interp_downsample_sspacing04.nc';

nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nx').itsDimsize;

LTHS_idt=1:time_spacing:Nh;
x_ref_all=nc{'velocity_x'}(1:Nt,1:Nh,1:Nh,LTHS_idt);
close(nc)

Nl=Nh/space_spacing;

size_l=6; 
patchsize_l = size_l; 
patchsize_h = size_l*space_spacing;

dim_h=patchsize_h^2;
dim_l=patchsize_l^2;

% boder to expand the field
border=patchsize_h/2;
left = border;
bottom = border;

nc2=netcdf(filename_interp,'r');
X_LR_interp_all = nc2{'X_LRinterp_enlarged'}(1:Nt,bottom+1:bottom+Nh,left+1:left+Nh,LTHS_idt);
X_LR_all = nc2{'X_LR_enlarged'}(1:Nt,bottom/space_spacing+1:bottom/space_spacing+Nl,left/space_spacing+1:left/space_spacing+Nl,LTHS_idt);
close(nc2);
%% EXTRACT PATCHES
num_planes=Nt*numel(LTHS_idt);

patches_LR_all = zeros(dim_l,num_patch*num_planes);
patches_LRinterp_all = zeros(dim_h,num_patch*num_planes);
patches_HR_all = zeros(dim_h,num_patch*num_planes);

plane_id=0;
for t=1:Nt
    t
    for i=1:numel(LTHS_idt)
        plane_id=plane_id+1;
        
        % random positions of patches
        [xrow_l,ycol_l,xrow_h,ycol_h] = random_indices(Nl,patchsize_l,space_spacing, num_patch);

        % prepare the field: load, enlarge, downsample and interpolate
        X_HR_ref = squeeze(x_ref_all(t,:,:,i));
        X_LR = squeeze(X_LR_all(t,:,:,i));
        X_LR_interp = squeeze(X_LR_interp_all(t,:,:,i));      
        
        % Extract random patches: HR and LR features
        patches_LR = zeros(dim_l, num_patch);
        patches_LRinterp = zeros(dim_h, num_patch);
        patches_HR = zeros(dim_h, num_patch);

        % loop all the patches
        for ii = 1:num_patch
            row_h = xrow_h(ii); col_h = ycol_h(ii);
            row_l = xrow_l(ii); col_l = ycol_l(ii);

            % LR
            patch_LR = X_LR(row_l:row_l+patchsize_l-1,col_l:col_l+patchsize_l-1);
            patches_LR(:,ii) = patch_LR(:);
            
            % LRinterp
            patch_LRinterp = X_LR_interp(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_LRinterp(:,ii) = patch_LRinterp(:);              

            % HR
            patch_HR = X_HR_ref(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_HR(:,ii) = patch_HR(:);              
        end 
        patches_LR_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id) = patches_LR; 
        patches_LRinterp_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id) = patches_LRinterp; 
        patches_HR_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id)= patches_HR;
    end
end 

filename=strcat('/data/PhDworks/isotropic/dictionarylearning/downsampling/trainingpatches_patchHR_patchLR_patchLRinterp_spaceratio'...
    ,num2str(space_spacing,'%.1d'),'_timeratio',num2str(time_spacing,'%.1d'),'_patchsize',num2str(size_l,'%.2d'),'_numpatch'...
    ,num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
save(filename,'patches_LR_all','patches_LRinterp_all','patches_HR_all','-v7.3');
