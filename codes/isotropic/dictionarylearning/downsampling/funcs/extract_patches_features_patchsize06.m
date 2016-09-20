%EXTRACT_PATCHES_FEATURES_PATCHSIZE06  extract features and residual (HR-interpolated LR) 
% patches (SR3). Patch size at LR is 6 x 6
%
%IN:
%   num_patch - number of patches extracted from each plane


%OUT:
%   filename - file name where data is stored

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function filename = extract_patches_features_patchsize06(num_patch)
fprintf('\n *************** EXTRACT PATCHES FEATURES  **************** \n');
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
dim_l=4*dim_h;

%% FILTERS 
O = zeros(1, space_spacing-1);
% O=0;
hf1 = [-1,O,1];
vf1 = hf1';
hf2 = [1,O,-2,O,1];
vf2 = hf2';


% boder to expand the field
border=patchsize_h/2;
left = border; right = border;
bottom=border; top = border;

nc2=netcdf(filename_interp,'r');
x_interp_all = nc2{'X_LRinterp_enlarged'}(1:Nt,:,:,LTHS_idt);
close(nc2);
%% EXTRACT PATCHES
num_planes=Nt*numel(LTHS_idt);

patches_LR_features_all = zeros(dim_l,num_patch*num_planes);
patches_HRLR_all = zeros(dim_h,num_patch*num_planes);

plane_id=0;
for t=1:Nt
    t
    for i=1:numel(LTHS_idt)
        plane_id=plane_id+1;
        
        % random positions of patches
        [xrow_l,ycol_l,xrow_h,ycol_h] = random_indices(Nl,patchsize_l,space_spacing, num_patch);

        % prepare the field: load, enlarge, downsample and interpolate
        X_HR_ref=enlarge_2D(squeeze(x_ref_all(t,:,:,i)),left, right, bottom, top);
        X_LR_interp = squeeze(x_interp_all(t,:,:,i));

        % extract features and truncate
        X_LR_Fea = zeros(size(X_HR_ref,1),size(X_HR_ref,2),4);
        X_LR_Fea(:, :, 1) = conv2(X_LR_interp, hf1, 'same');
        X_LR_Fea(:, :, 2) = conv2(X_LR_interp, vf1, 'same');
        X_LR_Fea(:, :, 3) = conv2(X_LR_interp,hf2,'same');
        X_LR_Fea(:, :, 4) = conv2(X_LR_interp,vf2,'same');
        
        % truncate
        X_LR_Fea = X_LR_Fea(bottom+1:bottom+Nh,left+1:left+Nh,:);
        X_HR_ref = X_HR_ref(bottom+1:bottom+Nh,left+1:left+Nh,:);
        X_LR_interp = X_LR_interp(bottom+1:bottom+Nh,left+1:left+Nh,:);
        X_HR_res = X_HR_ref - X_LR_interp;
        
        % Extract random patches: HR and LR features
        patches_LR_features = zeros(dim_l, num_patch);
        patches_HRLR = zeros(dim_h, num_patch);

        % loop all the patches
        for ii = 1:num_patch
            row_h = xrow_h(ii); col_h = ycol_h(ii);
            row_l = xrow_l(ii); col_l = ycol_l(ii);

            % LR (features)
            Lpatch_Fea1 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,1);
            Lpatch_Fea2 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,2);
            Lpatch_Fea3 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,3);
            Lpatch_Fea4 = X_LR_Fea(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1,4);
            patches_LR_features(:,ii) = [Lpatch_Fea1(:); Lpatch_Fea2(:); Lpatch_Fea3(:); Lpatch_Fea4(:)];

            % HR-LR (bicubic interp)
            Lpatch_HRLR = X_HR_res(row_h:row_h+patchsize_h-1,col_h:col_h+patchsize_h-1);
            patches_HRLR(:,ii) = Lpatch_HRLR(:);              
        end 
        patches_LR_features_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id) = patches_LR_features; 
        patches_HRLR_all(:,num_patch*(plane_id-1)+1:num_patch*plane_id)= patches_HRLR;
    end
end

%% Dimension reduction
start=tic();

patches_LR_features_all=patches_LR_features_all-repmat(mean(patches_LR_features_all),[dim_l 1]);
patches_LR_features_all=patches_LR_features_all ./ repmat(sqrt(sum(patches_LR_features_all.^2)),[dim_l 1]); 

C = patches_LR_features_all * patches_LR_features_all';
[V, D] = eig(C);
D = diag(D);
D = cumsum(D) / sum(D);

k = find(D >= 1e-3, 1); % ignore 0.1% energy
dim_l_pca = round(sqrt(numel(D)-k))^2; % just to make a number with integer sqrt, not important

V_pca = V(:, end-dim_l_pca+1:end); % ignore 0.23% energy
patches_LR_features_pca = V_pca' * patches_LR_features_all;

fprintf(['PCA in ',num2str(toc(start),'%.3f'),' seconds \n']);

clearvars patches_LR_features_all;
%% Save to file
filename=strcat('/data/PhDworks/isotropic/dictionarylearning/downsampling/trainingpatches_couplefeatures_spaceratio',...
    num2str(space_spacing,'%.1d'),'_timeratio',num2str(time_spacing,'%.1d'),'_patchsize',num2str(size_l,'%.2d'),...
    '_numpatch',num2str(Nt*numel(LTHS_idt)*num_patch,'%.6d'),'.mat');
save(filename,'patches_LR_features_pca','patches_HRLR_all','V_pca','dim_l_pca','-v7.3');
