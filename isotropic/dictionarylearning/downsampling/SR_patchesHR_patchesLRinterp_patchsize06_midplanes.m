%SR_PATCHESHRLR_PATCHESLRINTERP_PATCHSIZE06_MIDPLANES  using SR2 to
%super-resolve all midplanes (far from LTHS planes)
%
%IN:
%   filename_dict - file name to load learned dictionaries
%   params_train - parameters to learn the dictionary
%   params_rec - parameters to super-resolve
%OUT:
%   filename_SR - file name to save reconstructed planes

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function filename_SR = SR_patchesHR_patchesLRinterp_patchsize06_midplanes(filename_dict, params_train, params_rec)
fprintf('\n ************ SR - COUPLE DICTS HR-LRINTERP ************** \n');

%% INITIAL PARAMS
space_spacing=4;
time_spacing=6;
Nt = 37;
Nh = 96;
Nl=Nh/space_spacing;

size_l=6; 
patchsize_h = size_l*space_spacing;
dim_h=patchsize_h^2;
dim_l=dim_h;

overlap=size_l-1; % at LR

% boder to expand the field
border = patchsize_h/2;
left = border; right = border;
bottom=border; top = border;

%% LOAD DICTIONARY 
load(filename_dict);
D_HR = D_HR./repmat(sqrt(sum(D_HR.^2, 1)), dim_h, 1);
D_HR=sqrt(dim_h/dim_l)*D_HR; % scale the HR dictionary (very important!!!)
D_LR = D_LR./repmat(sqrt(sum(D_LR.^2, 1)), dim_l, 1);

%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_l = 1: size_l-overlap : Nl+1;
gridy_l = 1: size_l-overlap : Nl+1;
gridz_h = 1: space_spacing*(size_l-overlap) : Nh+1;
gridy_h = 1: space_spacing*(size_l-overlap) : Nh+1;  
grids={gridz_l,gridy_l,gridz_h,gridy_h};
sizes={Nl+size_l,Nl+size_l,Nh+patchsize_h,Nh+patchsize_h};

%% Loop all blocks
midplane_ids = 4:time_spacing:Nh;

filename_SR=strcat('/data/PhDworks/isotropic/dictionarylearning/downsampling/SR_HR_LRinterp_K',...
    num2str(params_train.K,'%.4d'),'_lambda',strrep(num2str(params_train.lambda,'%.2f'),'.',''),'_allmidplanes.nc');

fprintf('\n ********** START SUPERRESOLUTION ************ \n');
nc = netcdf(filename_SR,'clobber');
nc('Nt')=0;
nc('Nx')=numel(midplane_ids);
nc('Ny')=Nh;
nc('Nz')=Nh;
nc{'Zhat_all'}=ncfloat('Nt','Nz','Ny','Nx');

filename_ref = '/data/PhDworks/isotropic/refdata_downsampled4.nc';
nc1 = netcdf(filename_ref,'r');

filename_interp = '/data/PhDworks/isotropic/IsY_sspacing04.nc';
nc2=netcdf(filename_interp,'r');

for t=1:Nt
    fprintf('Superresolution: %.2d-th block.\n', t);
    Zhat_oneblock=zeros(Nh,Nh,numel(midplane_ids));
    for i=1:numel(midplane_ids)
        X_HR_org=nc1{'velocity_x'}(t,1:Nh,1:Nh,midplane_ids(i));
        X_LR_interp_enlarged = nc2{'X_LRinterp_enlarged'}(t,:,:,midplane_ids(i));

        % SR
        start=tic();

        X_HR_rec_enlarged = superresol(X_LR_interp_enlarged, D_HR, D_LR, params_rec,'SR2', grids, sizes);

        % REMOVE BOLDER (translation is not complete near bolder)
        X_HR_rec=X_HR_rec_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);
        X_LR_interp=X_LR_interp_enlarged(bottom+1:bottom+Nh,left+1:left+Nh);

        NRMSE1=sqrt(sum((X_LR_interp(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
        NRMSE2=sqrt(sum((X_HR_rec(:)-X_HR_org(:)).^2))/sqrt(sum(X_HR_org(:).^2));
       
        Zhat_oneblock(:,:,i)=X_HR_rec;
        fprintf(['SR snapshot t=',num2str(midplane_ids(i),'%.3d'),' in ',num2str(toc(start)/60,'%.3f'), ...
            ' minute(s), improves NRMSE by ', num2str((1-NRMSE2/NRMSE1)*100,'%.2f'),' percents\n']);
    end
    nc{'Zhat_all'}(t,:,:,:) = Zhat_oneblock;
end
close(nc); close(nc1); close(nc2);
