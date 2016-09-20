%SPARSITY_VS_ERROR estimate the NRMSE corresponding to each level of
%sparsity using PCA, wavelet, KSVD and ODL
%
%IN:
%   filename_dicts - name of file to load dictionaries (PCA, KSVD, ODL)
%   size_l - size of LR patches
%   num_snaps - number of fields to test

%OUT:
%   sparsity_(PCA,WL,KSVD,ODL) - sparsity of PCA, WL, KSVD and ODL (0 to 1)
%   NRMSE_(PCA,WL,KSVD,ODL) - equivalent NRMSE

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [sparsity_PCA, sparsity_WL, sparsity_KSVD, sparsity_ODL, NRMSE_PCA_all, ...
    NRMSE_WL_all, NRMSE_0DL_all, NRMSE_KSVD_all] = sparsity_vs_error(filename_dicts, size_l,num_snaps)

%% INITIAL PARAMS
Nh=96; % Number of grid-point in HR field
space_spacing=4;
patchsize_h = size_l*space_spacing;

overlap = size_l - 1; % at LR

% collect patches
border = 16; % to make the field of size 2^J
filename_ref='/data/PhDworks/isotropic/refdata_downsampled4.nc';
nc = netcdf(filename_ref,'r');

%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_h = 1: space_spacing*(size_l-overlap): Nh+1;
gridy_h = 1: space_spacing*(size_l-overlap): Nh+1;  

%% LOAD DICTIONARY 
load(filename_dicts);

lambdas = [1*1e-1,5*1e-2,1*1e-2,5*1e-3,1*1e-3,5*1e-4,1*1e-4,5*1e-5, 1e-5, 5*1e-6, 1e-6, 1e-7, 1e-8];
thres_WL = [exp(linspace(log(2),log(0.001),20)),0];
thres_pca = 0.05:0.05:1;

Ls = 8:16:patchsize_h^2;
NRMSE_0DL_all=zeros(length(lambdas),1);
NON_ZERO_COEFFS_0DL_all=zeros(length(lambdas),1);
NRMSE_KSVD_all=zeros(length(Ls),1);

NON_ZERO_COEFFS_WL_all=zeros(numel(thres_WL),1);
NRMSE_WL_all=zeros(numel(thres_WL),1);
NRMSE_PCA_all = zeros(numel(thres_pca),1);

params_ODL.mode=2; % LARS, min_{alpha} 0.5||x-Dalpha||_2^2 + lambda||alpha||_1 +0.5 lambda2||alpha||_2^2
params_ODL.lambda2=0; 

for t=1:num_snaps
    fprintf('\n Sparsity vs Error: snapshots  %.2d \n', t);

    %% HR field
    Xorg=nc{'velocity_x'}(t,:,:,4);
    
    Xorg_enlarged = [Xorg(Nh-border+1:Nh,Nh-border+1:Nh),Xorg(Nh-border+1:Nh,:),Xorg(Nh-border+1:Nh,1:border);...
        Xorg(:,Nh-border+1:Nh),Xorg,Xorg(:,1:border);...
        Xorg(1:border,Nh-border+1:Nh),Xorg(1:border,:),Xorg(1:border,1:border)];

    %% PCA
    patches_HR=extract_patches(Xorg_enlarged, {gridz_h,gridy_h}, patchsize_h);
    m_HR = mean(patches_HR,2);
    patches_HR = bsxfun(@minus, patches_HR, m_HR); 

    for i=1:length(thres_pca)
        K=round(thres_pca(i)*patchsize_h^2);
        Anew = patches_HR'*D_PCA(:,1:K);
        Anew = Anew./repmat(diag(D_PCA(:,1:K)'*D_PCA(:,1:K))',[size(patches_HR,2) 1]); % remember to normalize by norm of each basis
        patches_HR_reduced = Anew*D_PCA(:,1:K)';
        patches_HR_reduced = patches_HR_reduced + repmat(m_HR',[size(patches_HR,2) 1]);
        
        X_HR_rec_PCA=putback_patches(patches_HR_reduced', {gridz_h,gridy_h}, {Nh+2*border, Nh+2*border});
        X_HR_rec_PCA = X_HR_rec_PCA(border+1:end-border,border+1:end-border);

        NRMSE_PCA_all(i) = NRMSE_PCA_all(i) + 1/num_snaps*sqrt(sum((Xorg(:)-X_HR_rec_PCA(:)).^2))/sqrt(sum((Xorg(:)).^2));
    end
    
    %% ODL
    patches_HR=extract_patches(Xorg_enlarged, {gridz_h,gridy_h}, patchsize_h);
    m_HR=mean(patches_HR,1); norm_HR=sqrt(sum(patches_HR.^2,1));
    patches_HR=patches_HR - repmat(m_HR,patchsize_h^2,1);
    patches_HR=patches_HR./repmat(norm_HR,patchsize_h^2,1); 

    for i=1:length(lambdas)
        params_ODL.lambda=lambdas(i); 

        % solve sparse code and reconstruct HR patches
        A_ODL=mexLasso(patches_HR, D_ODL, params_ODL);

        patches_HR_rec = D_ODL*A_ODL;
        patches_HR_rec = patches_HR_rec.*repmat(norm_HR,patchsize_h^2,1);
        patches_HR_rec = patches_HR_rec + repmat(m_HR, patchsize_h^2, 1);

        X_HR_rec_DL=putback_patches(patches_HR_rec, {gridz_h,gridy_h}, {Nh+2*border, Nh+2*border});
        X_HR_rec_DL = X_HR_rec_DL(border+1:end-border,border+1:end-border);

        NRMSE_0DL_all(i) = NRMSE_0DL_all(i) + 1/num_snaps*sqrt(sum((X_HR_rec_DL(:)-Xorg(:)).^2))/sqrt(sum(Xorg(:).^2));
        NON_ZERO_COEFFS_0DL_all(i)= NON_ZERO_COEFFS_0DL_all(i) + 1/num_snaps*nnz(A_ODL)/size(A_ODL,2);
    end

   
    %% KSVD
    for i=1:length(Ls)
        params_KSVD.L = Ls(i); % set number of non-zero coefficients

        % solve sparse code and reconstruct HR patches
        A_KSVD=mexOMP(patches_HR,D_KSVD,params_KSVD);

        patches_HR_rec = D_KSVD*A_KSVD;
        patches_HR_rec = patches_HR_rec.*repmat(norm_HR,patchsize_h^2,1);
        patches_HR_rec = patches_HR_rec + repmat(m_HR, patchsize_h^2, 1);

        X_HR_rec_DL = putback_patches(patches_HR_rec, {gridz_h,gridy_h}, {Nh+2*border, Nh+2*border});
        X_HR_rec_DL = X_HR_rec_DL(border+1:end-border,border+1:end-border);

        NRMSE_KSVD_all(i) = NRMSE_KSVD_all(i) + 1/num_snaps*sqrt(sum((X_HR_rec_DL(:)-Xorg(:)).^2))/sqrt(sum(Xorg(:).^2));
    end


    %% Wavelet
    % Params
    L = 4 ; 
    J = 7; % 2D image of size 2^J
    jmax = 3; % number of octaves you want to filter at small scales

    % Wavelet decomposition
    wavelet = 'Daubechies';

    switch wavelet
        case 'Daubechies'
            qmf = MakeONFilter('Daubechies',4);
            wc = FWT2_PO(Xorg_enlarged,L,qmf);
            IWT = 'IWT2_PO';
        case 'Symmlet'
            qmf = MakeONFilter('Symmlet',4);
            plot(qmf)
            wc = FWT2_PO(Xorg_enlarged,L,qmf);
            IWT = 'IWT2_PO';
        case 'Battle'
            qmf = MakeONFilter('Battle',3);
            plot(qmf)
            wc = FWT2_PO(Xorg_enlarged,L,qmf);
            IWT = 'IWT2_PO';
        case 'Meyer'
            qmf=3;
            wc=FWT2_YM(Xorg_enlarged,2, L);
            IWT = 'IWT2_YM';
    end

    % extract only detail coefficients, don't touch 'approximation' coefficients
    W1 = wc(1:end,2^J-jmax+1:end);
    W2 = wc(2^(J-jmax)+1:end,1:2^(J-jmax));
    sigma = std([W1(:);W2(:)]);  % compute the variance of «small scale details » coeff.

    for i=1:numel(thres_WL)
        thres=thres_WL(i);
        wc_filt=zeros(size(wc));

        T = thres*sigma;    % choose threshold 

        % Inverse wavelet reconstruction
        wc_filt(abs(wc(:))>T) = wc(abs(wc(:))>T);
        wc_filt(1:2^(J-jmax),1:2^(J-jmax)) = wc(1:2^(J-jmax),1:2^(J-jmax)) ;

        eval(['X_HR_rec_WL = ',IWT,'(wc_filt,L,qmf);'])
        X_HR_rec_WL=X_HR_rec_WL(border+1:end-border,border+1:end-border);
        % Errors
        NRMSE_WL_all(i)=NRMSE_WL_all(i) + 1/num_snaps*sqrt(sum((Xorg(:)-X_HR_rec_WL(:)).^2))/sqrt(sum((Xorg(:)).^2));
        NON_ZERO_COEFFS_WL_all(i)= NON_ZERO_COEFFS_WL_all(i) + 1/num_snaps*nnz(wc_filt);
    end
end
close(nc);
NON_ZERO_COEFFS_0DL_all = round(NON_ZERO_COEFFS_0DL_all);
NON_ZERO_COEFFS_KSVD_all = Ls;

%% Convert to sparsity
sparsity_PCA = 1-thres_pca;
sparsity_WL = 1-NON_ZERO_COEFFS_WL_all./(Nh+2*border)^2;
sparsity_KSVD = 1-NON_ZERO_COEFFS_KSVD_all./(patchsize_h^2);
sparsity_ODL = 1-NON_ZERO_COEFFS_0DL_all./(patchsize_h^2);
