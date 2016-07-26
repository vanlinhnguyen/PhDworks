clear all; close all; clc;
warning off; 

%% INITIAL PARAMS
Nh=96; % Number of grid-point in HR field
space_spacing=4;
time_spacing=6;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;

Nl = Nh/space_spacing;

size_l=4; 
patchsize_l = size_l; 
patchsize_h = patchsize_l*space_spacing;
dim_h=patchsize_h^2;

overlap=3; % at LR
num_patch=113664; % total number of patches

params.lambda=0.1;
params.K=2*dim_h;  

% collect patches
bodersize=16;
filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');

N=5;

%% GRID OF FIRST NODES IN LR AND HR PATCHES (enlarged version)
gridz_h = 1: space_spacing*(patchsize_l-overlap): Nh+1;
gridy_h = 1: space_spacing*(patchsize_l-overlap): Nh+1;  
% gridz_h = 1: Nh+1;
% gridy_h = 1: Nh+1;

%% PCA

PATCHES_FILENAME=strcat('/data/ISOTROPIC/dictionary_learning/restart/sspacing04/trainingpatches_coupleHRLR_spaceratio4_timeratio6_patchsize04_numpatch085248.mat');
vars = load(PATCHES_FILENAME, 'patches_HR_all');
patches_HR_all = vars.patches_HR_all;

mu = mean(patches_HR_all,2);
patches_HR_all = bsxfun(@minus, patches_HR_all, mu); 
C = (patches_HR_all*patches_HR_all')./(size(patches_HR_all,2)-1); %'# cov(X)

[D, lambda] = eig(C);
[lambda, order] = sort(diag(lambda), 'descend');       %# sort cols high to low
D = D(:,order);

A = D'*patches_HR_all;

%% LOAD DICTIONARY 
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.K=2*dim_h; 
params.lambda=0.1;
params.lambda2=0; 
params.numThreads=4; % number of threads
params.iter=1000;  % max number of iterations.

FILENAME_ODL=strcat('/data/ISOTROPIC/dictionary_learning/restart/sspacing04/DICTIONARY_ODL_patchHR_K',num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'.mat');
data_ODL = load(FILENAME_ODL,'D','CoefMatrix');
D_ODL = data_ODL.D;
CoefMatrix_ODL = data_ODL.CoefMatrix;

ksvd_conf.dictsize = 2*dim_h;
ksvd_conf.Tdata = 16; % maximal sparsity: TBD
FILENAME_KSVD=strcat('/data/ISOTROPIC/dictionary_learning/restart/sspacing04/DICTIONARY_KSVD_patchHR_K',num2str(ksvd_conf.dictsize,'%.4d'),'_Tdata',strrep(num2str(ksvd_conf.Tdata,'%.2f'),'.',''),'.mat');
data_KSVD= load(FILENAME_KSVD,'D','CoefMatrix');
D_KSVD = data_KSVD.D;
CoefMatrix_KSVD = data_KSVD.CoefMatrix;

lambdas=[1*1e-1,5*1e-2,1*1e-2,5*1e-3,1*1e-3,5*1e-4,1*1e-4,5*1e-5, 1e-5, 5*1e-6, 1e-6, 1e-7, 1e-8,1e-9,1e-10];
thres_all=[2,1,0.5, 0.2, 0.1, 0.05, 0.02, 0.01];
thres_pca=0.05:0.05:1;

Ls=[16,32,48, 64,128,192, 256];
NRMSE_0DL_all=zeros(length(lambdas),1);
NON_ZERO_COEFFS_0DL_all=zeros(length(lambdas),1);
NRMSE_KSVD_all=zeros(length(Ls),1);
NON_ZERO_COEFFS_KSVD_all = Ls;
NON_ZERO_COEFFS_WL_all=zeros(numel(thres_all),1);
NRMSE_WL_all=zeros(numel(thres_all),1);
NRMSE_PCA_all = zeros(numel(thres_pca),1);

for t=1:N
    fprintf('\n Snapshots  %.2d \n', t);

    %% HR field
    Xorg=nc{'velocity_x'}(t,:,:,4);
    
    Xorg_enlarged = [Xorg(Nh-bodersize+1:Nh,Nh-bodersize+1:Nh),Xorg(Nh-bodersize+1:Nh,:),Xorg(Nh-bodersize+1:Nh,1:bodersize);...
        Xorg(:,Nh-bodersize+1:Nh),Xorg,Xorg(:,1:bodersize);...
        Xorg(1:bodersize,Nh-bodersize+1:Nh),Xorg(1:bodersize,:),Xorg(1:bodersize,1:bodersize)];

    %% PCA
    patches_HR=extract_patches(Xorg_enlarged, {gridz_h,gridy_h}, patchsize_h);
    m_HR = mean(patches_HR,2);
    patches_HR = bsxfun(@minus, patches_HR, m_HR); 

    for i=1:length(thres_pca)
        K=round(thres_pca(i)*16^2);
        Anew = patches_HR'*D(:,1:K);
        Anew = Anew./repmat(diag(D(:,1:K)'*D(:,1:K))',[size(patches_HR,2) 1]); % remember to normalize by norm of each basis
        patches_HR_reduced = Anew*D(:,1:K)';
        patches_HR_reduced = patches_HR_reduced + repmat(m_HR',[size(patches_HR,2) 1]);
        
        X_HR_rec_PCA=putback_patches(patches_HR_reduced', {gridz_h,gridy_h}, {Nh+2*bodersize, Nh+2*bodersize});
        X_HR_rec_PCA = X_HR_rec_PCA(bodersize+1:end-bodersize,bodersize+1:end-bodersize);

        NRMSE_PCA_all(i) = NRMSE_PCA_all(i) + 1/N*sqrt(sum((Xorg(:)-X_HR_rec_PCA(:)).^2))/sqrt(sum((Xorg(:)).^2));
    end
    
    %% ODL
    patches_HR=extract_patches(Xorg_enlarged, {gridz_h,gridy_h}, patchsize_h);
    m_HR=mean(patches_HR,1); norm_HR=sqrt(sum(patches_HR.^2,1));
    patches_HR=patches_HR - repmat(m_HR,patchsize_h^2,1);
    patches_HR=patches_HR./repmat(norm_HR,patchsize_h^2,1); 

    for i=1:length(lambdas)
        params.lambda=lambdas(i); 

        % solve sparse code and reconstruct HR patches
        CoefMatrix=mexLasso(patches_HR, D_ODL, params);

        patches_HR_rec = D_ODL*CoefMatrix;
        patches_HR_rec = patches_HR_rec.*repmat(norm_HR,patchsize_h^2,1);
        patches_HR_rec = patches_HR_rec + repmat(m_HR, patchsize_h^2, 1);

        X_HR_rec_DL=putback_patches(patches_HR_rec, {gridz_h,gridy_h}, {Nh+2*bodersize, Nh+2*bodersize});
        X_HR_rec_DL = X_HR_rec_DL(bodersize+1:end-bodersize,bodersize+1:end-bodersize);

        NRMSE_0DL_all(i) = NRMSE_0DL_all(i) + 1/N*sqrt(sum((X_HR_rec_DL(:)-Xorg(:)).^2))/sqrt(sum(Xorg(:).^2));
        NON_ZERO_COEFFS_0DL_all(i)= NON_ZERO_COEFFS_0DL_all(i) + 1/N*nnz(CoefMatrix)/size(CoefMatrix,2);
    end

   
    %% KSVM
    for i=1:length(Ls)
        params.L=Ls(i); 

        % solve sparse code and reconstruct HR patches
        CoefMatrix=mexLasso(patches_HR,D_KSVD,params);

        patches_HR_rec = D_KSVD*CoefMatrix;
        patches_HR_rec = patches_HR_rec.*repmat(norm_HR,patchsize_h^2,1);
        patches_HR_rec = patches_HR_rec + repmat(m_HR, patchsize_h^2, 1);

        X_HR_rec_DL=putback_patches(patches_HR_rec, {gridz_h,gridy_h}, {Nh+2*bodersize, Nh+2*bodersize});
        X_HR_rec_DL = X_HR_rec_DL(bodersize+1:end-bodersize,bodersize+1:end-bodersize);

        NRMSE_KSVD_all(i) = NRMSE_KSVD_all(i) + 1/N*sqrt(sum((X_HR_rec_DL(:)-Xorg(:)).^2))/sqrt(sum(Xorg(:).^2));
    end


    %% Wavelet
    % Params
    L = 4 ; 
    J = 7; % 2D image of size 2^J
    jmax = 3; % number of octaves you want to filter at small scales

    % Wavelet decomposition
    wavelet = 'Meyer';

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

    for i=1:numel(thres_all)
        thres=thres_all(i);
        wc_filt=zeros(size(wc));

        T = thres*sigma;    % choose threshold 

        % Inverse wavelet reconstruction
        wc_filt(abs(wc(:))>T) = wc(abs(wc(:))>T);
        wc_filt(1:2^(J-jmax),1:2^(J-jmax)) = wc(1:2^(J-jmax),1:2^(J-jmax)) ;

        eval(['X_HR_rec_WL = ',IWT,'(wc_filt,L,qmf);'])
        X_HR_rec_WL=X_HR_rec_WL(bodersize+1:end-bodersize,bodersize+1:end-bodersize);
        % Errors
        NRMSE_WL_all(i)=NRMSE_WL_all(i) + 1/N*sqrt(sum((Xorg(:)-X_HR_rec_WL(:)).^2))/sqrt(sum((Xorg(:)).^2));
        NON_ZERO_COEFFS_WL_all(i)= NON_ZERO_COEFFS_WL_all(i) + 1/N*nnz(wc_filt);
    end
end
close(nc);
NON_ZERO_COEFFS_0DL_all = round(NON_ZERO_COEFFS_0DL_all);

% save('sparsity_vs_error_testing.mat','thres_pca', 'NRMSE_PCA_all', 'NON_ZERO_COEFFS_0DL_all','NRMSE_0DL_all',...
%     'NON_ZERO_COEFFS_KSVD_all', 'NRMSE_KSVD_all', 'NON_ZERO_COEFFS_WL_all', 'NRMSE_WL_all');
