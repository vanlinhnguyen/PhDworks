%LEARNODL_PATCHESHRLR_PATCHESLRFEATURES_PATCHSIZE06  learn coupled dictionaries
%between residual and feature patches (SR3)
%
%IN:
%   filename_patch - file name to load patches
%   params - parameters to learn the dictionary

%OUT:
%   filename_dict - file name to save learned dictionaries

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function filename_dict = learnODL_patchesHRLR_patchesLRfeatures_patchsize06(filename_patch,params)
fprintf('\n ********** LEARN COUPLE DICTS RESIDUAL-FEATURES ********** \n');
%% INITIAL PARAMS
space_spacing=4;
size_l=6; 
dim_h=(size_l*space_spacing)^2;
dim_l=4*dim_h;

%% Load patches
load(filename_patch, 'patches_HRLR_all','patches_LR_features_pca','V_pca','dim_l_pca');

% SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
start=tic();

patches_HRLR_all=patches_HRLR_all-repmat(mean(patches_HRLR_all),[dim_h 1]);
patches_HRLR_all=patches_HRLR_all ./ repmat(sqrt(sum(patches_HRLR_all.^2)),[dim_h 1]); 

patches_LR_features_pca=patches_LR_features_pca-repmat(mean(patches_LR_features_pca),[dim_l_pca 1]);
patches_LR_features_pca=patches_LR_features_pca ./ repmat(sqrt(sum(patches_LR_features_pca.^2)),[dim_l_pca 1]); 

fprintf(['Processing data in ',num2str(toc(start),'%.3f'),' seconds \n']);



%% ONLINE DICTIONARY LEARNING
fprintf(['Train dictionary for lambda=',num2str(params.lambda,'%.3f'),', K=',num2str(params.K,'%.4d'),'\n']);

[D_LR,~] = mexTrainDL(patches_LR_features_pca, params);
CoefMatrix=mexLasso(patches_LR_features_pca,D_LR,params);
D_HR = (patches_HRLR_all * CoefMatrix')/(full(CoefMatrix * CoefMatrix'));

%%
filename_dict=strcat('/data/PhDworks/isotropic/dictionarylearning/downsampling/DICTIONARY_patchesHRLR_patchesLRfeatures_K',...
    num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'_patchsize',num2str(size_l,'%.2d'),'.mat');
save(filename_dict,'D_HR','D_LR','V_pca','CoefMatrix');