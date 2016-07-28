% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

close all; clc;
addpath('./funcs/')

%% INITIAL PARAMS
space_spacing=4;
time_spacing=6;
Nt = 37;
Nh = 96;
Nl=Nh/space_spacing;
size_l=6; 

%% Extract patches
fprintf('\n ***********************************************************\n');
fprintf('\n ******************* EXTRACT PATCHES ********************** \n');
fprintf('\n ***********************************************************\n');

num_patch = 12*12;
filename_patch_HRLR = extract_patches_LRHR_patchsize06(num_patch);
filename_patch_features = extract_patches_features_patchsize06(num_patch);

%% Learn coupled dicts
fprintf('\n ***********************************************************\n');
fprintf('\n ************ TRAIN COUPLED DICTIONARIES ****************** \n');
fprintf('\n ***********************************************************\n');

% Couple HR-LR
dim_h=(size_l*space_spacing)^2;
dim_l=size_l^2;
params_train_HRLR.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params_train_HRLR.lambda=0.1;
params_train_HRLR.lambda2=0;
params_train_HRLR.K=2*(dim_h+dim_l);
params_train_HRLR.numThreads=-1; % number of threads
params_train_HRLR.iter = 1000;  % max number of iterations.

filename_dict_HRLR = learnODL_patchesHR_patchesLR_patchsize06(filename_patch_HRLR,params_train_HRLR);

% Couple HR-LRinterp
dim_h = (size_l*space_spacing)^2;
dim_l = dim_h;
params_train_HRLRinterp.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params_train_HRLRinterp.lambda=0.1;
params_train_HRLRinterp.lambda2=0;
params_train_HRLRinterp.K=2*dim_h;
params_train_HRLRinterp.numThreads=-1; % number of threads
params_train_HRLRinterp.iter = 1000;  % max number of iterations.

filename_dict_HRLRinterp = learnODL_patchesHR_patchesLRinterp_patchsize06(filename_patch_HRLR,params_train_HRLR);


% Couple residual-features
dim_h = (size_l*space_spacing)^2;
load(filename_patch_features, 'dim_l_pca');

params_train_features.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params_train_features.lambda=0.1;
params_train_features.lambda2=0;
params_train_features.K=2*(dim_h+dim_l_pca);
params_train_features.numThreads=-1; % number of threads
params_train_features.iter = 1000;  % max number of iterations.

filename_dict_features = learnODL_patchesHRLR_patchesLRfeatures_patchsize06(filename_patch_features,params_train_features);

%% Super-resolution
fprintf('\n ***********************************************************\n');
fprintf('\n ***************** SUPER-RESOLUTION *********************** \n');
fprintf('\n ***********************************************************\n');

% filename_dict_HRLR = '/data/PhDworks/isotropic/dictionarylearning/downsampling/DICTIONARY_patchesHR_patchesLR_joint_K1224_lambda020_patchsize06.mat';
% filename_dict_HRLRinterp = '/data/PhDworks/isotropic/dictionarylearning/downsampling/DICTIONARY_patchesHR_patchesLRinterp_joint_K1152_lambda020_patchsize06.mat';
% filename_dict_features = '/data/PhDworks/isotropic/dictionarylearning/downsampling/DICTIONARY_patchesHRLR_patchesLRfeatures_K1314_lambda020_patchsize06.mat';

params_SR.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params_SR.lambda=1e-6; 
params_SR.lambda2=0;
params_SR.numThreads=-1; % number of threads

SR_patchesHR_patchesLR_patchsize06_midplanes(filename_dict_HRLR, params_train_HRLR, params_SR)
SR_patchesHR_patchesLRinterp_patchsize06_midplanes(filename_dict_HRLRinterp, params_train_HRLRinterp, params_SR)
SR_patchesHRLR_patchesLRfeatures_patchsize06_midplanes(filename_dict_features, params_train_features, params_SR)

