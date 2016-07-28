%LEARNODL_PATCHESHR_PATCHESLR_PATCHSIZE06  learn coupled dictionaries
%between HR and LR patches (SR1)
%
%IN:
%   filename_patch - file name to load patches
%   params - parameters to learn the dictionary

%OUT:
%   filename_dict - file name to save learned dictionaries

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function filename_dict = learnODL_patchesHR_patchesLR_patchsize06(filename_patch,params)
fprintf('\n *************** LEARN COUPLE DICTS HR-LR  **************** \n');

%% INITIAL PARAMS
space_spacing=4;
size_l=6; 
dim_h=(size_l*space_spacing)^2;
dim_l=size_l^2;

%% Load patches
load(filename_patch, 'patches_LR_all','patches_HR_all');

% SUBSTRACT MEAN AND PRE-NORMALIZE ALL PATCHES
start=tic();
patches_LR_all=patches_LR_all-repmat(mean(patches_LR_all),[dim_l 1]);
patches_LR_all=patches_LR_all ./ repmat(sqrt(sum(patches_LR_all.^2)),[dim_l 1]); 

patches_HR_all=patches_HR_all-repmat(mean(patches_HR_all),[dim_h 1]);
patches_HR_all=patches_HR_all ./ repmat(sqrt(sum(patches_HR_all.^2)),[dim_h 1]); % the normalization of dictionaries is from here: scales of sum(patches_HR_all.^2) and sum(patches_LR_all.^2 are different

% JOINT PATCHES
patches_all = [(1/sqrt(dim_h))*patches_HR_all; (1/sqrt(dim_l))*patches_LR_all];
patch_norm = sqrt(sum(patches_all.^2, 1));
patches_all = patches_all(:, patch_norm > 1e-5); 
patches_all = patches_all - repmat(mean(patches_all,1),dim_h+dim_l,1);
patches_all = patches_all./repmat(patch_norm, dim_h+dim_l, 1);

fprintf(['Processing data in ',num2str(toc(start),'%.3f'),' seconds \n']);

%% ONLINE DICTIONARY LEARNING

fprintf(['Train dictionary for lambda=',num2str(params.lambda,'%.3f'),', K=',num2str(params.K,'%.4d'),'\n']);

% JOINT LEARNING
[D,~] = mexTrainDL(patches_all,params);
CoefMatrix=mexLasso(patches_all,D,params);
D_HR = D(1:dim_h,:);
D_LR = D(dim_h+1:end,:);

%%
filename_dict=strcat('/data/PhDworks/isotropic/dictionarylearning/downsampling/DICTIONARY_patchesHR_patchesLR_joint_K',...
    num2str(params.K,'%.4d'),'_lambda',strrep(num2str(params.lambda,'%.2f'),'.',''),'_patchsize',num2str(size_l,'%.2d'),'.mat');
save(filename_dict,'D_HR','D_LR','CoefMatrix');