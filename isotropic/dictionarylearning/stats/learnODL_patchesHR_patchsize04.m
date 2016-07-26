%LEARNODL_PATCHESHR_PATCHSIZE04 learn PCA basis from HR patches of size 06
%at LR
%
%IN:
%   patches_HR_all - input HR patches

%OUT:
%   D - ODL dictionary
%   A - transformation coefficients

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [D,A] = learnODL_patchesHR_patchsize04(patches_HR_all)

%% INITIAL PARAMS
space_spacing=4;
size_l=4; 
patchsize_h = size_l*space_spacing;
dim_h=patchsize_h^2;

%% Load patches
mu = mean(patches_HR_all,2);
patches_HR_all = bsxfun(@minus, patches_HR_all, mu); 

m_HR=mean(patches_HR_all,1);
for i=1:size(patches_HR_all,2)
    patches_HR_all(:,i) = patches_HR_all(:,i) - m_HR(i);
end
norm_HR=sqrt(sum(patches_HR_all.^2,1));
for i=1:dim_h
    patches_HR_all(i,:) = patches_HR_all(i,:)./norm_HR;
end

%% ONLINE DICTIONARY LEARNING
params.mode=2; % min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
params.K=2*dim_h; 
params.lambda=0.2;
params.lambda2=0; 
params.numThreads=-1; % number of threads (using all cores if -1)
params.iter=1000;  % max number of iterations.

fprintf(['Train dictionary for lambda=',num2str(params.lambda,'%.3f'),', K=',num2str(params.K,'%.4d'),'\n']);

[D,~] = mexTrainDL(patches_HR_all,params);
A=mexLasso(patches_HR_all,D,params);
