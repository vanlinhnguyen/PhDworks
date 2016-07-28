%LEARNKSVD learn dictionary using KSVD
%
%IN:
%   patches_all - input HR patches

%OUT:
%   D - ODL dictionary
%   A - transformation coefficients

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [D,A] = learnKSVD(patches_all, params)

%% INITIAL PARAMS
dim_h = size(patches_all,1);
num_patch = size(patches_all,2);

%% Load patches
mu = mean(patches_all,2);
patches_all = bsxfun(@minus, patches_all, mu); 

m_HR=mean(patches_all,1);
for i=1:num_patch
    patches_all(:,i) = patches_all(:,i) - m_HR(i);
end
norm_HR=sqrt(sum(patches_all.^2,1));
for i=1:dim_h
    patches_all(i,:) = patches_all(i,:)./norm_HR;
end

%% KSVD
fprintf(['Train dictionary for T=',num2str(params.Tdata,'%.3f'),', K=',num2str(params.dictsize,'%.4d'),'\n']);
params.data = patches_all;
[D,A,~] = ksvd(params,'');