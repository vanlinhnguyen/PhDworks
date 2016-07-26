%LEARNPCA_PATCHESHR_PATCHSIZE04 learn PCA basis from HR patches of size 06
%at LR
%
%IN: 
%   patches_HR_all - input HR patches

%OUT:
%   D - PCA dictionary
%   A - transformation coefficients

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [D,A,lambda] = learnPCA_patchesHR_patchsize04(patches_HR_all)

C = (patches_HR_all*patches_HR_all')./(size(patches_HR_all,2)-1); % cov(X)

[D, lambda] = eig(C);
[lambda, order] = sort(diag(lambda), 'descend');       % sort cols high to low
D = D(:,order);

A = D'*patches_HR_all;