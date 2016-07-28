%LEARNPCA learn PCA basis
%
%IN: 
%   patches_all - input HR patches

%OUT:
%   D - PCA dictionary
%   A - transformation coefficients

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [D,A,lambda] = learnPCA(patches_all)

C = (patches_all*patches_all')./(size(patches_all,2)-1); % cov(X)

[D, lambda] = eig(C);
[lambda, order] = sort(diag(lambda), 'descend');       % sort cols high to low
D = D(:,order);

A = D'*patches_all;