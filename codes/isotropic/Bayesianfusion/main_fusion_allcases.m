% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016
close all; clc;
addpath('./funcs/')

%% Bayesian fusion model
fusion_Bayesian(3,8) % case 1
fusion_Bayesian(4,8) % case 2
fusion_Bayesian(6,4) % case 3
fusion_Bayesian(6,6) % case 4
fusion_Bayesian(3,4) % case 5
fusion_Bayesian(4,6) % case 6
fusion_Bayesian(6,8) % case 7

%% Linear Gaussian fusion model
fusion_linearGauss(3,8) % case 1
fusion_linearGauss(4,8) % case 2
fusion_linearGauss(6,4) % case 3
fusion_linearGauss(6,6) % case 4
fusion_linearGauss(3,4) % case 5
fusion_linearGauss(4,6) % case 6
fusion_linearGauss(6,8) % case 7
