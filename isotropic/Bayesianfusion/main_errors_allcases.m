% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016
close all; clc;
addpath('./funcs/')

err = @(xref,xrec) sqrt(sum((xref(:)-xrec(:)).^2))/sqrt(sum(xref(:).^2)); % NRMSE
% err = @(x1,x2) sum((x1(:)-x2(:)).^2)/sum(x2(:).^2); % NMSE

%% Bayesian fusion model 
[err_BF_mean1, err_BF_max1, ~, ~] = error_fusion_Bayesian(3,8,err); % case 1
fprintf(['BF- case 1- error (mean): ',num2str(err_BF_mean1,'%.2f'),', error max: ',num2str(err_BF_max1,'%.2f'),'\n']);

[err_BF_mean2,err_BF_max2] = error_fusion_Bayesian(4,8,err); % case 2  
fprintf(['BF- case 2- error (mean): ',num2str(err_BF_mean2,'%.2f'),', error max:',num2str(err_BF_max2,'%.2f'),'\n']);

[err_BF_mean3,err_BF_max3] = error_fusion_Bayesian(6,4,err); % case 3
fprintf(['BF- case 3- error (mean): ',num2str(err_BF_mean3,'%.2f'),', error max:',num2str(err_BF_max3,'%.2f'),'\n']);

[err_BF_mean4,err_BF_max4] = error_fusion_Bayesian(6,6,err); % case 4
fprintf(['BF- case 4- error (mean): ',num2str(err_BF_mean4,'%.2f'),', error max:',num2str(err_BF_max4,'%.2f'),'\n']);

[err_BF_mean5,err_BF_max5,err_BF_LF5, err_BF_HF5] = error_fusion_Bayesian(3,4,err); % case 5
fprintf(['BF- case 5- error (mean): ',num2str(err_BF_mean5,'%.2f'),', error max: ',num2str(err_BF_max5,'%.2f'),', large scales: ',num2str(err_BF_LF5,'%.2f'),', small scales: ',num2str(err_BF_HF5,'%.2f'),'\n']);

[err_BF_mean6,err_BF_max6,err_BF_LF6, err_BF_HF6] = error_fusion_Bayesian(4,6,err); % case 6
fprintf(['BF- case 6- error (mean): ',num2str(err_BF_mean6,'%.2f'),', error max: ',num2str(err_BF_max6,'%.2f'),', large scales: ',num2str(err_BF_LF6,'%.2f'),', small scales: ',num2str(err_BF_HF6,'%.2f'),'\n']);

[err_BF_mean7,err_BF_max7,err_BF_LF7, err_BF_HF7] = error_fusion_Bayesian(6,8,err); % case 7
fprintf(['BF- case 7- error (mean): ',num2str(err_BF_mean7,'%.2f'),', error max: ',num2str(err_BF_max7,'%.2f'),', large scales: ',num2str(err_BF_LF7,'%.2f'),', small scales: ',num2str(err_BF_HF7,'%.2f'),'\n']);

%% Linear Gaussian fusion model
[err_LG_mean1,err_LG_max1] = error_fusion_LG(3,8,err); % case 1
fprintf(['LG- case 1- error (mean): ',num2str(err_LG_mean1,'%.2f'),', error max: ',num2str(err_LG_max1,'%.2f'),'\n']);

[err_LG_mean2,err_LG_max2] = error_fusion_LG(4,8,err); % case 2
fprintf(['LG- case 2- error (mean): ',num2str(err_LG_mean2,'%.2f'),', error max: ',num2str(err_LG_max2,'%.2f'),'\n']);

[err_LG_mean3,err_LG_max3] = error_fusion_LG(6,4,err); % case 3
fprintf(['LG- case 3- error (mean): ',num2str(err_LG_mean3,'%.2f'),', error max: ',num2str(err_LG_max3,'%.2f'),'\n']);

[err_LG_mean4,err_LG_max4] = error_fusion_LG(6,6,err); % case 4
fprintf(['LG- case 4- error (mean): ',num2str(err_LG_mean4,'%.2f'),', error max: ',num2str(err_LG_max4,'%.2f'),'\n']);

[err_LG_mean5,err_LG_max5,err_LG_LF5, err_LG_HF5] = error_fusion_LG(3,4,err); % case 5
fprintf(['BF- case 5- error (mean): ',num2str(err_LG_mean5,'%.2f'),', error max: ',num2str(err_LG_max5,'%.2f'),', large scales: ',num2str(err_LG_LF5,'%.2f'),', small scales: ',num2str(err_LG_HF5,'%.2f'),'\n']);

[err_LG_mean6,err_LG_max6,err_LG_LF6, err_LG_HF6] = error_fusion_LG(4,6,err); % case 6
fprintf(['BF- case 6- error (mean): ',num2str(err_LG_mean6,'%.2f'),', error max: ',num2str(err_LG_max6,'%.2f'),', large scales: ',num2str(err_LG_LF6,'%.2f'),', small scales: ',num2str(err_LG_HF6,'%.2f'),'\n']);

[err_LG_mean7,err_LG_max7,err_LG_LF7, err_LG_HF7] = error_fusion_LG(6,8,err); % case 7
fprintf(['BF- case 7- error (mean): ',num2str(err_LG_mean7,'%.2f'),', error max: ',num2str(err_LG_max7,'%.2f'),', large scales: ',num2str(err_LG_LF7,'%.2f'),', small scales: ',num2str(err_LG_HF7,'%.2f'),'\n']);


%% Interpolation
% Space
[err_interp_space_mean1,err_interp_space_max1, err_interp_space_LF1, err_interp_space_HF1] = error_interp_space(3,err);
fprintf(['Interpolation (space) of ratio 3x3 - error (mean): ',num2str(err_interp_space_mean1,'%.2f'),', error max: ', ...
    num2str(err_interp_space_max1,'%.2f'),', large scales: ',num2str(err_interp_space_LF1,'%.2f'),', small scales: ',num2str(err_interp_space_HF1,'%.2f'),'\n']);

[err_interp_space_mean2, err_interp_space_max2, err_interp_space_LF2, err_interp_space_HF2] = error_interp_space(4,err);
fprintf(['Interpolation (space) of ratio 4x4 - error (mean): ',num2str(err_interp_space_mean2,'%.2f'),', error max: ', ...
    num2str(err_interp_space_max2,'%.2f'),', large scales: ',num2str(err_interp_space_LF2,'%.2f'),', small scales: ',num2str(err_interp_space_HF2,'%.2f'),'\n']);

[err_interp_space_mean3, err_interp_space_max3, err_interp_space_LF3, err_interp_space_HF3] = error_interp_space(6,err);
fprintf(['Interpolation (space) of ratio 6x6 - error (mean): ',num2str(err_interp_space_mean3,'%.2f'),', error max: ', ...
    num2str(err_interp_space_max3,'%.2f'),', large scales: ',num2str(err_interp_space_LF3,'%.2f'),', small scales: ',num2str(err_interp_space_HF3,'%.2f'),'\n']);

% Time
[err_interp_time_mean1, err_interp_time_max1, err_time_space_LF1, err_time_space_HF1] = error_interp_time(4,err);
fprintf(['Interpolation (time) of ratio 4 - error (mean): ',num2str(err_interp_time_mean1,'%.2f'),', error max: ', ...
    num2str(err_interp_time_max1,'%.2f'),', large scales: ',num2str(err_time_space_LF1,'%.2f'),', small scales: ',num2str(err_time_space_HF1,'%.2f'),'\n']);

[err_interp_time_mean2,err_interp_time_max2] = error_interp_time(6,err);
fprintf(['Interpolation (time) of ratio 6: ',num2str(err_interp_time_mean2,'%.2f'),', error max: ',num2str(err_interp_time_max2,'%.2f'),'\n']);

[err_interp_time_mean3,err_interp_time_max3] = error_interp_time(8,err);
fprintf(['Interpolation (time) of ratio 8: ',num2str(err_interp_time_mean3,'%.2f'),', error max: ',num2str(err_interp_time_max3,'%.2f'),'\n']);
