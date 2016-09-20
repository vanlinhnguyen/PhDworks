clear all; close all; clc;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';

%%  TIME: 04; SPACE: 03
space_spacing=03; % subsampling ration in space
time_spacing=04; % subsampling ration in time (from 40Hz to 4Hz)
LTHS_idt=1:time_spacing:96;

fprintf('compute errors TIME: 04  SPACE: 03\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_LG=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_NLM_greedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_NLM_nongreedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_RR=strcat('/data/ISOTROPIC/regression/RR_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 
filename_KRR_rbf=strcat('/data/ISOTROPIC/regression/KRR_rbf_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_BF,'r');
nc5=netcdf(filename_LG,'r');
nc6=netcdf(filename_NLM_greedy,'r');
nc7=netcdf(filename_NLM_nongreedy,'r');
nc8=netcdf(filename_RR,'r');
nc9=netcdf(filename_KRR_rbf,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,1:LTHS_idt(end));
x_interp_spatial=nc2{'Uinterp'}(:,:,:,1:LTHS_idt(end));
x_interp_temporal=nc3{'Uinterp'}(:,:,:,1:LTHS_idt(end));
x_fusion_BF=nc4{'Zhat_all'}(:,:,:,1:LTHS_idt(end));
x_fusion_LG=nc5{'Zhat_all'}(:,:,:,1:LTHS_idt(end));
x_NLM_greedy=permute(nc6{'x_HR_NLM_smallscales_all'}(:,1:LTHS_idt(end),:,:),[1,3,4,2])+x_interp_spatial;
x_NLM_nongreedy=permute(nc7{'x_rec_all'}(:,1:LTHS_idt(end),:,:),[1,3,4,2])+x_interp_spatial;
x_RR=nc8{'Urec'}(:,:,:,1:LTHS_idt(end));
x_KRR_rbf=nc9{'Urec'}(:,:,:,1:LTHS_idt(end));

x_RR(:,:,:,1:time_spacing:LTHS_idt(end))=x_nonoise(:,:,:,1:time_spacing:LTHS_idt(end));
x_KRR_rbf(:,:,:,1:time_spacing:LTHS_idt(end))=x_nonoise(:,:,:,1:time_spacing:LTHS_idt(end));

close(nc1); close(nc2); close(nc3); close(nc4);
close(nc5);  close(nc6); close(nc7);  close(nc8);  close(nc9);

err_interp_spatial_time4_space3=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time4_space3=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_BF_time4_space3=sqrt(sum((x_nonoise(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_LG_time4_space3=sqrt(sum((x_nonoise(:)-x_fusion_LG(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_NLM_greedy_time4_space3=sqrt(sum((x_nonoise(:)-x_NLM_greedy(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_NLM_nongreedy_time4_space3=sqrt(sum((x_nonoise(:)-x_NLM_nongreedy(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_RR_time4_space3=sqrt(sum((x_nonoise(:)-x_RR(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_KRR_rbf_time4_space3=sqrt(sum((x_nonoise(:)-x_KRR_rbf(:)).^2))/sqrt(sum((x_nonoise(:)).^2));

%%  TIME: 06; SPACE: 04
space_spacing=4; % subsampling ration in space
time_spacing=6; % subsampling ration in time (from 40Hz to 4Hz)
LTHS_idt=1:time_spacing:96;

fprintf('compute errors TIME: 06  SPACE: 04\n');
filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_LG=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_NLM_greedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_NLM_nongreedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_RR=strcat('/data/ISOTROPIC/regression/RR_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 
filename_KRR_rbf=strcat('/data/ISOTROPIC/regression/KRR_rbf_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_BF,'r');
nc5=netcdf(filename_LG,'r');
nc6=netcdf(filename_NLM_greedy,'r');
nc7=netcdf(filename_NLM_nongreedy,'r');
nc8=netcdf(filename_RR,'r');
nc9=netcdf(filename_KRR_rbf,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,1:LTHS_idt(end));
x_interp_spatial=nc2{'Uinterp'}(:,:,:,1:LTHS_idt(end));
x_interp_temporal=nc3{'Uinterp'}(:,:,:,1:LTHS_idt(end));
x_fusion_BF=nc4{'Zhat_all'}(:,:,:,1:LTHS_idt(end));
x_fusion_LG=nc5{'Zhat_all'}(:,:,:,1:LTHS_idt(end));
x_NLM_greedy=permute(nc6{'x_HR_NLM_smallscales_all'}(:,1:LTHS_idt(end),:,:),[1,3,4,2])+x_interp_spatial;
x_NLM_nongreedy=permute(nc7{'x_rec_all'}(:,1:LTHS_idt(end),:,:),[1,3,4,2])+x_interp_spatial;
x_RR=nc8{'Urec'}(:,:,:,1:LTHS_idt(end));
x_KRR_rbf=nc9{'Urec'}(:,:,:,1:LTHS_idt(end));

x_RR(:,:,:,1:time_spacing:LTHS_idt(end))=x_nonoise(:,:,:,1:time_spacing:LTHS_idt(end));
x_KRR_rbf(:,:,:,1:time_spacing:LTHS_idt(end))=x_nonoise(:,:,:,1:time_spacing:LTHS_idt(end));

close(nc1); close(nc2); close(nc3); close(nc4);
close(nc5);  close(nc6); close(nc7);  close(nc8);  close(nc9);

err_interp_spatial_time6_space4=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time6_space4=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_BF_time6_space4=sqrt(sum((x_nonoise(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_LG_time6_space4=sqrt(sum((x_nonoise(:)-x_fusion_LG(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_NLM_greedy_time6_space4=sqrt(sum((x_nonoise(:)-x_NLM_greedy(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_NLM_nongreedy_time6_space4=sqrt(sum((x_nonoise(:)-x_NLM_nongreedy(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_RR_time6_space4=sqrt(sum((x_nonoise(:)-x_RR(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_KRR_rbf_time6_space4=sqrt(sum((x_nonoise(:)-x_KRR_rbf(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
 

%%  TIME: 08; SPACE: 06
space_spacing=6; % subsampling ration in space
time_spacing=8; % subsampling ration in time (from 40Hz to 4Hz)
LTHS_idt=1:time_spacing:96;

fprintf('compute errors TIME: 08  SPACE: 06\n');

filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_LG=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
filename_NLM_greedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_NLM_nongreedy=strcat('/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'_sim12_acc12_neighbor5_tau0100.nc');
filename_RR=strcat('/data/ISOTROPIC/regression/RR_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 
filename_KRR_rbf=strcat('/data/ISOTROPIC/regression/KRR_rbf_sspacing',num2str(space_spacing,'%.1d'),'_tspacing',num2str(time_spacing,'%.1d'),'.nc'); 

nc1=netcdf(filename_ref,'r');
nc2 = netcdf(filename_interp_space,'r');
nc3=netcdf(filename_interp_time,'r');
nc4=netcdf(filename_BF,'r');
nc5=netcdf(filename_LG,'r');
nc6=netcdf(filename_NLM_greedy,'r');
nc7=netcdf(filename_NLM_nongreedy,'r');
nc8=netcdf(filename_RR,'r');
nc9=netcdf(filename_KRR_rbf,'r');

x_nonoise=nc1{'velocity_x'}(:,:,:,1:LTHS_idt(end));
x_interp_spatial=nc2{'Uinterp'}(:,:,:,1:LTHS_idt(end));
x_interp_temporal=nc3{'Uinterp'}(:,:,:,1:LTHS_idt(end));
x_fusion_BF=nc4{'Zhat_all'}(:,:,:,1:LTHS_idt(end));
x_fusion_LG=nc5{'Zhat_all'}(:,:,:,1:LTHS_idt(end));
x_NLM_greedy=permute(nc6{'x_HR_NLM_smallscales_all'}(:,1:LTHS_idt(end),:,:),[1,3,4,2])+x_interp_spatial;
x_NLM_nongreedy=permute(nc7{'x_rec_all'}(:,1:LTHS_idt(end),:,:),[1,3,4,2])+x_interp_spatial;
x_RR=nc8{'Urec'}(:,:,:,1:LTHS_idt(end));
x_KRR_rbf=nc9{'Urec'}(:,:,:,1:LTHS_idt(end));

x_RR(:,:,:,1:time_spacing:LTHS_idt(end))=x_nonoise(:,:,:,1:time_spacing:LTHS_idt(end));
x_KRR_rbf(:,:,:,1:time_spacing:LTHS_idt(end))=x_nonoise(:,:,:,1:time_spacing:LTHS_idt(end));

close(nc1); close(nc2); close(nc3); close(nc4);
close(nc5);  close(nc6); close(nc7);  close(nc8);  close(nc9);

err_interp_spatial_time8_space6=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_interp_temporal_time8_space6=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_BF_time8_space6=sqrt(sum((x_nonoise(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_fusion_LG_time8_space6=sqrt(sum((x_nonoise(:)-x_fusion_LG(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_NLM_greedy_time8_space6=sqrt(sum((x_nonoise(:)-x_NLM_greedy(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_NLM_nongreedy_time8_space6=sqrt(sum((x_nonoise(:)-x_NLM_nongreedy(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_RR_time8_space6=sqrt(sum((x_nonoise(:)-x_RR(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
err_KRR_rbf_time8_space6=sqrt(sum((x_nonoise(:)-x_KRR_rbf(:)).^2))/sqrt(sum((x_nonoise(:)).^2));

%%  TIME: 08; SPACE: 03 
% space_spacing=03; % subsampling ration in space
% time_spacing=08; % subsampling ration in time (from 40Hz to 4Hz)
% 
% fprintf('compute errors TIME: 08  SPACE: 03\n');
% 
% filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
% filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
% filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
% 
% nc1=netcdf(filename_ref,'r');
% nc2 = netcdf(filename_interp_space,'r');
% nc3=netcdf(filename_interp_time,'r');
% nc4=netcdf(filename_BF,'r');
% 
% x_nonoise=nc1{'velocity_x'}(:,:,:,:);
% x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
% x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
% x_fusion_BF=nc4{'Zhat_all'}(:,:,:,:);
% 
% close(nc1); close(nc2); close(nc3); close(nc4);
% 
% err_interp_spatial_time8_space3=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_interp_temporal_time8_space3=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_fusion_time8_space3=sqrt(sum((x_nonoise(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 08; SPACE: 04
% space_spacing=04; % subsampling ration in space
% time_spacing=08; % subsampling ration in time (from 40Hz to 4Hz)
% 
% fprintf('compute errors TIME: 08  SPACE: 04\n');
% 
% filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
% filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
% filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
% 
% nc1=netcdf(filename_ref,'r');
% nc2 = netcdf(filename_interp_space,'r');
% nc3=netcdf(filename_interp_time,'r');
% nc4=netcdf(filename_BF,'r');
% 
% x_nonoise=nc1{'velocity_x'}(:,:,:,:);
% x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
% x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
% x_fusion_BF=nc4{'Zhat_all'}(:,:,:,:);
% 
% close(nc1); close(nc2); close(nc3); close(nc4);
% 
% err_interp_spatial_time8_space4=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_interp_temporal_time8_space4=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_fusion_time8_space4=sqrt(sum((x_nonoise(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 04; SPACE: 06
% space_spacing=6; % subsampling ration in space
% time_spacing=4; % subsampling ration in time (from 40Hz to 4Hz)
% 
% fprintf('compute errors TIME: 04  SPACE: 06\n');
% 
% filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
% filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
% filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
% 
% nc1=netcdf(filename_ref,'r');
% nc2 = netcdf(filename_interp_space,'r');
% nc3=netcdf(filename_interp_time,'r');
% nc4=netcdf(filename_BF,'r');
% 
% x_nonoise=nc1{'velocity_x'}(:,:,:,:);
% x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
% x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
% x_fusion_BF=nc4{'Zhat_all'}(:,:,:,:);
% 
% close(nc1); close(nc2); close(nc3); close(nc4);
% 
% err_interp_spatial_time4_space6=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_interp_temporal_time4_space6=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_fusion_time4_space6=sqrt(sum((x_nonoise(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_nonoise(:)).^2));


%%  TIME: 06; SPACE: 06
% space_spacing=6; % subsampling ration in space
% time_spacing=6; % subsampling ration in time (from 40Hz to 4Hz)
% 
% fprintf('compute errors TIME: 06  SPACE: 06\n');
% 
% filename_interp_space=strcat('/data/ISOTROPIC/fusion/Uinterp_spatialspacing_',num2str(space_spacing,'%.2d'),'.nc');
% filename_interp_time=strcat('/data/ISOTROPIC/fusion/Uinterp_timespacing_',num2str(time_spacing,'%.2d'),'.nc');
% filename_BF=strcat('/data/ISOTROPIC/fusion/FusedData_diagCov_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'),'.nc');
% 
% nc1=netcdf(filename_ref,'r');
% nc2 = netcdf(filename_interp_space,'r');
% nc3=netcdf(filename_interp_time,'r');
% nc4=netcdf(filename_BF,'r');
% 
% x_nonoise=nc1{'velocity_x'}(:,:,:,:);
% x_interp_spatial=nc2{'Uinterp'}(:,:,:,:);
% x_interp_temporal=nc3{'Uinterp'}(:,:,:,:);
% x_fusion_BF=nc4{'Zhat_all'}(:,:,:,:);
% 
% close(nc1); close(nc2); close(nc3); close(nc4);
% 
% err_interp_spatial_time6_space6=sqrt(sum((x_nonoise(:)-x_interp_spatial(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_interp_temporal_time6_space6=sqrt(sum((x_nonoise(:)-x_interp_temporal(:)).^2))/sqrt(sum((x_nonoise(:)).^2));
% err_fusion_time6_space6=sqrt(sum((x_nonoise(:)-x_fusion_BF(:)).^2))/sqrt(sum((x_nonoise(:)).^2));