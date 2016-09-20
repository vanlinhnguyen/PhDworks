% clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
Nt = 37;
Nz = 96;
Ny = 96; 
Nx = 96;

LTHS_idt=1:time_spacing:Nx;
HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time
HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

sborder=4*space_spacing;
sleft=sborder; sright = sborder-(Ny-HTLS_idy(end));
sbottom=sborder; stop = sborder-(Nz-HTLS_idz(end));

tborder=4*time_spacing;
tleft=tborder; tright = tborder-(Nx-LTHS_idt(end));

HTHS_idt_enlarged = 1:Nx+tleft+tright;
LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);

k_max=Ny/2;
k_1D_HR=[0:k_max -k_max+1:1:-1];
[KZS,KYS] = meshgrid(k_1D_HR,k_1D_HR); 
k_cutoff = k_max/space_spacing;

NRMSE=@(x_ref,x_est) sqrt(sum((x_est(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));

%% Compute average NRMSE
% % ================================ Vary tau ==============================%
% filename_NLmean_tau1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0050.nc';
% filename_NLmean_tau2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% filename_NLmean_tau3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0200.nc';
% filename_NLmean_tau4='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0400.nc';
% 
% nc1  = netcdf(filename_ref,'r');
% nc21 = netcdf(filename_NLmean_tau1,'r');
% nc22 = netcdf(filename_NLmean_tau2,'r');
% nc23 = netcdf(filename_NLmean_tau3,'r');
% nc24 = netcdf(filename_NLmean_tau4,'r');
% 
% NRMSE_NLmeans_tau1=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_NLmeans_tau2=NRMSE_NLmeans_tau1;
% NRMSE_NLmeans_tau3=NRMSE_NLmeans_tau1;
% NRMSE_NLmeans_tau4=NRMSE_NLmeans_tau1;
% 
% for t=1:Nt
%     t
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:time_spacing-1           
%             t_current = t_PIV_prev + pos_t;
%             x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
% 
% %             x_rec = x_interp_space+nc2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_tau1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc22{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_tau2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc23{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_tau3(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc24{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_tau4(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);            
%         end 
%     end 
% end 
% close(nc1);close(nc21); close(nc22); close(nc23); close(nc24);
% 
% NRMSE_NLmeans_tau1_mean=mean(NRMSE_NLmeans_tau1,2);
% NRMSE_NLmeans_tau1_std=std(NRMSE_NLmeans_tau1,0,2);
% 
% NRMSE_NLmeans_tau2_mean=mean(NRMSE_NLmeans_tau2,2);
% NRMSE_NLmeans_tau2_std=std(NRMSE_NLmeans_tau2,0,2);
% 
% NRMSE_NLmeans_tau3_mean=mean(NRMSE_NLmeans_tau3,2);
% NRMSE_NLmeans_tau3_std=std(NRMSE_NLmeans_tau3,0,2);
% 
% NRMSE_NLmeans_tau4_mean=mean(NRMSE_NLmeans_tau4,2);
% NRMSE_NLmeans_tau4_std=std(NRMSE_NLmeans_tau4,0,2);
% 
% % ============================== Vary sim size ===========================%
% filename_NLmean_sim1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim04_acc04_neighbor5_tau0100.nc';
% filename_NLmean_sim2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim08_acc08_neighbor5_tau0100.nc';
% filename_NLmean_sim3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% 
% nc1  = netcdf(filename_ref,'r');
% nc21 = netcdf(filename_NLmean_sim1,'r');
% nc22 = netcdf(filename_NLmean_sim2,'r');
% nc23 = netcdf(filename_NLmean_sim3,'r');
% 
% NRMSE_NLmeans_sim1=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_NLmeans_sim2=NRMSE_NLmeans_sim1;
% NRMSE_NLmeans_sim3=NRMSE_NLmeans_sim1;
% 
% for t=1:Nt
%     t
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:time_spacing-1           
%             t_current = t_PIV_prev + pos_t;
%             x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
% 
% %             x_rec = x_interp_space+nc2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_sim1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc22{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_sim2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc23{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_sim3(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%         end 
%     end 
% end 
% close(nc1);close(nc21); close(nc22); close(nc23);
% 
% NRMSE_NLmeans_sim1_mean=mean(NRMSE_NLmeans_sim1,2);
% NRMSE_NLmeans_sim1_std=std(NRMSE_NLmeans_sim1,0,2);
% 
% NRMSE_NLmeans_sim2_mean=mean(NRMSE_NLmeans_sim2,2);
% NRMSE_NLmeans_sim2_std=std(NRMSE_NLmeans_sim2,0,2);
% 
% NRMSE_NLmeans_sim3_mean=mean(NRMSE_NLmeans_sim3,2);
% NRMSE_NLmeans_sim3_std=std(NRMSE_NLmeans_sim3,0,2);
% 
% % ============================= Vary acc size ============================%
% filename_NLmean_acc1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc00_neighbor5_tau0100.nc';
% filename_NLmean_acc2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc01_neighbor5_tau0100.nc';
% filename_NLmean_acc3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc06_neighbor5_tau0100.nc';
% filename_NLmean_acc4='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% 
% nc1  = netcdf(filename_ref,'r');
% nc21 = netcdf(filename_NLmean_acc1,'r');
% nc22 = netcdf(filename_NLmean_acc2,'r');
% nc23 = netcdf(filename_NLmean_acc3,'r');
% nc24 = netcdf(filename_NLmean_acc4,'r');
% 
% NRMSE_NLmeans_acc1=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_NLmeans_acc2=NRMSE_NLmeans_acc1;
% NRMSE_NLmeans_acc3=NRMSE_NLmeans_acc1;
% NRMSE_NLmeans_acc4=NRMSE_NLmeans_acc1;
% 
% for t=1:Nt
%     t
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:time_spacing-1           
%             t_current = t_PIV_prev + pos_t;
%             x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
% 
% %             x_rec = x_interp_space+nc2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_acc1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc22{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_acc2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc23{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_acc3(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc24{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_acc4(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);            
%         end 
%     end 
% end 
% close(nc1);close(nc21); close(nc22); close(nc23); close(nc24);
% 
% NRMSE_NLmeans_acc1_mean=mean(NRMSE_NLmeans_acc1,2);
% NRMSE_NLmeans_acc1_std=std(NRMSE_NLmeans_acc1,0,2);
% 
% NRMSE_NLmeans_acc2_mean=mean(NRMSE_NLmeans_acc2,2);
% NRMSE_NLmeans_acc2_std=std(NRMSE_NLmeans_acc2,0,2);
% 
% NRMSE_NLmeans_acc3_mean=mean(NRMSE_NLmeans_acc3,2);
% NRMSE_NLmeans_acc3_std=std(NRMSE_NLmeans_acc3,0,2);
% 
% NRMSE_NLmeans_acc4_mean=mean(NRMSE_NLmeans_acc4,2);
% NRMSE_NLmeans_acc4_std=std(NRMSE_NLmeans_acc4,0,2);
% 
% % =========================== Vary neighbor size =========================%
% filename_NLmean_Ns1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor2_tau0100.nc';
% filename_NLmean_Ns2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% filename_NLmean_Ns3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor8_tau0100.nc';
% 
% nc1  = netcdf(filename_ref,'r');
% nc21 = netcdf(filename_NLmean_Ns1,'r');
% nc22 = netcdf(filename_NLmean_Ns2,'r');
% nc23 = netcdf(filename_NLmean_Ns3,'r');
% 
% NRMSE_NLmeans_Ns1=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_NLmeans_Ns2=NRMSE_NLmeans_Ns1;
% NRMSE_NLmeans_Ns3=NRMSE_NLmeans_Ns1;
% 
% for t=1:Nt
%     t
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:time_spacing-1           
%             t_current = t_PIV_prev + pos_t;
%             x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
% 
% %             x_rec = x_interp_space+nc2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_Ns1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
% 
%             x_rec = x_interp_space+nc22{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_Ns2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%  
%             x_rec = x_interp_space+nc23{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans_Ns3(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%         end 
%     end 
% end 
% close(nc1);close(nc21); close(nc22); close(nc23);
% 
% NRMSE_NLmeans_Ns1_mean=mean(NRMSE_NLmeans_Ns1,2);
% NRMSE_NLmeans_Ns1_std=std(NRMSE_NLmeans_Ns1,0,2);
% 
% NRMSE_NLmeans_Ns2_mean=mean(NRMSE_NLmeans_Ns2,2);
% NRMSE_NLmeans_Ns2_std=std(NRMSE_NLmeans_Ns2,0,2);
% 
% NRMSE_NLmeans_Ns3_mean=mean(NRMSE_NLmeans_Ns3,2);
% NRMSE_NLmeans_Ns3_std=std(NRMSE_NLmeans_Ns3,0,2);

%% SPECTRA
% t_off=[2]; % close to PIV
% 
% nc_ref  = netcdf(filename_ref,'r');
% 
% filename_NLmean_tau1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0050.nc';
% filename_NLmean_tau2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% filename_NLmean_tau3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0200.nc';
% filename_NLmean_tau4='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0400.nc';
% nc_NLmean_tau1 = netcdf(filename_NLmean_tau1,'r');
% nc_NLmean_tau2 = netcdf(filename_NLmean_tau2,'r');
% nc_NLmean_tau3 = netcdf(filename_NLmean_tau3,'r');
% nc_NLmean_tau4 = netcdf(filename_NLmean_tau4,'r');
% 
% filename_NLmean_Ns1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor2_tau0100.nc';
% filename_NLmean_Ns2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% filename_NLmean_Ns3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor8_tau0100.nc';
% nc_NLmean_Ns1 = netcdf(filename_NLmean_Ns1,'r');
% nc_NLmean_Ns2 = netcdf(filename_NLmean_Ns2,'r');
% nc_NLmean_Ns3 = netcdf(filename_NLmean_Ns3,'r');
% 
% 
% filename_NLmean_sim1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim04_acc04_neighbor5_tau0100.nc';
% filename_NLmean_sim2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim08_acc08_neighbor5_tau0100.nc';
% filename_NLmean_sim3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% nc_NLmean_sim1 = netcdf(filename_NLmean_sim1,'r');
% nc_NLmean_sim2 = netcdf(filename_NLmean_sim2,'r');
% nc_NLmean_sim3 = netcdf(filename_NLmean_sim3,'r');
% 
% filename_NLmean_acc1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc00_neighbor5_tau0100.nc';
% filename_NLmean_acc2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc01_neighbor5_tau0100.nc';
% filename_NLmean_acc3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc06_neighbor5_tau0100.nc';
% filename_NLmean_acc4='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% nc_NLmean_acc1 = netcdf(filename_NLmean_acc1,'r');
% nc_NLmean_acc2 = netcdf(filename_NLmean_acc2,'r');
% nc_NLmean_acc3 = netcdf(filename_NLmean_acc3,'r');
% nc_NLmean_acc4 = netcdf(filename_NLmean_acc4,'r');
% 
% E_ref=zeros(Ny,1);
% E_err_NLmean_tau1=zeros(Ny,1); E_err_NLmean_tau2=zeros(Ny,1); E_err_NLmean_tau3=zeros(Ny,1); E_err_NLmean_tau4=zeros(Ny,1);
% E_err_NLmean_Ns1=zeros(Ny,1); E_err_NLmean_Ns2=zeros(Ny,1); E_err_NLmean_Ns3=zeros(Ny,1);
% E_err_NLmean_sim1=zeros(Ny,1); E_err_NLmean_sim2=zeros(Ny,1); E_err_NLmean_sim3=zeros(Ny,1);
% E_err_NLmean_acc1=zeros(Ny,1); E_err_NLmean_acc2=zeros(Ny,1); E_err_NLmean_acc3=zeros(Ny,1); E_err_NLmean_acc4=zeros(Ny,1);
% 
% norm_factor = 1/(Nt*(numel(LTHS_idt)-1)*numel(t_off));
% for t=1:Nt
%     t
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:numel(t_off)        
%             t_current = t_PIV_prev + t_off(pos_t);
%             x_ref = nc_ref{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
% 
%             E = estimate_spect_2D(x_ref); 
%             E_ref = E_ref + norm_factor*E;
%             
%             %======== vary tau =========%
%             x_NLmean_tau1 = x_interp_space+nc_NLmean_tau1{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_tau1 = E_err_NLmean_tau1 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_tau1);           
%             
%             x_NLmean_tau2 = x_interp_space+nc_NLmean_tau2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_tau2 = E_err_NLmean_tau2 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_tau2);  
%             
%             x_NLmean_tau3 = x_interp_space+nc_NLmean_tau3{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_tau3 = E_err_NLmean_tau3 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_tau3);  
%             
%             x_NLmean_tau4 = x_interp_space+nc_NLmean_tau4{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_tau4 = E_err_NLmean_tau4 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_tau4);  
%             
%             %======== vary Ns =========%
%             x_NLmean_Ns1 = x_interp_space+nc_NLmean_Ns1{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_Ns1 = E_err_NLmean_Ns1 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_Ns1);           
%             
%             x_NLmean_Ns2 = x_interp_space+nc_NLmean_Ns2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_Ns2 = E_err_NLmean_Ns2 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_Ns2);           
%             
%             x_NLmean_Ns3 = x_interp_space+nc_NLmean_Ns3{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_Ns3 = E_err_NLmean_Ns3 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_Ns3);           
%             
%             %======== vary sim patch size =========%
%             x_NLmean_sim1 = x_interp_space+nc_NLmean_sim1{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_sim1 = E_err_NLmean_sim1 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_sim1);           
%             
%             x_NLmean_sim2 = x_interp_space+nc_NLmean_sim2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_sim2 = E_err_NLmean_sim2 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_sim2);           
%             
%             x_NLmean_sim3 = x_interp_space+nc_NLmean_sim3{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_sim3 = E_err_NLmean_sim3 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_sim3);           
% 
%             %======== vary tau =========%
%             x_NLmean_acc1 = x_interp_space+nc_NLmean_acc1{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_acc1 = E_err_NLmean_acc1 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_acc1);           
%             
%             x_NLmean_acc2 = x_interp_space+nc_NLmean_acc2{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_acc2 = E_err_NLmean_acc2 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_acc2);           
%             
%             x_NLmean_acc3 = x_interp_space+nc_NLmean_acc3{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_acc3 = E_err_NLmean_acc3 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_acc3);           
%             
%             x_NLmean_acc4 = x_interp_space+nc_NLmean_acc4{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             E_err_NLmean_acc4 = E_err_NLmean_acc4 + norm_factor*estimate_spect_2D(x_ref - x_NLmean_acc4);           
%         end 
%     end 
% end 
% close(nc_ref);
% close(nc_NLmean_tau1); close(nc_NLmean_tau2); close(nc_NLmean_tau3); close(nc_NLmean_tau4);
% close(nc_NLmean_Ns1); close(nc_NLmean_Ns2); close(nc_NLmean_Ns3);
% close(nc_NLmean_sim1); close(nc_NLmean_sim2); close(nc_NLmean_sim3); 
% close(nc_NLmean_acc1); close(nc_NLmean_acc2); close(nc_NLmean_acc3); close(nc_NLmean_acc4);
% 
% save NLmean_propag1dir_sspacing3_tspacing4_allparams.mat ...
%      NRMSE_NLmeans_tau1_mean NRMSE_NLmeans_tau1_std ...
%      NRMSE_NLmeans_tau2_mean NRMSE_NLmeans_tau2_std ...
%      NRMSE_NLmeans_tau3_mean NRMSE_NLmeans_tau3_std ...
%      NRMSE_NLmeans_tau4_mean NRMSE_NLmeans_tau4_std ...
%      NRMSE_NLmeans_sim1_mean NRMSE_NLmeans_sim1_std ...
%      NRMSE_NLmeans_sim2_mean NRMSE_NLmeans_sim2_std ...
%      NRMSE_NLmeans_sim3_mean NRMSE_NLmeans_sim3_std ...
%      NRMSE_NLmeans_acc1_mean NRMSE_NLmeans_acc1_std ...
%      NRMSE_NLmeans_acc2_mean NRMSE_NLmeans_acc2_std ...
%      NRMSE_NLmeans_acc3_mean NRMSE_NLmeans_acc3_std ...
%      NRMSE_NLmeans_acc4_mean NRMSE_NLmeans_acc4_std ...
%      NRMSE_NLmeans_Ns1_mean NRMSE_NLmeans_Ns1_std ...
%      NRMSE_NLmeans_Ns2_mean NRMSE_NLmeans_Ns2_std ...
%      NRMSE_NLmeans_Ns3_mean NRMSE_NLmeans_Ns3_std ...
%      E_ref...
%      E_err_NLmean_tau1 E_err_NLmean_tau2 ...
%      E_err_NLmean_tau3 E_err_NLmean_tau4 ...
%      E_err_NLmean_Ns1 E_err_NLmean_Ns2 E_err_NLmean_Ns3...
%      E_err_NLmean_sim1 E_err_NLmean_sim2 E_err_NLmean_sim3 ...
%      E_err_NLmean_acc1 E_err_NLmean_acc2 ...
%      E_err_NLmean_acc3 E_err_NLmean_acc4;
 
%% Load precomputed data
load  NLmean_propag1dir_sspacing3_tspacing4_allparams.mat 

%%
NRMSE_mean_alltau = [NRMSE_NLmeans_tau1_mean NRMSE_NLmeans_tau2_mean NRMSE_NLmeans_tau3_mean NRMSE_NLmeans_tau4_mean];
NRMSE_std_alltau = [NRMSE_NLmeans_tau1_std NRMSE_NLmeans_tau2_std NRMSE_NLmeans_tau3_std NRMSE_NLmeans_tau4_std];

NRMSE_mean_allNs = [NRMSE_NLmeans_Ns1_mean NRMSE_NLmeans_Ns2_mean NRMSE_NLmeans_Ns3_mean];
NRMSE_std_allNs = [NRMSE_NLmeans_Ns1_std NRMSE_NLmeans_Ns2_std NRMSE_NLmeans_Ns3_std];

NRMSE_mean_allsim = [NRMSE_NLmeans_sim1_mean NRMSE_NLmeans_sim2_mean NRMSE_NLmeans_sim3_mean];
NRMSE_std_allsim = [NRMSE_NLmeans_sim1_std NRMSE_NLmeans_sim2_std NRMSE_NLmeans_sim3_std];

NRMSE_mean_allacc = [NRMSE_NLmeans_acc1_mean NRMSE_NLmeans_acc2_mean NRMSE_NLmeans_acc3_mean NRMSE_NLmeans_acc4_mean];
NRMSE_std_allacc = [NRMSE_NLmeans_acc1_std NRMSE_NLmeans_acc2_std NRMSE_NLmeans_acc3_std NRMSE_NLmeans_acc4_std];

%% Plot NRMSE
fsize=20;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 1200 1200]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% 
ax(1)=subplot(2,2,1,'position',[0.14 0.6 0.325 0.35]); % top left
h = barwitherr(NRMSE_std_alltau, NRMSE_mean_alltau);% Plot with errorbars

set(h(1), 'FaceColor',[0.9 0.9 0.9]) 
set(h(2), 'FaceColor',[0.7 0.7 0.7]) 
set(h(3), 'FaceColor',[0.5 0.5 0.5]) 
set(h(4), 'FaceColor',[0.3 0.3 0.3]) 

ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
ylim([0 0.22]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.2);
set (gca, 'YTickLabel', {'0.0','','0.1','','0.2'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 1:1:3);
set (gca, 'XTickLabel', {'','',''},'FontSize',fsize)

tex1 = text('String', '$\tau=0.05,0.1,0.2,0.4$', ...
        'interpreter', 'latex',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+6, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [2 0.2 0]);

% =======    
ax(3)=axes('position',[0.52 0.6 0.325 0.35]); % top right
h = barwitherr(NRMSE_std_allNs, NRMSE_mean_allNs);% Plot with errorbars
set(h(1), 'FaceColor',[0.7 0.7 0.7]) 
set(h(2), 'FaceColor',[0.5 0.5 0.5]) 
set(h(3), 'FaceColor',[0.3 0.3 0.3]) 

ylim([0 0.22]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.2);
set (gca, 'YTickLabel', {'','','','',''},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 1:1:3);
set (gca, 'XTickLabel', {'','',''},'FontSize',fsize)

box on;

tex3 = text('String', '$\sqrt{r}=5,11,17$', ...
        'interpreter', 'latex',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+6, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [2 0.2 0]);
% =======    
ax(2)=axes('position',[0.14 0.175 0.325 0.35]); % bottom left

h = barwitherr(NRMSE_std_allsim, NRMSE_mean_allsim);% Plot with errorbars

set(h(1), 'FaceColor',[0.7 0.7 0.7]) 
set(h(2), 'FaceColor',[0.5 0.5 0.5]) 
set(h(3), 'FaceColor',[0.3 0.3 0.3]) 

xlabel('$t/\delta t$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
ylim([0 0.22]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.2);
set (gca, 'YTickLabel', {'0.0','','0.1','','0.2'},'FontSize',fsize)
box on;

tex2 = text('String', '$\sqrt{p}=9,17,25$', ...
        'interpreter', 'latex',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+6, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [2 0.2 0]);
% =======    

ax(4)=axes('position',[0.52 0.175 0.325 0.35]);
h = barwitherr(NRMSE_std_allacc, NRMSE_mean_allacc);% Plot with errorbars

set(h(1), 'FaceColor',[0.9 0.9 0.9]) 
set(h(2), 'FaceColor',[0.7 0.7 0.7]) 
set(h(3), 'FaceColor',[0.5 0.5 0.5]) 
set(h(4), 'FaceColor',[0.3 0.3 0.3])

xlabel('$t/\delta t$','Interpreter','latex')
ylim([0 0.22]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.05:0.2);
set (gca, 'YTickLabel', {'','','','',''},'FontSize',fsize)

box on;

tex4 = text('String', '$\sqrt{q}=1,3,13,25$', ...
        'interpreter', 'latex',...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'Baseline', ...
        'FontUnits', 'pixels', ...
        'FontSize', fsize+6, ...
        'FontWeight', 'Bold', ...
        'FontName', fname, ...
        'Position', [2 0.2 0]);
 
export_fig('./figures/NLmean_propag1dir_spacing3_tspacing4_allparams_NRMSE', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()
 
%% Plot spectra
fsize=20;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 1200 1000]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% 
ax(1)=subplot(2,2,1,'position',[0.14 0.6 0.325 0.35]); % top left

h11=loglog(3:k_max,E_err_NLmean_tau1(3:k_max)./E_ref(3:k_max),'g-','LineWidth',2); 
hold on;
h12=loglog(3:k_max,E_err_NLmean_tau2(3:k_max)./E_ref(3:k_max),'m-','LineWidth',2); 
h13=loglog(3:k_max,E_err_NLmean_tau3(3:k_max)./E_ref(3:k_max),'r-','LineWidth',2); 
h14=loglog(3:k_max,E_err_NLmean_tau4(3:k_max)./E_ref(3:k_max),'b-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([1e-4,15]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-4,1e-3,1e-2,1e-1,1e0]);

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);
set (gca, 'XTickLabel', {''},'FontSize',fsize)

ylabel('$E^{\epsilon}(k)/E(k)$','Interpreter','Latex');

box on

% =======    
ax(3)=axes('position',[0.52 0.6 0.325 0.35]); % top right
h21=loglog(3:k_max,E_err_NLmean_Ns1(3:k_max)./E_ref(3:k_max),'m-','LineWidth',2); 
hold on;
h22=loglog(3:k_max,E_err_NLmean_Ns2(3:k_max)./E_ref(3:k_max),'r-','LineWidth',2); 
h23=loglog(3:k_max,E_err_NLmean_Ns3(3:k_max)./E_ref(3:k_max),'b-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([1e-4,15]);

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);
set (gca, 'XTickLabel', {''},'FontSize',fsize)

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-4,1e-3,1e-2,1e-1,1e0]);
set (gca, 'YTickLabel', {'','','','',''},'FontSize',fsize)


box on;

% =======    
ax(2)=axes('position',[0.14 0.2 0.325 0.35]); % bottom left
h31=loglog(3:k_max,E_err_NLmean_sim1(3:k_max)./E_ref(3:k_max),'m-','LineWidth',2); 
hold on;
h32=loglog(3:k_max,E_err_NLmean_sim2(3:k_max)./E_ref(3:k_max),'r-','LineWidth',2); 
h33=loglog(3:k_max,E_err_NLmean_sim3(3:k_max)./E_ref(3:k_max),'b-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([1e-4,15]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-4,1e-3,1e-2,1e-1,1e0]);

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);
% set (gca, 'XTickLabel', {''},'FontSize',fsize)

xlabel('$k$','Interpreter','Latex'); 
ylabel('$E^{\epsilon}(k)/E(k)$','interpreter','latex');
box on;


% =======    

ax(4)=axes('position',[0.52 0.2 0.325 0.35]);
h41=loglog(3:k_max,E_err_NLmean_acc1(3:k_max)./E_ref(3:k_max),'g-','LineWidth',2); 
hold on;
h42=loglog(3:k_max,E_err_NLmean_acc2(3:k_max)./E_ref(3:k_max),'m-','LineWidth',2); 
h43=loglog(3:k_max,E_err_NLmean_acc3(3:k_max)./E_ref(3:k_max),'r-','LineWidth',2); 
h44=loglog(3:k_max,E_err_NLmean_acc4(3:k_max)./E_ref(3:k_max),'b-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([1e-4,15]);

xlabel('$k$','Interpreter','Latex');  
set(gca, 'YTick', []);
box on;

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-4,1e-3,1e-2,1e-1,1e0]);
set (gca, 'YTickLabel', {'','','','',''},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);
% set (gca, 'XTickLabel', {''},'FontSize',fsize)

leg={'$\tau=0.05$','$\tau=0.1$','$\tau=0.2$','$\tau=0.4$'};
hleg=legendflex([h11,h12,h13,h14],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[-250 -60], ...
                       'nrow',2, ...
                       'fontsize',fsize,...
                       'Interpreter','Latex',...
                       'box','off');
                   
leg={'$\sqrt{r}=5$','$\sqrt{r}=11$','$\sqrt{r}=17$'};
legendflex([h21,h22,h23],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[220 -60], ...
                       'nrow',2, ...
                       'fontsize',fsize, ...
                       'Interpreter','Latex',...
                       'box','off');

leg={'$\sqrt{q}=1$','$\sqrt{q}=3$','$\sqrt{q}=13$','$\sqrt{q}=25$'};
legendflex([h41,h42,h43,h44],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[220 -460], ...
                       'nrow',2, ...
                       'fontsize',fsize, ...
                       'Interpreter','Latex',...
                       'box','off');    
                   
leg={'$\sqrt{p}=9$','$\sqrt{p}=17$','$\sqrt{p}=25$'};
legendflex([h31,h32,h33],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[-240 -460], ...
                       'nrow',2, ...
                       'fontsize',fsize, ...
                       'Interpreter','Latex',...
                       'box','off');                       
% stupid trick: comment exporting first, double click 2 times in the plot
% to make the legend appear, then export. Avoid this with single plot
% instead of subplot!
export_fig('./figures/NLmean_propag1dir_spacing3_tspacing4_allparams_spectra', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()
