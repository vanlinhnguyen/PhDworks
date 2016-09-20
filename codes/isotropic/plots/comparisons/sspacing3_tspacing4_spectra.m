clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;
Nt = 37;
Nz = 96; 
Ny = 96; 
Nx = 96;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
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

x_ref_all=nc1{'velocity_x'}(:,:,:,:);
x_interp_spatial_all=nc2{'Uinterp'}(:,:,:,:);
x_interp_time_all=nc3{'Uinterp'}(:,:,:,:);
x_BF_all=nc4{'Zhat_all'}(:,:,:,:);
x_LG_all=nc5{'Zhat_all'}(:,:,:,:);

x_NLM_greedy_all=permute(nc6{'x_HR_NLM_smallscales_all'}(:,:,:,:),[1,3,4,2]);
x_NLM_greedy_all=x_NLM_greedy_all+x_interp_spatial_all;

x_NLM_nongreedy_all=permute(nc7{'x_rec_all'}(:,:,:,:),[1,3,4,2]);
x_NLM_nongreedy_all=x_NLM_nongreedy_all+x_interp_spatial_all;

x_RR_all=nc8{'Urec'}(:,:,:,:);
x_KRR_all=nc9{'Urec'}(:,:,:,:);

close(nc1);  close(nc2);  close(nc3); close(nc4);
close(nc5);  close(nc6);  close(nc7); close(nc8); close(nc9);


LTHS_idt=1:time_spacing:Nx;
HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time
HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

NRMSE=@(x_ref,x_est) sqrt(sum((x_est(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));

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

%% SPECT
t_off1=[1,3]; % close to PIV

% FFT
E_ref1=zeros(Ny,1);
E_interp_space1=zeros(Ny,1);
E_interp_time1=zeros(Ny,1);
E_RR1=zeros(Ny,1);
E_KRR1=zeros(Ny,1);
E_NLM_greedy1=zeros(Ny,1);
E_NLM_nongreedy1=zeros(Ny,1);
E_BF1=zeros(Ny,1);
E_LG1=zeros(Ny,1);

E_err_interp_space1=zeros(Ny,1);
E_err_interp_time1=zeros(Ny,1);
E_err_RR1=zeros(Ny,1);
E_err_KRR1=zeros(Ny,1);
E_err_NLM_greedy1=zeros(Ny,1);
E_err_NLM_nongreedy1=zeros(Ny,1);
E_err_BF1=zeros(Ny,1);
E_err_LG1=zeros(Ny,1);

norm_factor = 1/(Nt*(numel(LTHS_idt)-1)*numel(t_off1));
for t=1:Nt
    t
    for blockid=1:numel(LTHS_idt)-1
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
               
        for pos_t=1:numel(t_off1)        
            t_current = t_PIV_prev + t_off1(pos_t);
            x_ref = squeeze(x_ref_all(t,:,:,t_current));
            x_interp_space = squeeze(x_interp_spatial_all(t,:,:,t_current));
            x_interp_time = squeeze(x_interp_time_all(t,:,:,t_current));
            x_RR = squeeze(x_RR_all(t,:,:,t_current));
            x_KRR = squeeze(x_KRR_all(t,:,:,t_current));
            x_NLM_greedy = squeeze(x_NLM_greedy_all(t,:,:,t_current));
            x_NLM_nongreedy = squeeze(x_NLM_nongreedy_all(t,:,:,t_current));
            x_BF = squeeze(x_BF_all(t,:,:,t_current));
            x_LG = squeeze(x_LG_all(t,:,:,t_current));
            
            % Spect
            E = estimate_spect_2D(x_ref); 
            E_ref1 = E_ref1 + norm_factor*E;

            E = estimate_spect_2D(x_interp_space);
            E_interp_space1 = E_interp_space1 + norm_factor*E;

            E = estimate_spect_2D(x_interp_time);
            E_interp_time1 = E_interp_time1 + norm_factor*E;

            E = estimate_spect_2D(x_RR);
            E_RR1 = E_RR1 + norm_factor*E;

            E = estimate_spect_2D(x_KRR);
            E_KRR1 = E_KRR1 + norm_factor*E;
            
            E = estimate_spect_2D(x_NLM_greedy);
            E_NLM_greedy1 = E_NLM_greedy1 + norm_factor*E;

            E = estimate_spect_2D(x_NLM_nongreedy);
            E_NLM_nongreedy1 = E_NLM_nongreedy1 + norm_factor*E;
 
            E = estimate_spect_2D(x_BF);
            E_BF1 = E_BF1 + norm_factor*E;

            E = estimate_spect_2D(x_LG);
            E_LG1 = E_LG1 + norm_factor*E;
            
            % Error
            E = estimate_spect_2D(x_ref - x_interp_space);
            E_err_interp_space1 = E_err_interp_space1 + norm_factor*E;

            E = estimate_spect_2D(x_ref - x_interp_time);
            E_err_interp_time1 = E_err_interp_time1 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_RR);
            E_err_RR1 = E_err_RR1 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_KRR);
            E_err_KRR1 = E_err_KRR1 + norm_factor*E;  
            
            E = estimate_spect_2D(x_ref - x_NLM_greedy);
            E_err_NLM_greedy1 = E_err_NLM_greedy1 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_NLM_nongreedy);
            E_err_NLM_nongreedy1 = E_err_NLM_nongreedy1 + norm_factor*E;              
            
            E = estimate_spect_2D(x_ref - x_BF);
            E_err_BF1 = E_err_BF1 + norm_factor*E;              

            E = estimate_spect_2D(x_ref - x_LG);
            E_err_LG1 = E_err_LG1 + norm_factor*E;                         
        end 
    end 
end 



%%
t_off2=[2]; % close to PIV

% FFT
E_ref2=zeros(Ny,1);
E_interp_space2=zeros(Ny,1);
E_interp_time2=zeros(Ny,1);
E_RR2=zeros(Ny,1);
E_KRR2=zeros(Ny,1);
E_NLM_greedy2=zeros(Ny,1);
E_NLM_nongreedy2=zeros(Ny,1);
E_BF2=zeros(Ny,1);
E_LG2=zeros(Ny,1);

E_err_interp_space2=zeros(Ny,1);
E_err_interp_time2=zeros(Ny,1);
E_err_RR2=zeros(Ny,1);
E_err_KRR2=zeros(Ny,1);
E_err_NLM_greedy2=zeros(Ny,1);
E_err_NLM_nongreedy2=zeros(Ny,1);
E_err_BF2=zeros(Ny,1);
E_err_LG2=zeros(Ny,1);

norm_factor = 1/(Nt*(numel(LTHS_idt)-1)*numel(t_off2));
for t=1:Nt
    t
    for blockid=1:numel(LTHS_idt)-1
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
               
        for pos_t=1:numel(t_off2)        
            t_current = t_PIV_prev + t_off2(pos_t);
            x_ref = squeeze(x_ref_all(t,:,:,t_current));
            x_interp_space = squeeze(x_interp_spatial_all(t,:,:,t_current));
            x_interp_time = squeeze(x_interp_time_all(t,:,:,t_current));
            x_RR = squeeze(x_RR_all(t,:,:,t_current));
            x_KRR = squeeze(x_KRR_all(t,:,:,t_current));
            x_NLM_greedy = squeeze(x_NLM_greedy_all(t,:,:,t_current));
            x_NLM_nongreedy = squeeze(x_NLM_nongreedy_all(t,:,:,t_current));
            x_BF = squeeze(x_BF_all(t,:,:,t_current));
            x_LG = squeeze(x_LG_all(t,:,:,t_current));
            
            % Spect
            E = estimate_spect_2D(x_ref); 
            E_ref2 = E_ref2 + norm_factor*E;

            E = estimate_spect_2D(x_interp_space);
            E_interp_space2 = E_interp_space2 + norm_factor*E;

            E = estimate_spect_2D(x_interp_time);
            E_interp_time2 = E_interp_time2 + norm_factor*E;

            E = estimate_spect_2D(x_RR);
            E_RR2 = E_RR2 + norm_factor*E;

            E = estimate_spect_2D(x_KRR);
            E_KRR2 = E_KRR2 + norm_factor*E;
            
            E = estimate_spect_2D(x_NLM_greedy);
            E_NLM_greedy2 = E_NLM_greedy2 + norm_factor*E;

            E = estimate_spect_2D(x_NLM_nongreedy);
            E_NLM_nongreedy2 = E_NLM_nongreedy2 + norm_factor*E;
 
            E = estimate_spect_2D(x_BF);
            E_BF2 = E_BF2 + norm_factor*E;

            E = estimate_spect_2D(x_LG);
            E_LG2 = E_LG2 + norm_factor*E;
            
            % Error
            E = estimate_spect_2D(x_ref - x_interp_space);
            E_err_interp_space2 = E_err_interp_space2 + norm_factor*E;

            E = estimate_spect_2D(x_ref - x_interp_time);
            E_err_interp_time2 = E_err_interp_time2 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_RR);
            E_err_RR2 = E_err_RR2 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_KRR);
            E_err_KRR2 = E_err_KRR2 + norm_factor*E;  
            
            E = estimate_spect_2D(x_ref - x_NLM_greedy);
            E_err_NLM_greedy2 = E_err_NLM_greedy2 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_NLM_nongreedy);
            E_err_NLM_nongreedy2 = E_err_NLM_nongreedy2 + norm_factor*E;              
            
            E = estimate_spect_2D(x_ref - x_BF);
            E_err_BF2 = E_err_BF2 + norm_factor*E;              

            E = estimate_spect_2D(x_ref - x_LG);
            E_err_LG2 = E_err_LG2 + norm_factor*E;              
        end 
    end 
end 


%% Plot spectra
fsize=20;
fname='CMU Serif';
cc=[96/255 96/255 96/255];

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
ax(1)=subplot(2,2,1,'position',[0.14 0.575 0.325 0.35]); % top left

% -5/3 line
xline=[4,20];
yline=5*[exp((-5/3)*log(4)),exp((-5/3)*log(20))];

h11=loglog(3:k_max,E_ref1(3:k_max),'k-','LineWidth',2); 
hold on;
h12=loglog(3:k_max,E_interp_space1(3:k_max),'b-','LineWidth',2); 
h13=loglog(3:k_max,E_interp_time1(3:k_max),'g-','LineWidth',2);
h14=loglog(3:k_max,E_KRR1(3:k_max),'m-','LineWidth',2);
h15=loglog(3:k_max,E_NLM_nongreedy1(3:k_max),'color',cc,'LineWidth',2); 
h16=loglog(3:k_max,E_BF1(3:k_max),'r-','LineWidth',2); 
loglog(xline,yline,'k-','LineWidth',2)
plot([k_cutoff,k_cutoff] ,[1e-8,1e-2],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([3*1e-6,1]);
text(15,0.5,'-5/3','HorizontalAlignment','right','FontSize',fsize)
text(6,1e-5,'$t/\delta t = 1$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);
set (gca, 'XTickLabel', {'',''},'FontSize',fsize)

ylabel('$E(k)$','interpreter','latex');

box on

% =======    
ax(3)=axes('position',[0.52 0.575 0.325 0.35]); % top right
% -5/3 line
xline=[4,20];
yline=5*[exp((-5/3)*log(4)),exp((-5/3)*log(20))];

h21=loglog(3:k_max,E_ref2(3:k_max),'k-','LineWidth',2); 
hold on;
h22=loglog(3:k_max,E_interp_space2(3:k_max),'b-','LineWidth',2); 
h23=loglog(3:k_max,E_interp_time2(3:k_max),'g-','LineWidth',2);
h24=loglog(3:k_max,E_KRR2(3:k_max),'m-','LineWidth',2);
h25=loglog(3:k_max,E_NLM_nongreedy2(3:k_max),'color',cc,'LineWidth',2); 
h26=loglog(3:k_max,E_BF2(3:k_max),'r-','LineWidth',2); 
loglog(xline,yline,'k-','LineWidth',2)
plot([k_cutoff,k_cutoff] ,[1e-8,1e-2],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([3*1e-6,1]);
text(15,0.5,'-5/3','HorizontalAlignment','right','FontSize',fsize)
text(6,1e-5,'$t/\delta t = 2$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);
set (gca, 'YTickLabel', {'','','','',''},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);
set (gca, 'XTickLabel', {'',''},'FontSize',fsize)

box on;

% =======    
ax(2)=axes('position',[0.14 0.175 0.325 0.35]); % bottom left

h31=loglog(3:k_max,E_err_interp_space1(3:k_max)./E_ref1(3:k_max),'b-','LineWidth',2); 
hold on;
h32=loglog(3:k_max,E_err_interp_time1(3:k_max)./E_ref1(3:k_max),'g-','LineWidth',2); 
h33=loglog(3:k_max,E_err_KRR1(3:k_max)./E_ref1(3:k_max),'m-','LineWidth',2); 
h34=loglog(3:k_max,E_err_NLM_nongreedy1(3:k_max)./E_ref1(3:k_max),'color',cc,'LineWidth',2); 
h35=loglog(3:k_max,E_err_BF1(3:k_max)./E_ref1(3:k_max),'r-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([1e-4,2]);

text(6,3*1e-1,'$t/\delta t = 1$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);

xlabel('$k$','Interpreter','Latex');
ylabel('$E^{\epsilon}(k)/E(k)$','Interpreter','Latex');

box on;

% =======    

ax(4)=axes('position',[0.52 0.175 0.325 0.35]);
h41=loglog(3:k_max,E_err_interp_space2(3:k_max)./E_ref2(3:k_max),'b-','LineWidth',2); 
hold on;
h42=loglog(3:k_max,E_err_interp_time2(3:k_max)./E_ref2(3:k_max),'g-','LineWidth',2); 
h43=loglog(3:k_max,E_err_KRR1(3:k_max)./E_ref2(3:k_max),'m-','LineWidth',2); 
h44=loglog(3:k_max,E_err_NLM_nongreedy2(3:k_max)./E_ref2(3:k_max),'color',cc,'LineWidth',2); 
h45=loglog(3:k_max,E_err_BF2(3:k_max)./E_ref2(3:k_max),'r-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([2,60]);
ylim([1e-4,2]);

text(6,3*1e-1,'$t/\delta t = 2$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', [1e0,1e1]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', [1e-5,1e-4,1e-3,1e-2,1e-1,1e0]);
set (gca, 'YTickLabel', {'','','','',''},'FontSize',fsize)

xlabel('$k$','Interpreter','Latex');
box on;

leg={'Reference','$\mathbf{I}_s \mathbf{y}$','$\mathbf{I}_t\mathbf{x}$','KRR','Nongreedy propag','Bayesian fusion'};
legendflex([h11,h12,h13,h14,h15,h16], leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 0], ...
                       'nrow',2, ...
                       'fontsize',fsize,...
                       'Interpreter','Latex',...
                       'box','off');
export_fig('./figures/spectra_sspacing3_tspacing4_spectra2d','-a1','-q101','-eps','-painters');
close()     