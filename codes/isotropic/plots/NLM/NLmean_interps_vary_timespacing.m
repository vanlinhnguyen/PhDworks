% clear all; close all; clc;

%% COMMON PARAMS

Nt = 37;
Nz = 96;
Ny = 96; 
Nx = 96;
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time

NRMSE=@(x_ref,x_est) sqrt(sum((x_est(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));

%% SPACE 3 TIME 4
% space_spacing=3;
% time_spacing=4;
% 
% filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
% filename_NLM1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% filename_NLM2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
% 
% HTHS_idy=1:Ny; % row indices of HTHS in space 
% HTHS_idz=1:Nz; % column indices of HTHS in space
% HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
% HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space
% 
% LTHS_idt=1:time_spacing:Nx;
% 
% sborder=4*space_spacing;
% tborder=4*time_spacing;
% tleft=tborder; tright = tborder-(Nx-LTHS_idt(end));
% HTHS_idt_enlarged = 1:Nx+tleft+tright;
% LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);
% 
% nc1=netcdf(filename_ref,'r');
% nc21=netcdf(filename_NLM1,'r');
% nc22=netcdf(filename_NLM2,'r');
% 
% NRMSE_space3_time4_interp_space=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_space3_time4_interp_time=NRMSE_space3_time4_interp_space;
% NRMSE_space3_time4_NLmeans1=NRMSE_space3_time4_interp_space;
% NRMSE_space3_time4_NLmeans2=NRMSE_space3_time4_interp_space;
% for t=1:Nt
%     t
%     PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
%     PIV_sampled=cat(3, PIV_sampled(:,:,Nx-tborder+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:tright));    
%     LTHS_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
%     LTHS_interp=permute(LTHS_interp(tborder+1:tborder+Nx,:,:),[2,3,1]);
% 
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:time_spacing-1           
%             t_current = t_PIV_prev + pos_t;
%             x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
%             x_interp_time = squeeze(LTHS_interp(:,:,t_current));
% 
%             x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_space3_time4_NLmeans1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%  
%             x_rec = x_interp_space+nc22{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_space3_time4_NLmeans2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%             
%             NRMSE_space3_time4_interp_space(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_space);
%             NRMSE_space3_time4_interp_time(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_time);
%         end 
%     end 
% end 
% close(nc1);close(nc21);close(nc22);
% 
% %
% NRMSE_space3_time4_interp_space_mean=mean(NRMSE_space3_time4_interp_space,2);
% NRMSE_space3_time4_interp_space_std=std(NRMSE_space3_time4_interp_space,0,2);
% 
% NRMSE_space3_time4_interp_time_mean=mean(NRMSE_space3_time4_interp_time,2);
% NRMSE_space3_time4_interp_time_std=std(NRMSE_space3_time4_interp_time,0,2);
% 
% NRMSE_space3_time4_NLmeans1_mean=mean(NRMSE_space3_time4_NLmeans1,2);
% NRMSE_space3_time4_NLmeans1_std=std(NRMSE_space3_time4_NLmeans1,0,2);
% 
% NRMSE_space3_time4_NLmeans2_mean=mean(NRMSE_space3_time4_NLmeans2,2);
% NRMSE_space3_time4_NLmeans2_std=std(NRMSE_space3_time4_NLmeans2,0,2);
% 
% 
% 
% %% SPACE 3 TIME 6
% space_spacing=3;
% time_spacing=6;
% 
% filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
% filename_NLM1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing6_sim12_acc12_neighbor5_tau0100.nc';
% filename_NLM2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing3_tspacing6_sim12_acc12_neighbor5_tau0100.nc';
% 
% HTHS_idy=1:Ny; % row indices of HTHS in space 
% HTHS_idz=1:Nz; % column indices of HTHS in space
% HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
% HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space
% LTHS_idt=1:time_spacing:Nx;
% 
% sborder=4*space_spacing;
% tborder=4*time_spacing;
% tleft=tborder; tright = tborder-(Nx-LTHS_idt(end));
% HTHS_idt_enlarged = 1:Nx+tleft+tright;
% LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);
% 
% nc1=netcdf(filename_ref,'r');
% nc21=netcdf(filename_NLM1,'r');
% nc22=netcdf(filename_NLM2,'r');
% 
% NRMSE_space3_time6_interp_space=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_space3_time6_interp_time=NRMSE_space3_time6_interp_space;
% NRMSE_space3_time6_NLmeans1=NRMSE_space3_time6_interp_space;
% NRMSE_space3_time6_NLmeans2=NRMSE_space3_time6_interp_space;
% for t=1:Nt
%     t
%     PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
%     PIV_sampled=cat(3, PIV_sampled(:,:,Nx-tborder+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:tright));    
%     LTHS_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
%     LTHS_interp=permute(LTHS_interp(tborder+1:tborder+Nx,:,:),[2,3,1]);
% 
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:time_spacing-1           
%             t_current = t_PIV_prev + pos_t;
%             x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
%             x_interp_time = squeeze(LTHS_interp(:,:,t_current));
% 
%             x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_space3_time6_NLmeans1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%  
%             x_rec = x_interp_space+nc22{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_space3_time6_NLmeans2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%             
%             NRMSE_space3_time6_interp_space(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_space);
%             NRMSE_space3_time6_interp_time(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_time);
%         end 
%     end 
% end 
% close(nc1);close(nc21);close(nc22);
% 
% %
% NRMSE_space3_time6_interp_space_mean=mean(NRMSE_space3_time6_interp_space,2);
% NRMSE_space3_time6_interp_space_std=std(NRMSE_space3_time6_interp_space,0,2);
% 
% NRMSE_space3_time6_interp_time_mean=mean(NRMSE_space3_time6_interp_time,2);
% NRMSE_space3_time6_interp_time_std=std(NRMSE_space3_time6_interp_time,0,2);
% 
% NRMSE_space3_time6_NLmeans1_mean=mean(NRMSE_space3_time6_NLmeans1,2);
% NRMSE_space3_time6_NLmeans1_std=std(NRMSE_space3_time6_NLmeans1,0,2);
% 
% NRMSE_space3_time6_NLmeans2_mean=mean(NRMSE_space3_time6_NLmeans2,2);
% NRMSE_space3_time6_NLmeans2_std=std(NRMSE_space3_time6_NLmeans2,0,2);
% 
% 
% %% SPACE 3 TIME 8
% space_spacing=3;
% time_spacing=8;
% 
% filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
% filename_NLM1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing8_sim12_acc12_neighbor5_tau0100.nc';
% filename_NLM2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing3_tspacing8_sim12_acc12_neighbor5_tau0100.nc';
% 
% HTHS_idy=1:Ny; % row indices of HTHS in space 
% HTHS_idz=1:Nz; % column indices of HTHS in space
% HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
% HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space
% LTHS_idt=1:time_spacing:Nx;
% 
% sborder=4*space_spacing;
% tborder=4*time_spacing;
% tleft=tborder; tright = tborder-(Nx-LTHS_idt(end));
% HTHS_idt_enlarged = 1:Nx+tleft+tright;
% LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);
% 
% nc1=netcdf(filename_ref,'r');
% nc21=netcdf(filename_NLM1,'r');
% nc22=netcdf(filename_NLM2,'r');
% 
% NRMSE_space3_time8_interp_space=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_space3_time8_interp_time=NRMSE_space3_time8_interp_space;
% NRMSE_space3_time8_NLmeans1=NRMSE_space3_time8_interp_space;
% NRMSE_space3_time8_NLmeans2=NRMSE_space3_time8_interp_space;
% for t=1:Nt
%     t
%     PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
%     PIV_sampled=cat(3, PIV_sampled(:,:,Nx-tborder+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:tright));    
%     LTHS_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
%     LTHS_interp=permute(LTHS_interp(tborder+1:tborder+Nx,:,:),[2,3,1]);
% 
%     for blockid=1:numel(LTHS_idt)-1
%         t_PIV_prev = LTHS_idt(blockid);
%         t_PIV_after = LTHS_idt(blockid+1);
%                
%         for pos_t=1:time_spacing-1           
%             t_current = t_PIV_prev + pos_t;
%             x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
%             [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
%             x_interp_time = squeeze(LTHS_interp(:,:,t_current));
% 
%             x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_space3_time8_NLmeans1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%  
%             x_rec = x_interp_space+nc22{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_space3_time8_NLmeans2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%             
%             NRMSE_space3_time8_interp_space(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_space);
%             NRMSE_space3_time8_interp_time(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_time);
%         end 
%     end 
% end 
% close(nc1);close(nc21);close(nc22);
% 
% %
% NRMSE_space3_time8_interp_space_mean=mean(NRMSE_space3_time8_interp_space,2);
% NRMSE_space3_time8_interp_space_std=std(NRMSE_space3_time8_interp_space,0,2);
% 
% NRMSE_space3_time8_interp_time_mean=mean(NRMSE_space3_time8_interp_time,2);
% NRMSE_space3_time8_interp_time_std=std(NRMSE_space3_time8_interp_time,0,2);
% 
% NRMSE_space3_time8_NLmeans1_mean=mean(NRMSE_space3_time8_NLmeans1,2);
% NRMSE_space3_time8_NLmeans1_std=std(NRMSE_space3_time8_NLmeans1,0,2);
% 
% NRMSE_space3_time8_NLmeans2_mean=mean(NRMSE_space3_time8_NLmeans2,2);
% NRMSE_space3_time8_NLmeans2_std=std(NRMSE_space3_time8_NLmeans2,0,2);
% 
% save NLmean_interps_vary_timespacing.mat ...
%      NRMSE_space3_time4_interp_space_mean NRMSE_space3_time4_interp_space_std ...
%      NRMSE_space3_time4_interp_time_mean NRMSE_space3_time4_interp_time_std ...
%      NRMSE_space3_time4_NLmeans1_mean NRMSE_space3_time4_NLmeans1_std ...
%      NRMSE_space3_time4_NLmeans2_mean NRMSE_space3_time4_NLmeans2_std ...
%      NRMSE_space3_time6_interp_space_mean NRMSE_space3_time6_interp_space_std ...
%      NRMSE_space3_time6_interp_time_mean NRMSE_space3_time6_interp_time_std ...
%      NRMSE_space3_time6_NLmeans1_mean NRMSE_space3_time6_NLmeans1_std ...
%      NRMSE_space3_time6_NLmeans2_mean NRMSE_space3_time6_NLmeans2_std ...
%      NRMSE_space3_time8_interp_space_mean NRMSE_space3_time8_interp_space_std ...
%      NRMSE_space3_time8_interp_time_mean NRMSE_space3_time8_interp_time_std ...
%      NRMSE_space3_time8_NLmeans1_mean NRMSE_space3_time8_NLmeans1_std ...
%      NRMSE_space3_time8_NLmeans2_mean NRMSE_space3_time8_NLmeans2_std;


%% NRMSE
load NLmean_interps_vary_timespacing.mat;

fsize=28;
fname='CMU Serif';

fig=figure();
set(fig, 'Position', [200 200 1000 800]);
set(fig,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.12 0.12 0.7 0.65]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

hold on;
plot(-4:4,repmat(mean(NRMSE_space3_time8_interp_space_mean),9),'g','LineWidth',3);
hLine11=plot(-2:2,[0;NRMSE_space3_time4_interp_time_mean;0],'m--s','LineWidth',1,'MarkerSize',17,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 1]); drawnow;
hLine12=plot(-2:2,[0;NRMSE_space3_time4_NLmeans1_mean;0],'r--s','LineWidth',1,'MarkerSize',17,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');drawnow;
hLine13=plot(-2:2,[0;NRMSE_space3_time4_NLmeans2_mean;0],'b--s','LineWidth',1,'MarkerSize',17,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');drawnow;

hLine21=plot(-3:3,[0;NRMSE_space3_time6_interp_time_mean;0],'m--o','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m');drawnow;
hLine22=plot(-3:3,[0;NRMSE_space3_time6_NLmeans1_mean;0],'r--o','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');drawnow;
hLine23=plot(-3:3,[0;NRMSE_space3_time6_NLmeans2_mean;0],'b--o','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');drawnow;

hLine31=plot(-4:4,[0;NRMSE_space3_time8_interp_time_mean;0],'m--d','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 1]); drawnow;
hLine32=plot(-4:4,[0;NRMSE_space3_time8_NLmeans1_mean;0],'r--d','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');drawnow;
hLine33=plot(-4:4,[0;NRMSE_space3_time8_NLmeans2_mean;0],'b--d','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');drawnow;


% hLine11=plot(-2:2,[0;NRMSE_space3_time4_interp_time_mean;0],'ms','LineWidth',1,'MarkerSize',17,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 1]); drawnow;
% hLine12=plot(-2:2,[0;NRMSE_space3_time4_NLmeans1_mean;0],'rs','LineWidth',1,'MarkerSize',17,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');drawnow;
% hLine13=plot(-2:2,[0;NRMSE_space3_time4_NLmeans2_mean;0],'bs','LineWidth',1,'MarkerSize',17,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');drawnow;
% 
% hLine21=plot(-3:3,[0;NRMSE_space3_time6_interp_time_mean;0],'mo','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','m');drawnow;
% hLine22=plot(-3:3,[0;NRMSE_space3_time6_NLmeans1_mean;0],'ro','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');drawnow;
% hLine23=plot(-3:3,[0;NRMSE_space3_time6_NLmeans2_mean;0],'bo','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');drawnow;
% 
% hLine31=plot(-4:4,[0;NRMSE_space3_time8_interp_time_mean;0],'md','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 1]); drawnow;
% hLine32=plot(-4:4,[0;NRMSE_space3_time8_NLmeans1_mean;0],'rd','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','r');drawnow;
% hLine33=plot(-4:4,[0;NRMSE_space3_time8_NLmeans2_mean;0],'bd','LineWidth',1,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor','b');drawnow;

% virtual plot for legend
h1=plot(-4,0,'s','LineWidth',2,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1,1,1]);
h2=plot(-4,0,'o','LineWidth',2,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1,1,1]);
h3=plot(-4,0,'d','LineWidth',2,'MarkerSize',15,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1,1,1]);
h4=plot(-4,0,'g-','LineWidth',2,'MarkerSize',15,'MarkerEdgeColor','g','MarkerFaceColor','g');
h5=plot(-4,0,'s','LineWidth',2,'MarkerSize',15,'MarkerEdgeColor','m','MarkerFaceColor','m');
h6=plot(-4,0,'s','LineWidth',2,'MarkerSize',15,'MarkerEdgeColor','r','MarkerFaceColor','r');
h7=plot(-4,0,'s','LineWidth',2,'MarkerSize',15,'MarkerEdgeColor','b','MarkerFaceColor','b');
hold off; 

box('on') ; 

xlim([-4,4])
% set(gca,'XTickLabel',{'1','2','3'})
xlabel('$\Delta t_c/\delta t$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

ylim([0 0.5]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.1:0.5);
set (gca, 'YTickLabel', {'0.0','0.1','0.2','0.3','0.4','0.5'},'FontSize',fsize)

leg={'$N/M=4$','$N/M=6$','$N/M=8$'};
legendflex([h1,h2,h3],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[-30 -120], ...
                       'nrow',1, ...
                       'fontsize',fsize-2,...
                       'Interpreter','Latex',...
                       'box','off');
leg={'$\mathbf{I}_s\mathbf{y}$','$\mathbf{I}_t\mathbf{x}$','$\mbox{Greedy propagation}$','$\mbox{Non-greedy propagation}$'};
legendflex([h4,h5,h6,h7],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[-30 -20], ...
                       'nrow',2, ...
                       'fontsize',fsize-2,...
                       'Interpreter','Latex',...
                       'box','off');

% % make transparent and plot
% % make transparent plot
% hMarkers = hLine11.MarkerHandle; hMarkers.FaceColorData = uint8(255*[1;0;1;0.4]);
% hMarkers = hLine12.MarkerHandle; hMarkers.FaceColorData = uint8(255*[1;0;0;0.4]);
% hMarkers = hLine13.MarkerHandle; hMarkers.FaceColorData = uint8(255*[0;0;1;0.4]);
% hMarkers = hLine21.MarkerHandle; hMarkers.FaceColorData = uint8(255*[1;0;1;0.4]);
% hMarkers = hLine22.MarkerHandle; hMarkers.FaceColorData = uint8(255*[1;0;0;0.4]);
% hMarkers = hLine23.MarkerHandle; hMarkers.FaceColorData = uint8(255*[0;0;1;0.4]);
% hMarkers = hLine31.MarkerHandle; hMarkers.FaceColorData = uint8(255*[1;0;1;0.4]);
% hMarkers = hLine32.MarkerHandle; hMarkers.FaceColorData = uint8(255*[1;0;0;0.4]);
% hMarkers = hLine33.MarkerHandle; hMarkers.FaceColorData = uint8(255*[0;0;1;0.4]);

export_fig('./figures/NLmean_interps_NRMSE_vary_timespacing', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()

