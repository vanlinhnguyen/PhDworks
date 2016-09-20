clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_NLM1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
filename_NLM2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag2dirs_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nz = nc('Nz').itsDimsize;
Ny = nc('Ny').itsDimsize; 
Nx = nc('Nx').itsDimsize;
close(nc)

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


%% Total error
% nc1=netcdf(filename_ref,'r');
% nc21=netcdf(filename_NLM1,'r');
% nc22=netcdf(filename_NLM2,'r');
% 
% NRMSE_interp_space=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
% NRMSE_interp_time=NRMSE_interp_space;
% NRMSE_NLmeans1=NRMSE_interp_space;
% NRMSE_NLmeans2=NRMSE_interp_space;
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
%             NRMSE_NLmeans1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%  
%             x_rec = x_interp_space+nc22{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);
%             NRMSE_NLmeans2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);
%             
%             NRMSE_interp_space(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_space);
%             NRMSE_interp_time(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_interp_time);
%         end 
%     end 
% end 
% close(nc1);close(nc21);close(nc22);
% 
% %
% NRMSE_interp_space_mean=mean(NRMSE_interp_space,2);
% NRMSE_interp_space_std=std(NRMSE_interp_space,0,2);
% 
% NRMSE_interp_time_mean=mean(NRMSE_interp_time,2);
% NRMSE_interp_time_std=std(NRMSE_interp_time,0,2);
% 
% NRMSE_NLmeans1_mean=mean(NRMSE_NLmeans1,2);
% NRMSE_NLmeans1_std=std(NRMSE_NLmeans1,0,2);
% 
% NRMSE_NLmeans2_mean=mean(NRMSE_NLmeans2,2);
% NRMSE_NLmeans2_std=std(NRMSE_NLmeans2,0,2);




% % % fig=figure();
% % % bar([NRMSE_interp_space_mean,NRMSE_interp_time_mean,NRMSE_NLmeans_mean]);
% % % hold on
% % % errorbar([NRMSE_interp_space_mean,NRMSE_interp_time_mean,NRMSE_NLmeans_mean],[NRMSE_interp_space_std,NRMSE_interp_time_std,NRMSE_NLmeans_std],'r.');
% % 
% % % temp=(time_spacing-2)*10+3;
% % % fig=figure();
% % % hold on;
% % % errorbar([1:10:temp],NRMSE_interp_space_mean,NRMSE_interp_space_std,'g-s');
% % % errorbar(2:10:temp,NRMSE_interp_time_mean,NRMSE_interp_time_std,'m-o');
% % % errorbar(3:10:temp,NRMSE_NLmeans_mean,NRMSE_NLmeans_std,'r-d');
% % % hold off
% % % ylim([0,0.3])
% % 

%% NRMSE
% fsize=28;
% fname='CMU Serif';
% 
% y = [NRMSE_interp_space_mean NRMSE_interp_time_mean NRMSE_NLmeans1_mean NRMSE_NLmeans2_mean];         % random y values (3 groups of 4 parameters) 
% errY = [NRMSE_interp_space_std NRMSE_interp_time_std NRMSE_NLmeans1_std NRMSE_NLmeans2_std];
% 
% fig=figure();
% set(fig, 'Position', [200 200 1000 800]);
% set(fig,'color','w')
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% % define how the figure inside the plot appear on the paper
% set(gcf,'Units','normal');
% set(gca,'Position',[0.2 0.2 0.75 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
% set(gcf,'Units','pixels');
% 
% h = barwitherr(errY, y);% Plot with errorbars
% 
% % gridLegend(h,10, [{'Spatial interp'};{'Time interp'};{'Propagation'}],'location','north');
% 
% % legend({'Spatial interp','Time interp','Greedy propagation','Non-greedy propagation'},'Location','northwest','FontSize',fsize-6)
% % legend('boxoff')
%                    
% box('on') ; 
%  
% % set(h(1), 'FaceColor',[0.9 0.9 0.9],'EdgeColor','g','LineWidth',3) 
% % set(h(2), 'FaceColor',[0.8 0.8 0.8],'EdgeColor','m','LineWidth',3) 
% % set(h(3), 'FaceColor',[0.7 0.7 0.7],'EdgeColor','r','LineWidth',3) 
% % set(h(4), 'FaceColor',[0.6 0.6 0.6],'EdgeColor','r','LineWidth',3) 
% 
% set(h(1), 'FaceColor',[0.8 0.8 0.8]) 
% set(h(2), 'FaceColor',[0.7 0.7 0.7]) 
% set(h(3), 'FaceColor',[0.6 0.6 0.6]) 
% set(h(4), 'FaceColor',[0.5 0.5 0.5]) 
% 
% set(gca,'XTickLabel',{'1','2','3'})
% xlabel('$t/\delta t$','Interpreter','latex')
% ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')
% 
% ylim([0 0.25]);
% 
% set(gca, 'YTickMode','manual');
% set(gca, 'YTick', 0:0.1:0.3);
% set (gca, 'YTickLabel', {'0.0','0.1','0.2','0.3'},'FontSize',fsize)
% 
% leg={'$\mathbf{I}_s\mathbf{y}$','$\mathbf{I}_s\mathbf{x}$','Greedy propagation','Non-greedy propagation'};
% hleg=legendflex(h,leg, 'ref', gcf, ...
%                        'anchor', {'n','n'}, ...
%                        'buffer',[60 -60], ...
%                        'nrow',2, ...
%                        'fontsize',fsize-2,...
%                        'Interpreter','Latex',...
%                        'box','off');
% 
% export_fig('./figures/NLmean_interps_NRMSE', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
% close()

%% SPECT
k_max=Ny/2;
k_1D_HR=[0:k_max -k_max+1:1:-1];
[KZS,KYS] = meshgrid(k_1D_HR,k_1D_HR); 
k_cutoff = k_max/space_spacing;

t_off1=[1,3]; % close to PIV

% FFT
nc1=netcdf(filename_ref,'r');
nc21=netcdf(filename_NLM1,'r');
nc22=netcdf(filename_NLM2,'r');

E_ref1=zeros(Ny,1);
E_interp_space1=zeros(Ny,1);
E_interp_time1=zeros(Ny,1);
E_NLmean11=zeros(Ny,1);
E_NLmean12=zeros(Ny,1);

E_err_interp_space1=zeros(Ny,1);
E_err_interp_time1=zeros(Ny,1);
E_err_NLmean11=zeros(Ny,1);
E_err_NLmean12=zeros(Ny,1);

norm_factor = 1/(Nt*(numel(LTHS_idt)-1)*numel(t_off1));
for t=1:Nt
    t
    PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
    PIV_sampled=cat(3, PIV_sampled(:,:,Nx-tborder+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:tright));    
    LTHS_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
    LTHS_interp=permute(LTHS_interp(tborder+1:tborder+Nx,:,:),[2,3,1]);

    for blockid=1:numel(LTHS_idt)-1
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
               
        for pos_t=1:numel(t_off1)        
            t_current = t_PIV_prev + t_off1(pos_t);
            x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
            [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
            x_interp_time = squeeze(LTHS_interp(:,:,t_current));

            x_NLmean1 = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
            x_NLmean2 = x_interp_space+nc22{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);
            
            % Spect
            E = estimate_spect_2D(x_ref); 
            E_ref1 = E_ref1 + norm_factor*E;

            E = estimate_spect_2D(x_interp_space);
            E_interp_space1 = E_interp_space1 + norm_factor*E;

            E = estimate_spect_2D(x_interp_time);
            E_interp_time1 = E_interp_time1 + norm_factor*E;

            E = estimate_spect_2D(x_NLmean1);
            E_NLmean11 = E_NLmean11 + norm_factor*E;

            E = estimate_spect_2D(x_NLmean2);
            E_NLmean12 = E_NLmean12 + norm_factor*E;
            
            % Error
            E = estimate_spect_2D(x_ref - x_interp_space);
            E_err_interp_space1 = E_err_interp_space1 + norm_factor*E;

            E = estimate_spect_2D(x_ref - x_interp_time);
            E_err_interp_time1 = E_err_interp_time1 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_NLmean1);
            E_err_NLmean11 = E_err_NLmean11 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_NLmean2);
            E_err_NLmean12 = E_err_NLmean12 + norm_factor*E;              
        end 
    end 
end 
close(nc1);close(nc21);close(nc22);




t_off2=[2]; % close to PIV

% FFT
nc1=netcdf(filename_ref,'r');
nc21=netcdf(filename_NLM1,'r');
nc22=netcdf(filename_NLM2,'r');

E_ref2=zeros(Ny,1);
E_interp_space2=zeros(Ny,1);
E_interp_time2=zeros(Ny,1);
E_NLmean21=zeros(Ny,1);
E_NLmean22=zeros(Ny,1);

E_err_interp_space2=zeros(Ny,1);
E_err_interp_time2=zeros(Ny,1);
E_err_NLmean21=zeros(Ny,1);
E_err_NLmean22=zeros(Ny,1);

norm_factor = 1/(Nt*(numel(LTHS_idt)-1)*numel(t_off2));
for t=1:Nt
    t
    PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
    PIV_sampled=cat(3, PIV_sampled(:,:,Nx-tborder+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:tright));    
    LTHS_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
    LTHS_interp=permute(LTHS_interp(tborder+1:tborder+Nx,:,:),[2,3,1]);

    for blockid=1:numel(LTHS_idt)-1
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
               
        for pos_t=1:numel(t_off2)        
            t_current = t_PIV_prev + t_off2(pos_t);
            x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
            [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);
            x_interp_time = squeeze(LTHS_interp(:,:,t_current));

            x_NLmean1 = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
            x_NLmean2 = x_interp_space+nc22{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);
            
            % Spect
            E = estimate_spect_2D(x_ref); 
            E_ref2 = E_ref2 + norm_factor*E;

            E = estimate_spect_2D(x_interp_space);
            E_interp_space2 = E_interp_space2 + norm_factor*E;

            E = estimate_spect_2D(x_interp_time);
            E_interp_time2 = E_interp_time2 + norm_factor*E;

            E = estimate_spect_2D(x_NLmean1);
            E_NLmean21 = E_NLmean21 + norm_factor*E;

            E = estimate_spect_2D(x_NLmean2);
            E_NLmean22 = E_NLmean22 + norm_factor*E;
            
            % Error
            E = estimate_spect_2D(x_ref - x_interp_space);
            E_err_interp_space2 = E_err_interp_space2 + norm_factor*E;

            E = estimate_spect_2D(x_ref - x_interp_time);
            E_err_interp_time2 = E_err_interp_time2 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_NLmean1);
            E_err_NLmean21 = E_err_NLmean21 + norm_factor*E;  

            E = estimate_spect_2D(x_ref - x_NLmean2);
            E_err_NLmean22 = E_err_NLmean22 + norm_factor*E;              
        end 
    end 
end 
close(nc1);close(nc21);close(nc22);

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
ax(1)=subplot(2,2,1,'position',[0.14 0.55 0.325 0.35]); % top left

% -5/3 line
xline=[4,20];
yline=5*[exp((-5/3)*log(4)),exp((-5/3)*log(20))];

h11=loglog(3:k_max,E_ref1(3:k_max),'k-','LineWidth',2); 
hold on;
h12=loglog(3:k_max,E_interp_space1(3:k_max),'g-','LineWidth',2); 
h13=loglog(3:k_max,E_interp_time1(3:k_max),'m-','LineWidth',2);
h14=loglog(3:k_max,E_NLmean11(3:k_max),'r-','LineWidth',2); 
h15=loglog(3:k_max,E_NLmean12(3:k_max),'b-','LineWidth',2); 
loglog(xline,yline,'k-','LineWidth',2)
plot([k_cutoff,k_cutoff] ,[1e-8,1e-2],'k--','LineWidth',2)
hold off
xlim([1,70]);
ylim([3*1e-6,1]);
text(15,0.5,'-5/3','HorizontalAlignment','right','FontSize',fsize)
text(5,1e-5,'$t/\delta t = 1$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

ylabel('$E(k)$','interpreter','latex');

box on
set(gca, 'XTick', []);

% =======    
ax(3)=axes('position',[0.52 0.55 0.325 0.35]); % top right
% -5/3 line
xline=[4,20];
yline=5*[exp((-5/3)*log(4)),exp((-5/3)*log(20))];

h21=loglog(3:k_max,E_ref2(3:k_max),'k-','LineWidth',2); 
hold on;
h22=loglog(3:k_max,E_interp_space2(3:k_max),'g-','LineWidth',2); 
h23=loglog(3:k_max,E_interp_time2(3:k_max),'m-','LineWidth',2);
h24=loglog(3:k_max,E_NLmean21(3:k_max),'r-','LineWidth',2); 
h25=loglog(3:k_max,E_NLmean22(3:k_max),'b-','LineWidth',2); 
loglog(xline,yline,'k-','LineWidth',2)
plot([k_cutoff,k_cutoff] ,[1e-8,1e-2],'k--','LineWidth',2)
hold off
xlim([1,70]);
ylim([3*1e-6,1]);
text(15,0.5,'-5/3','HorizontalAlignment','right','FontSize',fsize)
text(5,1e-5,'$t/\delta t = 2$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

% leg=legend([h1 h2 h3],{'$Ns=\;5$','$Ns=11$','$Ns=17$'},'location','northwest');
% set(leg,'FontSize',fsize-6);
% set(leg,'Interpreter','Latex')
% legend('boxoff')

set(gca, 'YTick', []);
set(gca, 'XTick', []);

box on;

% =======    
ax(2)=axes('position',[0.14 0.15 0.325 0.35]); % bottom left

h31=loglog(3:k_max,E_err_interp_space1(3:k_max)./E_ref1(3:k_max),'g-','LineWidth',2); 
hold on;
h32=loglog(3:k_max,E_err_interp_time1(3:k_max)./E_ref1(3:k_max),'m-','LineWidth',2); 
h33=loglog(3:k_max,E_err_NLmean11(3:k_max)./E_ref1(3:k_max),'r-','LineWidth',2); 
h34=loglog(3:k_max,E_err_NLmean12(3:k_max)./E_ref1(3:k_max),'b-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([0,60]);
ylim([1e-4,2]);

text(5,3*1e-1,'$t/\delta t = 1$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

xlabel('$k$','Interpreter','Latex');
ylabel('$E^{\epsilon}(k)/E(k)$','Interpreter','Latex');

box on;

% =======    

ax(4)=axes('position',[0.52 0.15 0.325 0.35]);
h41=loglog(3:k_max,E_err_interp_space2(3:k_max)./E_ref2(3:k_max),'g-','LineWidth',2); 
hold on;
h42=loglog(3:k_max,E_err_interp_time2(3:k_max)./E_ref2(3:k_max),'m-','LineWidth',2); 
h43=loglog(3:k_max,E_err_NLmean21(3:k_max)./E_ref2(3:k_max),'r-','LineWidth',2); 
h44=loglog(3:k_max,E_err_NLmean22(3:k_max)./E_ref2(3:k_max),'b-','LineWidth',2); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([0,60]);
ylim([1e-4,2]);

text(5,3*1e-1,'$t/\delta t = 2$','Interpreter','Latex','HorizontalAlignment','right','FontSize',fsize+5)

xlabel('$k$','Interpreter','Latex');
set(gca, 'YTick', []);
box on;

leg={'Reference','$\mathbf{I}_s \mathbf{y}$','$\mathbf{I}_s\mathbf{x}$','Greedy propagation','Non-greedy propagation'};
hleg=legendflex([h11,h12,h13,h14,h15],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 -50], ...
                       'nrow',1, ...
                       'fontsize',fsize,...
                       'Interpreter','Latex',...
                       'box','off');
filename_fig_spect=strcat('./figures/NLmean_interps_spectra2d');
export_fig(filename_fig_spect,'-a1','-q101','-eps','-painters');
close()     
%%
% fsize=32;
% fname='CMU Serif';
% 
% h=figure();
% set(h, 'Position', [200 200 1000 800]);
% set(h,'color','w')
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% % define how the figure inside the plot appear on the paper
% set(gcf,'Units','normal');
% set(gca,'Position',[0.2 0.2 0.75 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
% set(gcf,'Units','pixels');
% 
% % -5/3 line
% xline=[4,20];
% yline=5*[exp((-5/3)*log(4)),exp((-5/3)*log(20))];
% 
% h1=loglog(3:k_max,E_ref(3:k_max),'k-','LineWidth',2); 
% hold on;
% h2=loglog(3:k_max,E_interp_space(3:k_max),'g-','LineWidth',2); 
% h3=loglog(3:k_max,E_interp_time(3:k_max),'m-','LineWidth',2);
% h4=loglog(3:k_max,E_NLmean1(3:k_max),'r-','LineWidth',2); 
% h5=loglog(3:k_max,E_NLmean2(3:k_max),'b-','LineWidth',2); 
% loglog(xline,yline,'k-','LineWidth',2)
% plot([k_cutoff,k_cutoff] ,[1e-8,1e-2],'k--','LineWidth',2)
% hold off
% xlim([1,70]);
% ylim([3*1e-6,1]);
% 
% leg=legend([h1 h2 h3 h4 h5],{'Reference','$\mathbf{I}_s \mathbf{y}$','$\mathbf{I}_s\mathbf{x}$','Greedy propagation','Non-greedy propagation'},'location','southwest');
% set(leg,'FontSize',fsize-2);
% set(leg,'Interpreter','Latex')
% legend('boxoff')
% 
% text(10,0.5,'-5/3','HorizontalAlignment','right','FontSize',fsize+2)
% 
% xlabel('k'); ylabel('E(k)');
% 
% filename_fig_spect=strcat('./figures/NLmean_interps_spectra2d_ttd2');
% export_fig(filename_fig_spect,'-a1','-q101','-eps','-painters');
% close()

%%
% h=figure();
% set(h, 'Position', [100 100 1000 800]);
% set(h,'color','w')
% 
% % Change default axes fonts.
% set(0,'DefaultAxesFontName', fname)
% set(0,'DefaultAxesFontSize', fsize)
% 
% % Change default text fonts.
% set(0,'DefaultTextFontname', fname)
% set(0,'DefaultTextFontSize', fsize)
% 
% % define how the figure inside the plot appear on the paper
% set(gcf,'Units','normal');
% set(gca,'Position',[0.2 0.2 0.75 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
% set(gcf,'Units','pixels');
% 
% h1=loglog(3:k_max,E_err_interp_space(3:k_max)./E_ref(3:k_max),'g-','LineWidth',2); 
% hold on;
% h2=loglog(3:k_max,E_err_interp_time(3:k_max)./E_ref(3:k_max),'m-','LineWidth',2); 
% h3=loglog(3:k_max,E_err_NLmean1(3:k_max)./E_ref(3:k_max),'r-','LineWidth',2); 
% h4=loglog(3:k_max,E_err_NLmean2(3:k_max)./E_ref(3:k_max),'b-','LineWidth',2); 
% plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
% plot([1,100] ,[1,1],'k--','LineWidth',2)
% hold off
% xlim([0,60]);
% ylim([1e-4,15]);
% 
% xlabel('k'); ylabel('E(k)');
% 
% leg=legend([h1 h2 h3 h4],{'$\mathbf{I}_s \mathbf{y}$','$\mathbf{I}_s\mathbf{x}$','Greedy propagation','Non-greedy propagation'},'location','northwest');
% set(leg,'FontSize',fsize-6);
% set(leg,'Interpreter','Latex')
% legend('boxoff')
% 
% leg={'$\mathbf{I}_s \mathbf{y}$','$\mathbf{I}_s\mathbf{x}$','Greedy propagation','Non-greedy propagation'};
% hleg=legendflex([h1,h2,h3,h4],leg, 'ref', gcf, ...
%                        'anchor', {'n','n'}, ...
%                        'buffer',[60 -60], ...
%                        'nrow',2, ...
%                        'fontsize',fsize-2,...
%                        'Interpreter','Latex',...
%                        'box','off');
% 
% filename_fig_spect=strcat('./figures/NLmean_interps_errspectra2d_normalized_ttd2');
% export_fig(filename_fig_spect,'-a1','-q101','-eps','-painters');
% close()

