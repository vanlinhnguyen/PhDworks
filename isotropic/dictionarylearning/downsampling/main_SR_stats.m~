% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

close all; clc;

%% INITIAL PARAMS
space_spacing=4;
time_spacing=6;
Nt = 37;
Nh = 96;
Nl=Nh/space_spacing;
size_l=6; 

%% NRMSE

fprintf('\n ***********************************************************\n');
fprintf('\n ************************ NRMSE *************************** \n');
fprintf('\n ***********************************************************\n');

err = @(x1,x2) sqrt(sum((x1(:)-x2(:)).^2))/sqrt(sum(x2(:).^2)); % NRMSE
% err = @(x1,x2) sum((x1(:)-x2(:)).^2)/sum(x2(:).^2); % NMSE
[err_interp_mean,err_SR1_HRLR_mean,err_SR2_HRLRinterp_mean,err_SR3_feas_mean,...
    err_interp_std,err_SR1_HRLR_std,err_SR2_HRLRinterp_std,err_SR3_feas_std] = estimate_errors(err);

fprintf(['SR(feature) improves ',num2str(100-100*err_SR3_feas_mean/err_interp_mean,'%.3f'),' percents compared to interpolation\n']);
fprintf(['SR(interpolated LR) improves ',num2str(100-100*err_SR2_HRLRinterp_mean/err_interp_mean,'%.3f'),' percents compared to interpolation\n']);
fprintf(['SR(LR) improves ',num2str(100-100*err_SR1_HRLR_mean/err_interp_mean,'%.3f'),' percents compared to interpolation\n']);
%%
fsize=24;
fname='CMU Serif';

figure();

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [400 100 1200 400]);
set(gcf, 'Color', 'w');


y = [err_interp_mean err_SR1_HRLR_mean err_SR2_HRLRinterp_mean err_SR3_feas_mean];
errY = [err_interp_std err_SR1_HRLR_std err_SR2_HRLRinterp_std err_SR3_feas_std];

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.225 0.2 0.7 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

barh(1, y(1),0.5,'g');
hold on
barh(2, y(2),0.5,'b');
barh(3, y(3),0.5,'m');
barh(4, y(4),0.5,'r');
h1=herrorbar(y,1:4, errY,'k.');
plot([0.2754,0.2754],[0,5],'k--','linewidth',2) % level of error with subsampling
hold off

box on ; 
set(h1,'linewidth',1.5)

xlabel('$\bar{\epsilon}$', 'interpreter', 'latex') 

xlim([0 0.3]);
ylim([0.5 4.5]);

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.1:0.3);
set (gca, 'XTickLabel', {'0.0','0.1','0.2','0.3'},'FontSize',fsize)

set(gca, 'YTickMode','manual');
set (gca, 'YTickLabel', '');

text(-0.033,1,'Interp','Interpreter','latex','FontSize',fsize)
text(-0.025,2,'SR1','Interpreter','latex','FontSize',fsize)
text(-0.025,3,'SR2','Interpreter','latex','FontSize',fsize)
text(-0.025,4,'SR3','Interpreter','latex','FontSize',fsize)

filename_fig_spect=strcat('./figures/err_compare_all_spacespacing_',num2str(space_spacing,'%.2d'));
export_fig(filename_fig_spect,'-a1','-q101','-eps','-painters');
close() 



%% 2D SPECTRA
[E_ref,E_interp,E_SR1_HRLR,E_SR2_HRLRinterp, E_SR3_feas, ...
    E_err_interp, E_err_SR1_HRLR, E_err_SR2_HRLRinterp, E_err_SR3_feas] = estimate_2Dspect();

save estimate_2Dspect.mat E_ref E_interp E_SR1_HRLR E_SR2_HRLRinterp E_SR3_feas ...
    E_err_interp E_err_SR1_HRLR E_err_SR2_HRLRinterp E_err_SR3_feas;

% load estimate_2Dspect.mat

k_max = Nh/2;
k_1D_HR=[0:k_max -k_max+1:1:-1];
[KZS,KYS] = meshgrid(k_1D_HR,k_1D_HR); 
k_cutoff = k_max/space_spacing;

%% ENERGY SPECTRA
fsize=24;
fname='CMU Serif';

% -5/3 line
xline=[4,15];
yline=10^1*[exp((-5/3)*log(4)),exp((-5/3)*log(15))];

figure();

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

set(gcf, 'Position', [400 100 1000 800]);
set(gcf, 'Color', 'w');

h1=loglog(2:k_max,E_ref(2:k_max),'k-*','LineWidth',1.5); 
hold on;
h2=loglog(2:k_max,E_interp(2:k_max),'g-*','LineWidth',1.5); 
h3=loglog(2:k_max,E_SR3_feas(2:k_max),'r-*','LineWidth',1.5); 
h4=loglog(2:k_max,E_SR2_HRLRinterp(2:k_max),'m-*','LineWidth',1.5); 
h5=loglog(2:k_max,E_SR1_HRLR(2:k_max),'b-*','LineWidth',1.5); 

loglog(xline,yline,'k-','LineWidth',2)
plot([k_cutoff,k_cutoff] ,[1e-8,1e-2],'k--','LineWidth',2)
hold off
xlim([2,50]);
ylim([1e-8,5]);

text(10,0.6,'-5/3','HorizontalAlignment','right','FontSize',fsize+2)
text(14,2*1e-8,'$\mathbf{k_c}$','HorizontalAlignment','right','FontSize',fsize,'Interpreter','Latex')

xlabel('$k$','Interpreter','Latex'); ylabel('$E(k)$','Interpreter','Latex');

leg=legend([h1 h2 h5 h4 h3],{'Reference','$Interp$','$SR1$','$SR2$','$SR3$'},'location','southwest');
set(leg,'FontSize',fsize);
set(leg,'Interpreter','Latex')
legend boxoff
% grid on

filename_fig_spect=strcat('./figures/spectra2d_spacespacing_',num2str(space_spacing,'%.2d'));
export_fig(filename_fig_spect,'-eps','-q101','-a4');
close()

%% SPECTRA OF ERRORS
fsize=24;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 1400 600]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% ======= spectra of errors   
ax(1)=subplot(2,2,1,'position',[0.1 0.15 0.375 0.75]); % top left
loglog(2:k_max,E_err_interp(2:k_max),'g-*','LineWidth',1.5); 
hold on;
loglog(2:k_max,E_err_SR3_feas(2:k_max),'r-*','LineWidth',1.5); 
loglog(2:k_max,E_err_SR2_HRLRinterp(2:k_max),'m-*','LineWidth',1.5); 
loglog(2:k_max,E_err_SR1_HRLR(2:k_max),'b-*','LineWidth',1.5); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^-2],'k--','LineWidth',2)
hold off
text(15,8*1e-7,'$\mathbf{k_c}$','HorizontalAlignment','right','FontSize',fsize,'Interpreter','Latex')

xlim([2,50]);
ylim([5*1e-7 1e-2]);

xlabel('$k$','Interpreter','Latex'); ylabel('$E_\epsilon(k)$','Interpreter','Latex');
legend boxoff
% grid on 


% ======= spectra of errors (normalized by the energy spectrum of reference)
ax(2)=axes('position',[0.595 0.15 0.375 0.75]); % top right
h1=loglog(2:k_max,E_err_interp(2:k_max)./E_ref(2:k_max),'g-*','LineWidth',1.5); 
hold on;
h2=loglog(2:k_max,E_err_SR3_feas(2:k_max)./E_ref(2:k_max),'r-*','LineWidth',1.5); 
h3=loglog(2:k_max,E_err_SR2_HRLRinterp(2:k_max)./E_ref(2:k_max),'m-*','LineWidth',1.5); 
h4=loglog(2:k_max,E_err_SR1_HRLR(2:k_max)./E_ref(2:k_max),'b-*','LineWidth',1.5); 
plot([k_cutoff,k_cutoff] ,[10^-8,10^0],'k--','LineWidth',2)
plot([1,100] ,[1,1],'k--','LineWidth',2)
hold off
xlim([2,50]);
ylim([1e-6 2]);

text(15,2*1e-6,'$\mathbf{k_c}$','HorizontalAlignment','right','FontSize',fsize,'Interpreter','Latex')

xlabel('$k$','Interpreter','Latex'); ylabel('$E_\epsilon(k)/E(k)$','Interpreter','Latex');
% grid on


% plot in R2014
leg={'$Interp$','$SR1$','$SR2$','$SR3$'};
legendflex([h1,h4,h3,h2],leg, 'ref', gcf, ...
                       'anchor', {'n','n'}, ...
                       'buffer',[0 0], ...
                       'nrow',1, ...
                       'fontsize',fsize,...
                       'Interpreter','Latex',...
                       'box','off');
filename_fig_spect=strcat('./figures/errspectra2d_nonnormalized_normalized_timespacing_',num2str(time_spacing,'%.2d'),'_spacespacing_',num2str(space_spacing,'%.2d'));
export_fig(filename_fig_spect,'-a1','-q101','-eps','-painters');
close()     
