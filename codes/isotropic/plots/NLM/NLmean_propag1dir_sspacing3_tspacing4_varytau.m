clear all; close all; clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
filename_NLM1='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0050.nc';
filename_NLM2='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0075.nc';
filename_NLM3='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0100.nc';
filename_NLM4='/data/ISOTROPIC/NLM/interpdiff/server/NLmean_propag1dir_sspacing3_tspacing4_sim12_acc12_neighbor5_tau0150.nc';
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

%% Total error
sborder=4*space_spacing;
sleft=sborder; sright = sborder-(Ny-HTLS_idy(end));
sbottom=sborder; stop = sborder-(Nz-HTLS_idz(end));

tborder=4*time_spacing;
tleft=tborder; tright = tborder-(Nx-LTHS_idt(end));

HTHS_idt_enlarged = 1:Nx+tleft+tright;
LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);

nc1  = netcdf(filename_ref,'r');
nc21 = netcdf(filename_NLM1,'r');
nc22 = netcdf(filename_NLM2,'r');
nc23 = netcdf(filename_NLM3,'r');
nc24 = netcdf(filename_NLM4,'r');

NRMSE_NLmeans1=zeros(time_spacing-1,Nt*(numel(LTHS_idt)-1));
NRMSE_NLmeans2=NRMSE_NLmeans1;
NRMSE_NLmeans3=NRMSE_NLmeans1;
NRMSE_NLmeans4=NRMSE_NLmeans1;

for t=1:Nt
    t
    for blockid=1:numel(LTHS_idt)-1
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
               
        for pos_t=1:time_spacing-1           
            t_current = t_PIV_prev + pos_t;
            x_ref = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
            [x_interp_space, x_diff] = interp_border(x_ref, space_spacing, sborder, 1);

%             x_rec = x_interp_space+nc2{'x_rec_all'}(t,t_current,HTHS_idz,HTHS_idy);
            x_rec = x_interp_space+nc21{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
            NRMSE_NLmeans1(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);

            x_rec = x_interp_space+nc22{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
            NRMSE_NLmeans2(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);

            x_rec = x_interp_space+nc23{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
            NRMSE_NLmeans3(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);

            x_rec = x_interp_space+nc24{'x_HR_NLM_smallscales_all'}(t,t_current,HTHS_idz,HTHS_idy);
            NRMSE_NLmeans4(pos_t,(t-1)*(numel(LTHS_idt)-1)+blockid) = NRMSE(x_ref, x_rec);            
        end 
    end 
end 
close(nc1);close(nc21); close(nc22); close(nc23); close(nc24);

%%
NRMSE_NLmeans1_mean=mean(NRMSE_NLmeans1,2);
NRMSE_NLmeans1_std=std(NRMSE_NLmeans1,0,2);

NRMSE_NLmeans2_mean=mean(NRMSE_NLmeans2,2);
NRMSE_NLmeans2_std=std(NRMSE_NLmeans2,0,2);

NRMSE_NLmeans3_mean=mean(NRMSE_NLmeans3,2);
NRMSE_NLmeans3_std=std(NRMSE_NLmeans3,0,2);

NRMSE_NLmeans4_mean=mean(NRMSE_NLmeans4,2);
NRMSE_NLmeans4_std=std(NRMSE_NLmeans4,0,2);

%%
fsize=32;
fname='CMU Serif';

%%
y = [NRMSE_NLmeans1_mean NRMSE_NLmeans2_mean NRMSE_NLmeans3_mean NRMSE_NLmeans4_mean];
errY = [NRMSE_NLmeans1_std NRMSE_NLmeans2_std NRMSE_NLmeans3_std NRMSE_NLmeans4_std];

h=figure();
set(h, 'Position', [200 200 1000 800]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.2 0.2 0.75 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

h = barwitherr(errY, y);% Plot with errorbars

legend({'$q=01$','q=03','q=13','q=25'},'Intepreter','Latex','Location','northwest','FontSize',fsize-6)
legend('boxoff')
box on ; 

set(h(1), 'FaceColor',[0.8 0.8 0.8]) 
set(h(2), 'FaceColor',[0.7 0.7 0.7]) 
set(h(3), 'FaceColor',[0.6 0.6 0.6]) 
set(h(4), 'FaceColor',[0.5 0.5 0.5]) 

set(gca,'XTickLabel',{'1','2','3'})
xlabel('$t/\delta t$','Interpreter','latex')
ylabel('$\bar{\epsilon}$', 'interpreter', 'latex')

ylim([0 0.2]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.1:0.2);
set (gca, 'YTickLabel', {'0.0','0.1','0.2'},'FontSize',fsize)
export_fig('./figures/NLmean_propag1dir_spacing3_tspacing4_varytau_NRMSE', '-a1','-q101','-eps','-painters'); % leave '-painters' can cause quality problems
close()

