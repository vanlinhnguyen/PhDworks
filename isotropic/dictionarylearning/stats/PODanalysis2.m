clear all; close all; clc;

%% INITIAL PARAMS
Nh=96; % Number of grid-point in HR field
Nt=37;
space_spacing=4;
time_spacing=6;

Nl = Nh/space_spacing;

size_l=4; 
patchsize_l = size_l; 
patchsize_h = patchsize_l*space_spacing;
dim_h=patchsize_h^2;

overlap=3; % at LR
num_patch=113664; % total number of patches

params.lambda=0.1;
params.K=2*dim_h;  


%% SR
Nx=numel(1:time_spacing:Nh);
Xall=zeros(Nh*Nh,Nt*Nx);
filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
count=1;
for t=1:Nt
    for i=1:time_spacing:Nh
        Xall(:,count)=reshape(nc{'velocity_x'}(t,:,:,i), Nh*Nh,1);
        count=count+1;
    end
end
close(nc);

%%
mu = mean(Xall,2);
Xall = bsxfun(@minus, Xall, mu); 
% Xall = bsxfun(@rdivide, Xall, sigma);
%[D, A, lambda] = POD(Xall);

%% POD
[K, A, ~, D, lambda] = PODGenerateFunc(Xall');
lambda(length(lambda)) = 0;                                                 % last eigenvalue should be zero                                               
rel_ener = lambda/sum(lambda);                                              % relative energy associated with mode m
acc_ener=cumsum(rel_ener);

%%
fsize=20;
fname='CMU Serif';

h=figure;
set(h, 'Position', [200 200 800 600]);
set(h,'color','w')

% Change default axes fonts.
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)

% Change default text fonts.
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

% define how the figure inside the plot appear on the paper
set(gcf,'Units','normal');
set(gca,'Position',[0.15 0.2 0.725 0.75]); % [x_leftlowcorner y_leftlowcorner width height]
set(gcf,'Units','pixels');

% plot here
plot((1:K+1)./(K+1),acc_ener,'r-','LineWidth',3);

xlabel('Sparsity')
ylabel('$\displaystyle \frac{\sum_{i=1}^{M} \left(\lambda_i\right)}{\sum_{i=1}^{D_h} \left(\lambda_i\right)}$','Interpreter','latex')

% ylabel('Accumulative variance')
ylim([0 1]);
% xlim([0 600]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.2:1);
set (gca, 'YTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.2:1);
set (gca, 'XTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)

box on
grid on

export_fig('./figures/POD_accumulative_energy','-eps');
% close()
