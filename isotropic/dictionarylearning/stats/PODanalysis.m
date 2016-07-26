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
C = (Xall*Xall')./(size(Xall,2)-1); %'# cov(X)

[D, lambda] = eig(C);
[lambda, order] = sort(diag(lambda), 'descend');       %# sort cols high to low
D = D(:,order);

% A = Xall*D(:,1:end);
A = D'*Xall;
%%
K = 100;
Z=D(:,1:K)'*Xall;
Xrec = D(:,1:K)*Z;

% fig1=figure(); imagesc(reshape(Xall(:,1),Nh,Nh)); caxis([-3,3]);
% fig2=figure(); imagesc(reshape(Xrec(:,1),Nh,Nh)); caxis([-3,3]);


%%
nc = netcdf(filename_ref,'r');
Xorg=reshape(nc{'velocity_x'}(1,:,:,11),Nh*Nh,1);
close(nc);
Xorg = bsxfun(@minus, Xorg, mu); 

%% 
K=round(0.1*Nh*Nh);
Anew = D(:,1:K)'*Xorg;
Anew = Anew./diag(D(:,1:K)'*D(:,1:K)); % remember to normalize by norm of each basis
Xreduced = D(:,1:K)*Anew;
% fig3=figure(); imagesc(reshape(Xorg,Nh,Nh)); caxis([-3,3]);
% fig4=figure(); imagesc(reshape(Xreduced,Nh,Nh)); caxis([-3,3]);
NRMSE = sqrt(sum((Xorg(:)-Xreduced(:)).^2))/sqrt(sum((Xorg(:)).^2))

%% 
thres_pca=1e-7:0.05:1;
Ks=round(thres_pca*Nh*Nh);
Ks(1)=1;

nc = netcdf(filename_ref,'r');
Xorg=reshape(nc{'velocity_x'}(1,:,:,11),Nh*Nh,1);
close(nc);
Xorg = bsxfun(@minus, Xorg, mu); 

NRMSE_PCA_all=zeros(numel(Ks),1);
for i=1:length(thres_pca)
    K=Ks(i);
    Anew = D(:,1:K)'*Xorg;
    Anew = Anew./(diag(D(:,1:K)'*D(:,1:K))); % remember to normalize by norm of each basis
    Xreduced = D(:,1:K)*Anew;
    %Xreduced = Xreduced + mu;
    NRMSE_PCA_all(i) = NRMSE_PCA_all(i) + sqrt(sum((Xorg(:)-Xreduced(:)).^2))/sqrt(sum((Xorg(:)).^2));
end

fig1 = figure(); imagesc(reshape(Xreduced(:,1),Nh,Nh)); caxis([-3,3]);
fig2 = figure(); imagesc(reshape(Xorg(:,1),Nh,Nh)); caxis([-3,3]);
fig3 = figure(); imagesc(reshape(mu(:,1),Nh,Nh)); caxis([-3,3]);

%%
acc_ener=(cumsum(lambda)./sum(lambda));
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
plot((1:920)./96^2,acc_ener(1:920),'r-','LineWidth',3);

xlabel('Sparsity')
ylabel('$\displaystyle \frac{\sum_{i=1}^{M} \left(\lambda_i\right)}{\sum_{i=1}^{D_h} \left(\lambda_i\right)}$','Interpreter','latex')

% ylabel('Accumulative variance')
ylim([0 1]);
% xlim([0 600]);

set(gca, 'YTickMode','manual');
set(gca, 'YTick', 0:0.2:1);
set (gca, 'YTickLabel', {'0.0','0.2','0.4','0.6','0.8','1.0'},'FontSize',fsize)

set(gca, 'XTickMode','manual');
set(gca, 'XTick', 0:0.02:0.1);
set (gca, 'XTickLabel', {'0.00','0.02','0.04','0.06','0.08','0.10'},'FontSize',fsize)

box on
grid on

export_fig('./figures/POD_accumulative_energy','-eps');
close();
