clear all; close all; clc;

%% INITIAL PARAMETERS
N_HR=384; % Number of grid-point in streamwise direction

% cutoff
scale_factor_space=2;
k_max=N_HR/2;
k_max_real=2*k_max/3;
k_cutoff=k_max/2;
N_LR=2*k_cutoff;

%% Load 3D field
nc=netcdf('/data/ISOTROPIC/data/FIELD-020.nc','r'); 
velocity_x_resolved=nc{'velocity_x'}(:,:,:);
close(nc)

%% 3D spectra of full resolution
% fftn to frequency space
F = fftn(velocity_x_resolved);

k_3D_HR=zeros(N_HR,1);
E_HR=zeros(N_HR,1);
k_1D_HR=[0:N_HR/2 -N_HR/2+1:1:-1];

for m=1:N_HR
    m
    for n=1:N_HR
        for p=1:N_HR
            temp=sqrt(k_1D_HR(m)^2+k_1D_HR(n)^2+k_1D_HR(p)^2);
            k_3D_HR(round(temp)+1,1)=round(temp); % first wave number is zero
            E_HR(round(temp)+1,1) = E_HR(round(temp)+1,1) + abs(1/N_HR^3*F(m,n,p))^2;
        end
    end
end

%% 3D spectra of low resolution
% remove high frequency
ids=1:2*k_max;
ids([k_cutoff+1:k_max k_max+1:2*k_max-k_cutoff])=[];
F2=F(ids,ids,ids);


k_1D_LR=[0:N_LR/2 -N_LR/2+1:1:-1];
for m=1:N_LR
    m
    for n=1:N_LR
        for p=1:N_LR
            if round(sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2))==0
                F2(m,n,p)=0;
            end
        end
    end
end

% ifftn back to physical space (normalize by size)
velocity_x_filter=numel(F2)/numel(velocity_x_resolved)*ifftn(F2);

F_LR=fftn(velocity_x_filter);
k_3D_LR=zeros(N_HR,1);
E_LR=zeros(N_HR,1);

for m=1:N_LR
    m
    for n=1:N_LR
        for p=1:N_LR
            temp=sqrt(k_1D_LR(m)^2+k_1D_LR(n)^2+k_1D_LR(p)^2);
            k_3D_LR(round(temp)+1,1)=round(temp); % first wave number is zero
            E_LR(round(temp)+1,1) = E_LR(round(temp)+1,1) + abs(1/N_LR^3*F_LR(m,n,p))^2;
        end
    end
end


%% Plot
% -5/3 line
xline=[4,50];
yline=10^1*[exp((-5/3)*log(4)),exp((-5/3)*log(50))];

%
fsize=25;
fname='CMU Serif';
set(0,'DefaultAxesFontName', fname)
set(0,'DefaultAxesFontSize', fsize)
set(0,'DefaultTextFontname', fname)
set(0,'DefaultTextFontSize', fsize)

nc=netcdf('/data/ISOTROPIC/data/STAT.nc','r');
energy_spectra = nc{'energy_spectra'}(:,:);
close(nc)

fig1=figure(); 
set(gcf, 'Position', [400 100 1000 800]);
set(fig1, 'Color', 'w');

loglog(k_3D_HR(2:k_max_real),1.5*E_HR(2:k_max_real),'k-','LineWidth',2); 
hold on;
% loglog(k_3D_HR(2:128),energy_spectra(550,2:128),'r-'); % compare
loglog(k_3D_LR(2:k_max_real/scale_factor_space),1.5*E_LR(2:k_max_real/scale_factor_space),'r-','LineWidth',2); 
plot([k_max_real/scale_factor_space,k_max_real/scale_factor_space] ,[1*10^-8,2.5*10^-5],'b--','LineWidth',2)
plot([k_max_real/(2*scale_factor_space),k_max_real/(2*scale_factor_space)] ,[1*10^-8,1.2*10^-3],'b--','LineWidth',2)
plot([k_max_real/(3*scale_factor_space),k_max_real/(3*scale_factor_space)] ,[1*10^-8,6*10^-3],'b--','LineWidth',2)
plot([k_max_real/(6*scale_factor_space),k_max_real/(6*scale_factor_space)] ,[1*10^-8,4*10^-2],'b--','LineWidth',2)
loglog(xline,yline,'k-','LineWidth',2)
hold off
xlim([1,200]);
ylim([6*1e-8,5]);

text(20,0.5,'$-5/3$','HorizontalAlignment','right','FontSize',fsize+2,'Interpreter','Latex')
text(95,1.2*1e-7,'$\frac{k_{max}}{2}$','HorizontalAlignment','right','FontSize',fsize+2,'Interpreter','Latex')
text(48,1.2*1e-7,'$\frac{k_{max}}{4}$','HorizontalAlignment','right','FontSize',fsize+2,'Interpreter','Latex')
text(32,1.2*1e-7,'$\frac{k_{max}}{6}$','HorizontalAlignment','right','FontSize',fsize+2,'Interpreter','Latex')
text(16,1.2*1e-7,'$\frac{k_{max}}{12}$','HorizontalAlignment','right','FontSize',fsize+2,'Interpreter','Latex')

grid on;
xlabel('k'); ylabel('E(k)');

filename=strcat('./figures/spectra3d_ref');
export_fig(filename,'-eps','-q101','-a4','-nocrop');
close();
