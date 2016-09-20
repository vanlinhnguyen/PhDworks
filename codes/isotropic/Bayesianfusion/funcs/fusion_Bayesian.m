%FUSION_BAYESIAN fusion using MAP estimate
%
%IN:
%    space_spacing - subsampling in space
%    time_spacing - subsampling in time

%OUT:

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function fusion_Bayesian(space_spacing, time_spacing)

fprintf('\n ========================================================= \n'); 
fprintf('\n BAYESIAN FUSION: (SSPACING,TSPACING) = (%.1d,%.1d) \n', space_spacing,time_spacing); 
fprintf('\n ========================================================= \n'); 

%% Define file locations (to load or to save)
filename_ref='/data/PhDworks/isotropic/refdata_downsampled4.nc';
filename_fusion=strcat('/data/PhDworks/isotropic/Bayesianfusion/FusedData_BF_tspacing',num2str(time_spacing,'%.1d'),'_sspacing',num2str(space_spacing,'%.1d'),'.nc');

Nt = 37;
Nx = 96;
Ny = 96;
Nz = 96;

HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

LTHS_idt=1:time_spacing:Nx;

%% ============================ Interpolation step ========================
%  ========================================================================
%  ========================================================================
% ==== Space interp
start=tic();

left=4*space_spacing; right = 4*space_spacing-(Ny-HTLS_idy(end));
bottom=4*space_spacing; top = 4*space_spacing-(Nz-HTLS_idz(end));
[gridz_HS_enlarged, gridy_HS_enlarged] = meshgrid(1:Nz+left+right,1:Ny+bottom+top);
gridz_LS_enlarged=gridz_HS_enlarged(1:space_spacing:end,1:space_spacing:end);
gridy_LS_enlarged=gridy_HS_enlarged(1:space_spacing:end,1:space_spacing:end);

nc1=netcdf(filename_ref,'r');
IsY_all = zeros(Nt,Nz,Ny,Nx);
for t=1:Nt
    fprintf('Interpolating %.4d-th snapshot.\n', t);    
    for i=1:Nx
        % interp
        asnap_LTHS=nc1{'velocity_x'}(t,:,:,i);
        asnap_LTHS_enlarged=enlarge_2D(asnap_LTHS,left, right, bottom, top);
        
        asnap_LTLS=asnap_LTHS_enlarged(1:space_spacing:end,1:space_spacing:end);
        asnap_LTHS_interp=interp2(gridz_LS_enlarged, gridy_LS_enlarged, asnap_LTLS, gridz_HS_enlarged, gridy_HS_enlarged,'spline');
               
        IsY_all(t,:,:,i) = asnap_LTHS_interp(bottom+1:bottom+Ny,left+1:left+Nz);
    end
end 
close(nc1);

fprintf('Spatial interpolation in %.2f seconds.\n', toc(start)); 

% ==== Time interp
start=tic();

left=4*time_spacing; right = 4*time_spacing-(Nx-LTHS_idt(end));
HTHS_idt_enlarged = 1:Nx+left+right;
LTHS_idt_enlarged = 1:time_spacing:numel(HTHS_idt_enlarged);

nc1=netcdf(filename_ref,'r');

ItX_all = zeros(Nt,Nz,Ny,Nx);

for t=1:Nt
    fprintf('Interpolating %.4d-th snapshot.\n', t);    
    PIV_sampled=nc1{'velocity_x'}(t,:,:,:);
    PIV_sampled=cat(3, PIV_sampled(:,:,Nx-left+1:1:Nx), PIV_sampled, PIV_sampled(:,:,1:1:right));
    
    PIV_interp=interp1(LTHS_idt_enlarged, permute(PIV_sampled(:,:,LTHS_idt_enlarged),[3 1 2]),HTHS_idt_enlarged,'spline');
    ItX_all(t,:,:,:)=permute(PIV_interp(left+1:left+Nx,:,:),[2,3,1]);
end

close(nc1);

fprintf('Time interpolation in %.2f seconds.\n', toc(start)); 


%% ======================= Estimate fusion parameters =====================
%  ========================================================================
%  ========================================================================

% ==== Compute Cns 
start=tic();
nc1 = netcdf(filename_ref,'r');
hs=nc1{'velocity_x'}(:,:,:,LTHS_idt) - IsY_all(:,:,:,LTHS_idt);
close(nc1);
var_hs=hs.^2; clearvars hs;

S=squeeze(sum(var_hs,1));
S=squeeze(sum(S,3));
S=S./(Nt*numel(LTHS_idt));

% average
S_ave=zeros(Ny,Nz);
for i=1:space_spacing
    for j=1:space_spacing
        S_ave(i:space_spacing:Ny,j:space_spacing:Nz) = mean(mean(S(i:space_spacing:Ny,j:space_spacing:Nz)));
    end  
end

% ==== Compute Cnt
T_ave=zeros(time_spacing,1);

nc1 = netcdf(filename_ref,'r');
for i=1:time_spacing
    ht=nc1{'velocity_x'}(:,HTLS_idz,HTLS_idy,i:time_spacing:Nx) - ItX_all(:,HTLS_idz,HTLS_idy,i:time_spacing:Nx);
    T_ave(i)=mean(ht(:).^2);
end
close(nc1);

fprintf('Learning parameters in %.2f seconds.\n', toc(start)); 


%% =============================== Fusion =================================
%  ========================================================================
%  ========================================================================
start=tic();

nc2 = netcdf(filename_fusion,'clobber');
nc2('Nt')=0;
nc2('Nx')=Nx;
nc2('Ny')=Ny;
nc2('Nz')=Nz;
nc2{'Zhat_all'}=ncfloat('Nt','Nz','Ny','Nx');

for t=1:Nt
    fprintf('Fusion %.2d-th block.\n', t);
    Zhat=zeros(Nz,Ny,Nx);
    for i=1:Nx
        pos_t=rem(i-1,time_spacing)+1;

        % estimate weights
        if(pos_t==1)
            Wt=ones(Nz,Ny);
            Ws=zeros(Nz,Ny);
        else
            Wt=S_ave./(S_ave+T_ave(pos_t));
            Ws=T_ave(pos_t)./(S_ave+T_ave(pos_t));
        end

        % fusion
        IsY =  squeeze(IsY_all(t,:,:,i));
        ItX =  squeeze(ItX_all(t,:,:,i));

        Zhat(:,:,i) = Wt.*ItX + Ws.*IsY;
    end 
    nc2{'Zhat_all'}(t,:,:,:) = Zhat;
end
close(nc2);
fprintf('Fusion in %.2f seconds.\n', toc(start)); 