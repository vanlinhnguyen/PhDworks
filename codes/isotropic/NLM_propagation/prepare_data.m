clear all; close all; 
clc;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nz = nc('Nz').itsDimsize;
Ny = nc('Ny').itsDimsize;
Nx = nc('Nx').itsDimsize;
close(nc)

HTHS_idy=1:Ny; % row indices of HTHS in space 
HTHS_idz=1:Nz; % column indices of HTHS in space
HTLS_idt=1:1:Nx; % indices of all HTHS snapshots in time

%% Ratio 3x3
% Params
ratios = [3,4,6];

for i = 1:numel(ratios)
    space_spacing=ratios(i);

    HTLS_idy= 1:space_spacing:Ny; % row indices of HTLS in space
    HTLS_idz= 1:space_spacing:Nz; % column indices of HTLS in space

    % Space interp
    nc1=netcdf(filename_ref,'r');
    border=6*space_spacing;

    left=border; right = border-(Ny-HTLS_idy(end));
    bottom=border; top = border-(Nz-HTLS_idz(end));

    filename_interp = strcat('/data/PhDworks/isotropic/NLM/Udiff_spacespacing',num2str(ratios(i),'%.1d'),'.nc');
    nc = netcdf(filename_interp,'clobber'); 
    nc('Nt')=0;
    nc('Nx')=Nx;
    nc('Ny')=Ny;
    nc('Nz')=Nz;
    nc{'Uinterp_all'}=ncfloat('Nt','Nx','Ny','Nz');
    nc{'Udiff_all'}=ncfloat('Nt','Nx','Ny','Nz');

    for t=1:Nt
        fprintf('Estimating block: %.2d-th snapshot.\n', t);    
        for i=1:Nx
            x_ref= nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,i);
            [x_interp, x_diff] = interp_border (x_ref, space_spacing, border, 1);

            nc{'Uinterp_all'}(t,i,:,:) = x_interp;
            nc{'Udiff_all'}(t,i,:,:) = x_diff;
        end
    end 
    close(nc1); close(nc); 
end
