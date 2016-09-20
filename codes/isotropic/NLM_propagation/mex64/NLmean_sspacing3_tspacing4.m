% clear all; close all; 
clc;

%% INITIAL PARAMS
space_spacing=3;
time_spacing=4;

simpatch_haftsize = 3; % the half width of the patche
simpatch_fullsize = 2*simpatch_haftsize+1; % the full width.
accpatch_haftsize = simpatch_haftsize; 
accpatch_fullsize = 2*accpatch_haftsize + 1; 
neighbor_haftsize = 2; % the half width of the patch
neighbor_fullsize = 2*neighbor_haftsize+1; % the full width.
dim_simpatch=simpatch_fullsize*simpatch_fullsize;
dim_accpatch=accpatch_fullsize*accpatch_fullsize;
dim_neighbor=neighbor_fullsize*neighbor_fullsize;

filename_ref='/data/ISOTROPIC/data/data_downsampled4.nc';
nc = netcdf(filename_ref,'r');
Nt = nc('Nt').itsDimsize;
Nh = nc('Nz').itsDimsize;
close(nc)

dim_im=Nh*Nh;

tau = .1;

LTHS_idt=1:time_spacing:Nh;
HTHS_idy=1:Nh; % row indices of HTHS in space 
HTHS_idz=1:Nh; % column indices of HTHS in space
HTLS_idy= 1:space_spacing:Nh; % row indices of HTLS in space
HTLS_idz= 1:space_spacing:Nh; % column indices of HTLS in space

NRMSE=@(x_ref,x_est) sqrt(sum((x_est(:)-x_ref(:)).^2))/sqrt(sum(x_ref(:).^2));

%% Space interp
[gridpatches_y, gridpatches_z] = grid_gen_Ccodes(int16(Nh),int16(Nh),int16(simpatch_haftsize),int16(neighbor_haftsize));
% [gridpatches_y, gridpatches_z] = grid_gen_Ccodes(Nh,Nh,simpatch_haftsize,neighbor_haftsize);
% gridpatches_y=gridpatches_y(:);
% gridpatches_z=gridpatches_z(:);
     
acc_ids_1D = simpatch_haftsize+1-accpatch_haftsize:simpatch_haftsize+1+accpatch_haftsize;
acc_ids = reshape(0:dim_simpatch-1,[simpatch_fullsize,simpatch_fullsize]);
acc_ids = acc_ids(acc_ids_1D, acc_ids_1D);
acc_ids = acc_ids(:);

mex -O CFLAGS="\$CFLAGS -std=c99" NLmean.c
nc1=netcdf(filename_ref,'r');

for t=1:1
    fprintf('Estimating block: %.2d-th snapshot.\n', t);    
%     for blockid=1:numel(LTHS_idt)-1
    for blockid=1:1
        t_PIV_prev = LTHS_idt(blockid);
        t_PIV_after = LTHS_idt(blockid+1);
        
        x_PIV_prev= nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_PIV_prev);
        x_PIV_interp_prev = interp_border (x_PIV_prev, space_spacing);
        x_PIV_diff_prev = x_PIV_prev - x_PIV_interp_prev;
        
        x_PIV_after= nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_PIV_after);
        x_PIV_interp_after = interp_border (x_PIV_after, space_spacing);
        x_PIV_diff_after = x_PIV_after - x_PIV_interp_after;
        
        for pos_t=1:time_spacing-1
            start=tic();
            
            t_current = t_PIV_prev + pos_t;
            x_current = nc1{'velocity_x'}(t,HTHS_idz,HTHS_idy,t_current);
            x_current_interp = interp_border(x_current, space_spacing);
          
            [x_rec1,weights1] = NLmean(x_current_interp, x_PIV_interp_prev, x_PIV_diff_prev, gridpatches_y, gridpatches_z, simpatch_haftsize, neighbor_haftsize, acc_ids, tau);
%             [x_rec2,weights2] = NLmean(x_current_interp, x_PIV_interp_after, x_PIV_diff_after, gridpatches_y, gridpatches_z, simpatch_haftsize,accpatch_haftsize, neighbor_haftsize, tau);
%             x_rec = x_current_interp + (x_rec1 + x_rec2)./(weights1+weights2);
            x_rec = x_current_interp + x_rec1./weights1;
            
            NRMSE1 = NRMSE(x_current, x_current_interp);
            NRMSE2 = NRMSE(x_current, x_rec);
            fprintf(['\n Reconstruction in ',num2str(toc(start),'%.3f'),' seconds, with improvement of NRMSE is ', num2str(100*(1-NRMSE2/NRMSE1),'%.3f'), ' \n'])
        end
    end
end 
close(nc1); 

%%
% imagesc(weights1)   