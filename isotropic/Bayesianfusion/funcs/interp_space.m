function interp_space(gridy_HS, gridz_HS, HTHS_idy, HTHS_idz, HTLS_idy,HTLS_idz,filename_ref, filename_interp_space)
%INTERP_SPACE interpolate HTHS velocity fields from HTLS measurements for
%each snapshot and save to filename_interp_space
% interp_space(HTHS_idy, HTHS_idz, HTLS_idy,HTLS_idz,filename_interp_space) 

fprintf('\n ********** START SPATIAL INTERPOLATION ************ \n');
gridy_LS=gridy_HS(HTLS_idy,HTLS_idz);
gridz_LS=gridz_HS(HTLS_idy,HTLS_idz);

nc1=netcdf(filename_ref,'r');
Nt = nc1('Nt').itsDimsize;
Ny = nc1('Ny').itsDimsize;
Nz = nc1('Nz').itsDimsize;

nc2 = netcdf(filename_interp_space,'clobber');
nc2('Nt')=0; 
nc2('Ny')=Ny;
nc2('Nz')=Nz;
nc2{'Uinterp'}=ncfloat('Nt','Ny','Nz');

for t=1:Nt
    fprintf('Interpolating %.4d-th snapshot.\n', t);
    asnap_LTHS=nc1{'Uall'}(t,HTHS_idy,HTHS_idz);
    asnap_LTLS=asnap_LTHS(HTLS_idy,HTLS_idz);
    nc2{'Uinterp'}(t,:,:)=interp2(gridz_LS,gridy_LS,asnap_LTLS,gridz_HS,gridy_HS,'spline');
end 
close(nc1); close(nc2);

fprintf('\n ********** COMPLETED SPATIAL INTERPOLATION ************ \n');
end

