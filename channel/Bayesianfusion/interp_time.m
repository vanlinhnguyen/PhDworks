function interp_time(LTHS_idt, filename_ref, filename_interp_time)
%INTERP_TIME interpolate HTHS velocity fields from HTLS measurements for
%each snapshot and save to filename_interp_space
% interp_time(filename_ref, filename_interp_time)

fprintf('\n ********** START TIME INTERPOLATION ************ \n');

nc1=netcdf(filename_ref,'r');
PIV_sampled=nc1{'Uall'}(LTHS_idt,:,:);

Nt = nc1('Nt').itsDimsize;
Ny = nc1('Ny').itsDimsize;
Nz = nc1('Nz').itsDimsize;
close(nc1); 

t_all=1:Nt;

PIV_interp=zeros(Nt,Ny,Nz);

% interpolating by range will be faster than point by point
% considering do all runs at once if memory allows
ranges=1:50:Ny;
ranges=[ranges,Ny];
for i=1:numel(ranges)-1
    yids=ranges(i):ranges(i+1)-1;
    start=tic();
    PIV_interp(:,yids,:)=interp1(LTHS_idt,PIV_sampled(:,yids,:),t_all,'spline');
    fprintf('Estimating row %.2d-th to %.2d-th in %.2f seconds.\n', yids(1), yids(end), toc(start));
end

fprintf('\n ********** START WRITTING TO FILE ************ \n');
nc2 = netcdf(filename_interp_time,'clobber');
nc2('Nt')=0;
nc2('Ny')=Ny;
nc2('Nz')=Nz;
nc2{'Uinterp'}=ncfloat('Nt','Ny','Nz');

for t=1:Nt
    nc2{'Uinterp'}(t,:,:)=squeeze(PIV_interp(t,:,:));
end
close(nc2);

fprintf('\n ********** COMPLETED TIME INTERPOLATION ************ \n');

end