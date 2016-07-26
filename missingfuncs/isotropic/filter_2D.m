%FILTER_2D  Return low-pass filter of input x_2D using a ideal Fourier
%filter, with scaling factor scale_y and scale_z in spanwise and vertical
%direction
%
%IN:
%   x_2D - an input 2D field
%   scale_y - scale factor in vertical direction
%   scale_z - scale factor in spanwise


%OUT:
%   x_2D_fil - fa filtered version of x_2D

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016


function x_2D_fil = filter_2D(x_2D, scale_y, scale_z)

[nrows,ncols] = size(x_2D);

kys = (mod(1/2 + (0:(nrows-1))/nrows, 1) - 1/2); 
kzs = (mod(1/2 + (0:(ncols-1))/ncols, 1) - 1/2); 
kcz = max(kzs)/scale_z;
kcy = max(kys)/scale_y;
kc=sqrt(kcz^2+kcy^2);

[KZS,KYS] = meshgrid(kzs,kys); 

LPF = (KZS.*KZS + KYS.*KYS < kc^2); 

F = fftn(x_2D);
x_2D_fil=ifftn(LPF.*F);

end