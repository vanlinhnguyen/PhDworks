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