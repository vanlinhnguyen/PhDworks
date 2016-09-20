%FILTER_2D  Return low-pass filter of input x_1D using a ideal Fourier
%filter, with scaling factor scale
%
%IN:
%   x_1D - an input 1D signal
%   scale - scale factor (cutoff freq k_c = k_{max}/scale)


%OUT:
%   x_1D_fil - fa filtered version of x_1D

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016
function x_1D_fil = filter_1D(x_1D, scale)

if iscolumn(x_1D)
    x_1D=x_1D';
end

N = numel(x_1D);
ks = (mod(1/2 + (0:(N-1))/N, 1) - 1/2); 
kc = max(ks)/scale;
LPF = (ks.^2 < kc^2); 

F = fft(x_1D);
x_1D_fil=ifft(LPF.*F);

end