%ENLARGE_1D  enlarge a 1D signal using periodic boundary condition
%
%IN:
%   x_1D - input signal
%   left - number of point to extend at the begining
%   right - number of point to extend at the end

%OUT:
%   x_1D_expand - output signal

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function x_1D_expand = enlarge_1D(x_1D, left, right)
%EXPAND_1D expand a 1D signal by [left,right]

N=numel(x_1D);

gridx = -left+1:N+right;

gridx(gridx<1) = N+gridx(gridx<1); 
gridx(gridx>N) = gridx(gridx>N) - N; 

x_1D_expand = x_1D(gridx);

end