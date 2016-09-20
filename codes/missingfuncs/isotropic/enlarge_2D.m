%ENLARGE_2D  enlarge a 2D signal using periodic boundary condition
%
%IN:
%   x_2D - input 2D field
%   [left, right, bottom, top] - number of point to extend

%OUT:
%   x_2D_expand - output 2D field

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function x_2D_expand = enlarge_2D(x_2D, left, right, bottom, top)
[nrow,ncol]=size(x_2D);

[gridz,gridy] = meshgrid(-left+1:ncol+right,-bottom+1:nrow+top);

gridy(gridy<1) = nrow+gridy(gridy<1); 
gridz(gridz<1) = ncol+gridz(gridz<1);
gridy(gridy>nrow) = gridy(gridy>nrow) - nrow; 
gridz(gridz>ncol) = gridz(gridz>ncol) - ncol;

x_2D_expand = x_2D(gridy + (gridz-1)*nrow);
end