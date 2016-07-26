%RANDOM_INDICES  return starting point of a random patch
%
%IN:
%   Nl - size of LR field (Nl x Nl)
%   patchsize_l - size of LR patches
%   scale - the ratio to scale-up
%   num-patch - number of patches to extract from each field (maximum (Nl-patchsize_l+1)x(Nl-patchsize_l+1))


%OUT:
%   xrow_l,ycol_l - coordinate of starting point at LR
%   xrow_h,ycol_h - coordinate of starting point at HR

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [xrow_l,ycol_l,xrow_h,ycol_h] = random_indices(Nl,patchsize_l,scale, num_patch)
xl = randperm(Nl-patchsize_l); % taking care of bolder (convolution)
yl = randperm(Nl-patchsize_l); % taking care of bolder (convolution)
xl=xl(1:sqrt(num_patch));
yl=yl(1:sqrt(num_patch));

[Xl,Yl] = meshgrid(xl,yl);
 
xh=zeros(size(xl)); yh=xh; % starting point of HR patches
for i=1:numel(xl)
    % (a,b) in LR corresponding to ((a-1)*spacing+1:a*spacing,(b-1)*spacing+1:b*spacing) in HR
    xh(1,i) = (xl(i)-1)*scale+1;
    yh(1,i) = (yl(i)-1)*scale+1;
end
[Xh,Yh] = meshgrid(xh,yh);

xrow_l = Xl(:); ycol_l = Yl(:); 
xrow_h = Xh(:); ycol_h = Yh(:); 