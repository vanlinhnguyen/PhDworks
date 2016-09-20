%PUTBACK_PATCHES to put back patches in the whole scene
%
%IN:
%   patches - patches of size [num_of_patches x dim]
%   grids - grids of starting points for each patch
%   sizes - size of each patch

%OUT:
%   X - output  field by overlapping the patches

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [X] = putback_patches(patches, grids, sizes)
    gridz=grids{1}; gridy=grids{2};
    Nz=sizes{1}; Ny=sizes{2};
    
    patchsize = sqrt(size(patches,1));
    X=zeros(Nz, Ny); cntMat=zeros(Nz, Ny);
          
    for ii = 1:length(gridz) 
        for jj = 1:length(gridy)
            zz = gridz(ii);
            yy = gridy(jj);

            patch_HR = reshape(patches(:,(ii-1)*length(gridy)+jj),[patchsize, patchsize]);
            X(yy:yy+patchsize-1, zz:zz+patchsize-1) = X(yy:yy+patchsize-1, zz:zz+patchsize-1) + patch_HR;
            cntMat(yy:yy+patchsize-1, zz:zz+patchsize-1) = cntMat(yy:yy+patchsize-1, zz:zz+patchsize-1) + 1;
        end
    end
    idx = (cntMat < 1); 
    cntMat(idx) = 1; X(idx) = NaN;
    X = X./cntMat;  
end
