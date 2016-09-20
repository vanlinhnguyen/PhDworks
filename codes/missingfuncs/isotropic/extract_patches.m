%EXTRACT_PATCHES to extract patches from a 2D image
%
%IN:
%   X_2D - input 2D field
%   grids - grids of starting points for each patch
%   patchsize - size of each patch

%OUT:
%   patches_all - output  patches of size [num_of_patches x patchsize]

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function [patches_all] = extract_patches(X_2D, grids, patchsize)
    gridz=grids{1}; gridy=grids{2};
    dims=size(X_2D);
    if length(dims)==3 % features, 3rd dim is the number of filters
        noFils=dims(3);
        patches_all=zeros(noFils*patchsize*patchsize,length(gridz)*length(gridy));
        for ii = 1:length(gridz)
            for jj = 1:length(gridy)
                zz = gridz(ii);
                yy = gridy(jj);
                patch = zeros(noFils*patchsize*patchsize,1);
                for kk=1:noFils
                    patch((kk-1)*patchsize*patchsize+1:kk*patchsize*patchsize) = X_2D(yy:yy+patchsize-1, zz:zz+patchsize-1,kk);
                end              
                patches_all(:,(ii-1)*length(gridy)+jj)=patch;
            end
        end                    
    else
        % collect patches
        patches_all=zeros(patchsize^2,length(gridz)*length(gridy));
        for ii = 1:length(gridz)
            for jj = 1:length(gridy)
                zz = gridz(ii);
                yy = gridy(jj);
                patch = X_2D(yy:yy+patchsize-1, zz:zz+patchsize-1);
                patches_all(:,(ii-1)*length(gridy)+jj)=patch(:);
            end
        end        
    end
end
