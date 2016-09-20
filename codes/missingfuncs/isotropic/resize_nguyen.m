%RESIZE_NGUYEN see imresize (matlab)
%
%IN:
%   Xin - input 2D field
%   scale - ratio to scale Xin
%   kernel - interpolation kernel

%OUT:
%   Xout - scaled output field

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016

function Xout = resize_nguyen(Xin, scale,kernel)
%RESIZE_NGUYEN Resize a 2D matrix.
    switch(kernel)
        case 'bicubic'
            % Calculate interpolation weights and indices for each dimension.
            kernel_width=4; % size of kernel using bicubic
            antialiasing=1; % using antialiasing
            weights = cell(1,2);
            indices = cell(1,2);
            for k = 1:2
                [weights{k}, indices{k}] = contributions(size(Xin, k), ...
                    size(Xin, k)*scale, scale, @cubic, ...
                    kernel_width, antialiasing);
            end

            Xout = Xin;
            for dim = 1:2
                Xout = imresizemex(Xout, weights{dim}', indices{dim}', dim);
            end
        otherwise
            fprintf('bicubic only!\n');
    end


function [weights, indices] = contributions(in_length, out_length, ...
                                            scale, kernel, ...
                                            kernel_width, antialiasing)

    if (scale < 1) && (antialiasing)
        % Use a modified kernel to simultaneously interpolate and
        % antialias.
        h = @(x) scale * kernel(scale * x);
        kernel_width = kernel_width / scale;
    else
        % No antialiasing; use unmodified kernel.
        h = kernel;
    end

    % Output-space coordinates.
    x = (1:out_length)';

    % Input-space coordinates. Calculate the inverse mapping such that 0.5
    % in output space maps to 0.5 in input space, and 0.5+scale in output
    % space maps to 1.5 in input space.
    u = x/scale + 0.5 * (1 - 1/scale);

    % What is the left-most pixel that can be involved in the computation?
    left = floor(u - kernel_width/2);

    % What is the maximum number of pixels that can be involved in the
    % computation?  Note: it's OK to use an extra pixel here; if the
    % corresponding weights are all zero, it will be eliminated at the end
    % of this function.
    P = ceil(kernel_width) + 2;

    % The indices of the input pixels involved in computing the k-th output
    % pixel are in row k of the indices matrix.
    indices = bsxfun(@plus, left, 0:P-1);

    % The weights used to compute the k-th output pixel are in row k of the
    % weights matrix.
    weights = h(bsxfun(@minus, u, indices));

    % Normalize the weights matrix so that each row sums to 1.
    weights = bsxfun(@rdivide, weights, sum(weights, 2));

    % Clamp out-of-range indices; has the effect of replicating end-points.
    indices = min(max(1, indices), in_length);

    % If a column in weights is all zero, get rid of it.
    kill = find(~any(weights, 1));
    if ~isempty(kill)
        weights(:,kill) = [];
        indices(:,kill) = [];
    end
 

%=====================================================================
function f = cubic(x)
% See Keys, "Cubic Convolution Interpolation for Digital Image
% Processing," IEEE Transactions on Acoustics, Speech, and Signal
% Processing, Vol. ASSP-29, No. 6, December 1981, p. 1155.

    absx = abs(x);
    absx2 = absx.^2;
    absx3 = absx.^3;

    f = (1.5*absx3 - 2.5*absx2 + 1) .* (absx <= 1) + ...
                    (-0.5*absx3 + 2.5*absx2 - 4*absx + 2) .* ...
                    ((1 < absx) & (absx <= 2));
%---------------------------------------------------------------------