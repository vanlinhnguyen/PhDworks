%ESTIMATE_SPECT_2D  estimate a 2D energy spectrum of a square 2D field
%
%IN:
%   x_2D - input 2D field

%OUT:
%   E_HR - output energy spectrum

% Copyright (C) Linh Van Nguyen (linh.van.nguyen@hotmail.com) 2016
function E_HR = estimate_spect_2D(x_2D)
%ESTIMATE_SPECT_2D estimate 2D energy spectra

N = size(x_2D,1);
k_max=N/2;
k_1D_HR=[0:k_max -k_max+1:1:-1];

k_2D_HR=zeros(N,1);
E_HR=zeros(N,1);

F = fftn(x_2D);
for i=1:N
    for j=1:N
        temp=sqrt(k_1D_HR(i)^2+k_1D_HR(j)^2);
        k_2D_HR(round(temp)+1,1)=round(temp); % first wave number is zero
        E_HR(round(temp)+1,1) = E_HR(round(temp)+1,1) + abs(1/N^2*F(i,j))^2;
    end
end