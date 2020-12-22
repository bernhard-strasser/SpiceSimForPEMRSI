function k = mfft2d(x,norm_on)
% perform 2d fft with correct normalization and shift
% Chao Ma
% 

temp = size(x);
Ny = temp(1);
Nx = temp(2);

% by default, we do normalization
if nargin<2
	norm_on = 1;
end

% fft with correct shift
k = fftshift(fft(ifftshift(...
    fftshift(fft(ifftshift(...
    x,...
    1),[],1),1),...
    2),[],2),2);

% fft with normalization
if norm_on == 1
	k = k./sqrt(Ny)./sqrt(Nx);
end