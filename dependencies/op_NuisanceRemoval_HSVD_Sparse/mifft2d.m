function x = mifft2d(k,norm_on)
	% perform inverse 2d fft with correct normalization and shift
	% Chao Ma
	%

temp = size(k);
Ny = temp(1);
Nx = temp(2);

% by default, we do normalization
if nargin<2
	norm_on = 1;
end

% ifft with correct shift
x = ifftshift(ifft(fftshift(...
    ifftshift(ifft(fftshift(...
    k,...
    1),[],1),1),...
    2),[],2),2); 

% ifft with normalization
if norm_on == 1
	x = x.*sqrt(Ny).*sqrt(Nx);
end