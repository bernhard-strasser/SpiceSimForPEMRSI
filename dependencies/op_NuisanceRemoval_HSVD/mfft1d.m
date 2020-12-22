function xf   = mfft1d(x, N, dim, norm_on, shift_on)
% perform 1d fft with correct normalization and shift
% Chao Ma
% 

% by default, we do normalization
if ~exist('dim','var') || isempty(dim)
    dim       = 3;
end
if ~exist('norm_on','var') || isempty(norm_on)
	norm_on   = 1;
end
if ~exist('shift_on','var') || isempty(shift_on)
    shift_on  = 0;
end
if ~exist('N','var') || isempty(N)
    N         = size(x,dim);
end

% fft with correct shift
if shift_on
    xf        = fftshift(fft(ifftshift(x,dim),N,dim),dim);
else
    xf        = fftshift(fft(x,N,dim),dim);
end

% fft with normalization
if norm_on    == 1
	xf        = xf./sqrt(N);
end