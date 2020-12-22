function xzp = mzprec2d(x,NewSize,wHandle)
% 2D zero-padding reconstruction
% Chao Ma
%
if nargin < 3
	window_on = 0; 
else
	window_on = 1;
end

NyLr = size(x,1);
NxLr = size(x,2);
Nt = size(x,3);
NyHr = NewSize(1);
NxHr = NewSize(2);
Nos = NyHr/NyLr*NxHr/NxLr;

% x to k
kLr = mfft2d(x,0);

% windowing
if window_on
	wx = window(wHandle,NxLr);
	wy = window(wHandle,NyLr);
	w = wy(:)*reshape(wx,1,NxLr);
	w = repmat(w,[1,1,Nt]);
	kLr = kLr.*w;
end

% zero-padding
kHrzp = mzp2d(kLr,NewSize);

% k to x
xzp = mifft2d(kHrzp,0)*Nos;
