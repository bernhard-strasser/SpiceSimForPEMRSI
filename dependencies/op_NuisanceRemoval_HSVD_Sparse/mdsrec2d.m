function xd = mdsrec2d(x,NewSize,norm_on)
% 2D reconstruction while truncating half of the k-space data
% Chao Ma
%

if ~exist('norm_on','var')
    norm_on = 1;
end
temp        = size(x);
NyHr        = temp(1);
NxHr        = temp(2);
NyLr        = NewSize(1);
NxLr        = NewSize(2);
Nos         = NyHr/NyLr*NxHr/NxLr;

% x to k
kHr         = mfft2d(x,0);
% truncation
kLr         = mtr2d(kHr,[NyLr,NxLr]);
% k to x
xd          = mifft2d(kLr,0);

if norm_on
    xd      = xd/Nos;
end