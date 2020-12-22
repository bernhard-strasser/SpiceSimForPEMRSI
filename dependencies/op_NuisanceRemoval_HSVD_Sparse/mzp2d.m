function kzp = mzp2d(k,NewSize)
% 2D zero-padding
% Chao Ma
%

Ny                                     = size(k,1);
Nx                                     = size(k,2);
NyInterp                               = NewSize(1);
NxInterp                               = NewSize(2);

kzp                                    = zeros(NyInterp,NxInterp,size(k,3),size(k,4));
NyFloor                                = NyInterp/2 + 1 - Ny/2;
NyCeil                                 = NyInterp/2 + Ny/2;
NxFloor                                = NxInterp/2 + 1 - Nx/2;
NxCeil                                 = NxInterp/2 + Nx/2;

kzp(NyFloor:NyCeil,NxFloor:NxCeil,:,:) = k;