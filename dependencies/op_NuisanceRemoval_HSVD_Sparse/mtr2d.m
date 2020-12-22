function ktr = mtr2d(k,NewSize,Shifted)
% perform 2d truncation
% Chao Ma
%

if nargin<3
	Shifted = 1; % truncation for shifted data
end
temp = size(k);
Ny = temp(1);
Nx = temp(2);
NyLr = NewSize(1);
NxLr = NewSize(2);

if Shifted
ktr = k((Ny/2+1-NyLr/2):(Ny/2+NyLr/2),...
        (Nx/2+1-NxLr/2):(Nx/2+NxLr/2),:,:);
else
	Nt = size(k,3);
    Nc = size(k,4);
	ktr = zeros(NyLr,NxLr,Nt,Nc);
	% top left corner
	ktr(1:NyLr/2,1:NxLr/2,:,:) = k(1:NyLr/2,1:NxLr/2,:,:);
	% top right corner
	ktr(1:NyLr/2,NxLr/2+1:NxLr,:,:) = k(1:NyLr/2,(Nx-NxLr/2+1):Nx,:,:);
	% bottom left corner
	ktr(NyLr/2+1:NyLr,1:NxLr/2,:,:) = k((Ny-NyLr/2+1):Ny,1:NxLr/2,:,:);
	% bottom right corner
	ktr(NyLr/2+1:NyLr,NxLr/2+1:NxLr,:,:) = k((Ny-NyLr/2+1):Ny,(Nx-NxLr/2+1):Nx,:,:);
end