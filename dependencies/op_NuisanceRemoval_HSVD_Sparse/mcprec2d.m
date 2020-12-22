function xCp = mcprec2d(x,Psi,Nos,downsample)
% conjugate phase reconstruction
% Chao Ma
%

if nargin < 3
	Nos        = 1;
end
if nargin < 4
	downsample = 1;
end
if length(Nos) == 1
	Nos        = [Nos,Nos];
end
Ny 	           = size(x,1);
Nx 	           = size(x,2);

xZp            = mzprec2d(x,[Ny*Nos(1),Nx*Nos(2)]);
xZpBc          = xZp.*conj(Psi);
if downsample
	xCp		   = mdsrec2d(xZpBc,[Ny,Nx]);
else
	xCp        = xZpBc;
end