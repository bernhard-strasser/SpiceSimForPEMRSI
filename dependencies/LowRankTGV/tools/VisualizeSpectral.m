function VisualizeSpectral( V ,S,name_fig)
warning('off','MATLAB:DELETE:FileNotFound')
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here

nDims=ndims(V);
SizeData=size(V);
s=sprintf('%s_SpectralComp.ps',name_fig);
if exist(s);delete(s);end
ImSiC=2048;ImSiR=2048;
timept=1:size(V,1);

figs(1)=figure('visible', 'off');
dS=diag(S);
if(numel(dS)>64)
   dS = dS(1:64);
end
semilogy(1:numel(dS),abs(dS(:)),'-*k');
title('Singular Value Spectra')
print(figs(1), '-append', '-dpsc2', s);
close(figs(1))

for k=1:SizeData(end);
    figs(k+1)=figure('visible', 'off');
    plot(timept,real(fftshift(fft(V(:,k)'))),timept,imag(fftshift(fft(V(:,k)'))),timept,abs(fftshift(fft(V(:,k)'))),'r');
    print(figs(k+1), '-append', '-dpsc2', s);
    close(figs(k+1))
end


end

