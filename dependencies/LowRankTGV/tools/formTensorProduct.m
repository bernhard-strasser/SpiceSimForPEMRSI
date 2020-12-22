function Xnew = formTensorProduct(X, Y,Spatial_dim)
% formTensorProduct - 
% Forms the mode-K tensor product between array X and vector/matrix Y,
% where K is the last dimension of X

nDimsX = ndims(X);
sizeX  = size(X);

% Handle degenerate cases for cases when there is only a single slice and /
% or a single component
if (nDimsX < (Spatial_dim+1) )
    sizeX = [sizeX, ones(1, (Spatial_dim+1) - nDimsX)];
end

Xnew = reshape(reshape(X, [], size(Y, 2)) * Y', [sizeX(1:(end - 1)), ...
        size(Y, 1)]);

end




