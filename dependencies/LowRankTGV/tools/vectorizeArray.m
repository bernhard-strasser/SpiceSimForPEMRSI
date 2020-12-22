function out = vectorizeArray(x)
% vectorizeArray - 
% Vectorizes an input array by stacking the input array columns. 
% This function produces the same effect as calling the MATLAB colon
% operator ':' on an array. For example, if X is an [M x N] array, then
% vectorizeArray(A) = A(:)

out = x(:);

end

