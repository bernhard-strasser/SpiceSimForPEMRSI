function [ m, n ] = ind2sub2d( M, i )
%IND2SUB3D fast 2d version of ind2sub
%
% ---- INPUTS ----
% M     = size of 1st array dimension
% i     = column vector of array index values
%
% ---- OUTPUTS ----
% m,n   = column vectors of subscripts for 1st and 2nd dimensions
%-------------------------------------------------------------------------------
n = floor((i-1)/M) + 1;
% m = mod((i-1),M) + 1;
m = i - M*(n-1);

end

