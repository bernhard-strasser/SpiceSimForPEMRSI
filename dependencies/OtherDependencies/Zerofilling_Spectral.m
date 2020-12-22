function OutArray = Zerofilling_Spectral(OutArray,Zerofill_To,PerformFFT_flag)
%
% EllipticalFilter_x_y Apply an elliptical filter to k-space data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [OutArray,mask] = EllipticalFilter_x_y(OutArray,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims               ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients        ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     PerformFFT_flag              ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
%
% Output:
% -     OutArray                     ...     The filtered/masked output array
% -     mask                         ...     The mask of the filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat_1_0

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations
if(nargin < 1)
	OutArray = 0;
	return
end
if(nargin < 2)
    return
end
if(~exist('PerformFFT_flag','var'))
   PerformFFT_flag = false; 
end 


% 0.2 Declarations


% 0.3 Definitions
    
size_OutArray = size(OutArray);
size_OutArray = [size_OutArray ones([1 numel(Zerofill_To)-numel(size_OutArray)])];
AppendZeros = round(Zerofill_To - size_OutArray);
%AppendZeros(AppendZeros < 0) = 0;
ApplyAlongDims = find(ne(AppendZeros,0));



 





%% 1. FFT to k-space

if(PerformFFT_flag)

    for filter_dim = ApplyAlongDims
        OutArray = fft(OutArray,[],filter_dim);
        OutArray = fftshift(OutArray,filter_dim);
    end  

end




%% 2. Compute Mask


for dummy_dim = ApplyAlongDims

	if(AppendZeros(dummy_dim) > 0)
		AppendZeros_AtEnd = zeros([size_OutArray(1:dummy_dim-1) floor(AppendZeros(dummy_dim)) size_OutArray(dummy_dim+1:end)]);
		OutArray = cat(dummy_dim,OutArray,AppendZeros_AtEnd);
	else
		AccessPart = {'1', num2str(Zerofill_To(dummy_dim))};
		AccessString = [repmat(':,',[1 dummy_dim-1]) AccessPart{1} ':' AccessPart{2} ',' repmat(':,',[1 numel(size_OutArray)-dummy_dim])];
		AccessString(end) = [];
		OutArray = eval(['OutArray(' AccessString ');']);
	end
    size_OutArray = size(OutArray);
    
end


%% 4. FFT to ImageSpace


if(PerformFFT_flag)
    
    for filter_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,filter_dim);
        OutArray = ifft(OutArray,[],filter_dim);
    end  
    
end




%% 5. Postparations







