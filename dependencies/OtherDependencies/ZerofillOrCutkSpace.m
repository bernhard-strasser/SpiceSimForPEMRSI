function OutArray = ZerofillOrCutkSpace(OutArray,Zerofill_To,PerformFFT_flag)
%
% ZerofillOrCutkSpace Zerofill or Cut k-Space Data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function either zerofills or cuts data in k-space. For zerofilling time data 
% (zerofilling only one-sided instead of on both sides like in k-space), use "Zerofilling_Spectral".
%
%
% OutArray = ZerofillOrCutkSpace(OutArray,Zerofill_To,PerformFFT_flag)
%
% Input: 
% -     OutArray                     ...    Input array to which the filter should be applied.
% -     Zerofill_To                  ...    Array to which data should be zerofilled or cut. E.g. size(OutArray) = [32 64 64 512], Zerofill_To = [32 128 128 512]. 
% -     PerformFFT_flag              ...    If it is 1, the image gets Fourier transformed to k-space before applying the filter, 
%                                           and transformed back to image domain afterwards
%
% Output:
% -     OutArray                     ...    The filtered/masked output array
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
        OutArray = ifftshift(OutArray,filter_dim);
        OutArray = ifft(OutArray,[],filter_dim);
        OutArray = fftshift(OutArray,filter_dim);
    end  

end




%% 2. Compute Mask


for dummy_dim = ApplyAlongDims

	if(AppendZeros(dummy_dim) > 0)
		AppendZeros_AtBeginning = zeros([size_OutArray(1:dummy_dim-1) ceil(AppendZeros(dummy_dim)/2) size_OutArray(dummy_dim+1:end)]);
		AppendZeros_AtEnd = zeros([size_OutArray(1:dummy_dim-1) floor(AppendZeros(dummy_dim)/2) size_OutArray(dummy_dim+1:end)]);
		OutArray = cat(dummy_dim,AppendZeros_AtBeginning,OutArray,AppendZeros_AtEnd);
	else
		center = floor(size(OutArray,dummy_dim)/2) + 1;
		AccessPart = {num2str(center - floor(Zerofill_To(dummy_dim)/2)), num2str(center + ceil(Zerofill_To(dummy_dim)/2) - 1)};
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
        OutArray = fft(OutArray,[],filter_dim);
        OutArray = fftshift(OutArray,filter_dim);
    end  
    
end




%% 5. Postparations







