function [OutArray,ShiftMap,RefSpec] = FrequencyAlignment(InArray,PeakSearchSettingsOrShiftMap,ApplyAlongDim,ZerofillingFactor,RefSpecIn)
%
% FrequencyAlignment Align frequencies of csi spectra.
%
% This function was written by Bernhard Strasser, April 2015.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [OutArray,mask] = FrequencyAlignment(OutArray,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
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
size_InArray = size(InArray);

if(~exist('RefSpecIn','var'))
	RefSpecIn = [];
end
RefSpec = RefSpecIn;

% Handle mask
MaskWasProvided = true;
if(isstruct(InArray) && isfield(InArray,'mask'))
	mask = InArray.mask;
	InArray = InArray.csi;
	size_InArray = size(InArray);
else
	mask = ones(size_InArray(setdiff(1:numel(size_InArray),ApplyAlongDim)));
	MaskWasProvided = false;	
end

if(isfield(PeakSearchSettingsOrShiftMap,'PeakSearchPPM'))
	Settings = PeakSearchSettingsOrShiftMap;
	if(~isfield(Settings,'PolyfitRegion'))
		Settings.PolyfitRegion = [3.4 2.01];
	end
	Settings.PolyfitRegion = sort(Settings.PolyfitRegion,2,'descend');
	OnlyApplyShiftMap_flag = false;
else
	ShiftMap = PeakSearchSettingsOrShiftMap;
	OnlyApplyShiftMap_flag = true;
    RefSpecWasProvided = false;
end


% Define Reference Voxel if necessary
if(~OnlyApplyShiftMap_flag)
    RefSpecWasProvided = true;
    if(isempty(RefSpec))
        if(MaskWasProvided)
            % Define Reference Voxel as Center of Mass Voxel
            RefVox = regionprops(mask, 'centroid');
            RefVox = round(RefVox.Centroid);
        else
            RefVox = floor(size_InArray/2)+1;
            fprintf('\n\nWARNING: No mask input for FrequencyAlignment. Reference voxel will be set as (%d,%d,%d). Might be wrong!',RefVox(1),RefVox(2),RefVox(3))
        end

        if(numel(RefVox) < 3)
            RefVox(3) = 1;
        end

        % Set RefSpec
        RefSpec = InArray(RefVox(1),RefVox(2),RefVox(3),:);
        RefSpecWasProvided = false;
    end
end






% 0.2 Declarations


% 0.3 Definitions
    
size_SearchArray = size_InArray; size_SearchArray(ApplyAlongDim) = size_SearchArray(ApplyAlongDim)*ZerofillingFactor;
OutArray = InArray;
SearchArray = InArray;




%% 1. Zerofilling & FFT SearchArray

if(~OnlyApplyShiftMap_flag)

    SearchArray = Zerofilling_Spectral(SearchArray,size_SearchArray,0);
    SearchArray = fftshift(fft(SearchArray,[],ApplyAlongDim),ApplyAlongDim);

    RefSpec2 = Zerofilling_Spectral(RefSpec,[ones([1 numel(size_SearchArray)-1]) size_SearchArray(ApplyAlongDim)],0);
    RefSpec2 = fftshift(fft(RefSpec2,[],ApplyAlongDim),ApplyAlongDim);
end

%% 2. Shift RefSpec by a Certain Amount of Points.

% Instead of shifting the individual spectra of all voxels N times, computing the DotProduct and finding the maximum,
% do instead: Shift the RefSpec once N times, save that, and calculate the DotProducts for all shifts of the RefSpec.
% In the fist case, we would need to do NumbOfVox x N shifts, whereas in the latter we only need N shifts.
if(~OnlyApplyShiftMap_flag)
	
	CS_vec_zf = compute_chemshift_vector_1_1(Settings.LarmorFreq,Settings.Dwelltime/10^9,Settings.vecsize*ZerofillingFactor); 
	ShiftMap = zeros([size(OutArray,1) size(OutArray,2) size(OutArray,3)]);
	SearchForPeak_LeftPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM));
	SearchForPeak_RightPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM));
	SearchForPeak_Center_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM));

	ReferenceSpecMat_Spec = squeeze(RefSpec2(1,1,1,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts));
	ReferenceSpecMat = zeros(size(ReferenceSpecMat_Spec,1));

	CircShiftVec = SearchForPeak_LeftPt_Pts-SearchForPeak_Center_Pts :1: SearchForPeak_RightPt_Pts-SearchForPeak_Center_Pts;
	for i = 1:abs(SearchForPeak_RightPt_Pts-SearchForPeak_LeftPt_Pts+1)
		ReferenceSpecMat(i,:) = circshift(ReferenceSpecMat_Spec,CircShiftVec(i));
	end
	
end


%% 3. Calculate & Apply ShiftMap


OutArray = fftshift(fft(OutArray,[],ApplyAlongDim),ApplyAlongDim);

for x = 1:size(OutArray,1)
	for y = 1:size(OutArray,2)
		for z = 1:size(OutArray,3)
			
			if(mask(x,y,z) == 0 || (RefSpecWasProvided && x==RefVox(1) && y == RefVox(2) && z == RefVox(3)))			% Dont process if mask=0 or reference voxel
				continue
			end
			
			
			if(~OnlyApplyShiftMap_flag)
				% Calculate ShiftMap
				DotProd = abs(ReferenceSpecMat) * abs(squeeze(SearchArray(x,y,z,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts)));
                MaxDotProdMatch = find(DotProd == max(DotProd)); MaxDotProdMatch = MaxDotProdMatch(1);
				ShiftMap(x,y,z) = -round( CircShiftVec(MaxDotProdMatch) / ZerofillingFactor); % - because we shifted the reference, but below we want to shift the other spectra
			end
			
			% Apply ShiftMap
			OutArray(x,y,z,:) = circshift(squeeze(OutArray(x,y,z,:)),[ShiftMap(x,y,z) 1]);

			
		end
	end
end

OutArray = ifft(fftshift(OutArray,ApplyAlongDim),[],ApplyAlongDim);








%% 5. Postparations







