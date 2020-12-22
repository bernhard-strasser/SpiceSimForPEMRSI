function [MRStruct,AdditionalOut] = op_FrequencyAlign2(MRStruct,Settings,AdditionalIn)
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
% [MRStruct,mask] = FrequencyAlignment(MRStruct,ApplyAlongDims,EllipsoidCoefficients,PerformFFT_flag)
%
% Input: 
% -     MRStruct                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims               ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients        ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     PerformFFT_flag              ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
%
% Output:
% -     MRStruct                     ...     The filtered/masked output array
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
	MRStruct = 0;
	return
end
% if(nargin < 2)
%     return
% end
if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
    if(isfield(MRStruct,'Data_file'))
        MRStruct.Par = read_ascconv(MRStruct.Data_file); 
    end
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end

size_InArray = size(MRStruct.Data);

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'ApplyAlongDim'))
    Settings.ApplyAlongDim = find(size_InArray == max(size_InArray));
end
if(~isfield(Settings,'ZerofillingFactor'))
    Settings.ZerofillingFactor = 2;
end
if(~isfield(Settings,'PeakSearchPPM'))
    Settings.PeakSearchPPM = 3.02;  % Default: Cr
end
if(~isfield(Settings,'PeakSearchRangePPM'))
    Settings.PeakSearchRangePPM = 0.2;
end
if(~isfield(Settings,'UseSVDForRef_flag'))
    Settings.UseSVDForRef_flag = true;
end
if(~isfield(Settings,'UseSVDForRef_NoSingVals'))
    Settings.UseSVDForRef_NoSingVals = 1;
end
if(~isfield(Settings,'MaxShifts_ppm'))
    Settings.MaxShifts_ppm = [-0.1 0.1];     % to left and right direction
end
if(~isfield(Settings,'AlignRefPeak_flag'))
    Settings.AlignRefPeak_flag = 1;
end
if(~isfield(Settings,'PeakDetection'))
    Settings.PeakDetection = struct;
end
if(~isfield(Settings.PeakDetection,'MinSNR'))
    Settings.PeakDetection.MinSNR = 2;
end
if(~isfield(Settings.PeakDetection,'MinDotProd'))
    Settings.PeakDetection.MinDotProd = 0.55;
end
if(~isfield(Settings,'FindShiftAlgorithm'))
    Settings.FindShiftAlgorithm = 'DotProduct';
end
if(~isfield(Settings,'UseThisInStructMask'))
    Settings.UseThisInStructMask = 'BrainMask';
end
if(~exist('AdditionalIn','var'))
    AdditionalIn = struct();
end

if(~isfield(AdditionalIn,'RefSpecIn'))
	RefSpec = [];
else
    RefSpec = AdditionalIn.RefSpecIn;
end

% Handle mask
MaskWasProvided = true;
if(isfield(AdditionalIn,'Mask'))
    mask = AdditionalIn.Mask;
elseif(isfield(Settings,'UseThisInStructMask') && isfield(MRStruct,(Settings.UseThisInStructMask)))
    mask = MRStruct.(Settings.UseThisInStructMask);
elseif(isfield(MRStruct,'BrainMask'))
	mask = MRStruct.BrainMask;
elseif(isfield(MRStruct,'Mask'))
	mask = MRStruct.Mask;    
else
	mask = ones(size_InArray(setdiff(1:numel(size_InArray),Settings.ApplyAlongDim)));
	MaskWasProvided = false;	
end

if(exist('AdditionalIn','var') && isfield(AdditionalIn,'ShiftMap'))
    AdditionalOut.ShiftMap = AdditionalIn.ShiftMap;
	OnlyApplyShiftMap_flag = true;
    RefSpecWasProvided = false;
else
	OnlyApplyShiftMap_flag = false;

end


%% Define Reference Spectrum if necessary

if(~OnlyApplyShiftMap_flag)
    RefSpecWasProvided = true;
    if(isempty(RefSpec))
        RefSpecWasProvided = false;
        
        if(Settings.UseSVDForRef_flag)
            CurData = MRStruct.Data .* mask;    % We don't want to take all the lipids etc.
            [~, ~, V] = svd(reshape(CurData,[numel(CurData)/size(CurData,4) size(CurData,4)]));
            V = V';
            RefSpec = reshape(sum(V(1:Settings.UseSVDForRef_NoSingVals,:),1),[ones([1 Settings.ApplyAlongDim-1]) size(V,1)]); RefVox = [0 0 0];
            clear CurData
        else
        
            if(isfield(Settings,'RefVox') && numel(Settings.RefVox) > 1)
                RefVox = Settings.RefVox;
            elseif(MaskWasProvided)
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
            RefSpec = MRStruct.Data(RefVox(1),RefVox(2),RefVox(3),:);
        end
    end
end






% 0.2 Declarations


% 0.3 Definitions
    
size_SearchArray = size_InArray; size_SearchArray(Settings.ApplyAlongDim) = size_SearchArray(Settings.ApplyAlongDim)*Settings.ZerofillingFactor;
SearchArray = MRStruct.Data;

if(~OnlyApplyShiftMap_flag)
    CS_vec_zf = compute_chemshift_vector(MRStruct.RecoPar.LarmorFreq,MRStruct.RecoPar.Dwelltimes(1)/10^9,MRStruct.RecoPar.vecSize*Settings.ZerofillingFactor);
    SearchForPeak_Center_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM));
  	SearchForPeak_LeftPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM - Settings.PeakSearchRangePPM));
	SearchForPeak_RightPt_Pts = find(min(abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM)) == abs(CS_vec_zf - Settings.PeakSearchPPM + Settings.PeakSearchRangePPM));
end


%% 1. Zerofilling & FFT SearchArray

if(~OnlyApplyShiftMap_flag)
    SearchArray = Zerofilling_Spectral(SearchArray,size_SearchArray,0);
    SearchArray = fftshift(fft(SearchArray,[],Settings.ApplyAlongDim),Settings.ApplyAlongDim);

    RefSpec2 = Zerofilling_Spectral(RefSpec,[ones([1 numel(size_SearchArray)-1]) size_SearchArray(Settings.ApplyAlongDim)],0);
    RefSpec2 = fftshift(fft(RefSpec2,[],Settings.ApplyAlongDim),Settings.ApplyAlongDim);
end


%% Align Ref Spectrum


if(Settings.AlignRefPeak_flag)
    % Find closest local maximum to Settings.PeakSearchPPM
    
    % Find all peaks with prominence higher than 1
    [PeakHght,PeakPos] = findpeaks(abs(squeeze(RefSpec2))/max(abs(RefSpec2(SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts))),'MinPeakProminence',0.3);
%     figure; plot(abs(squeeze(RefSpec2))); hold on, scatter(PeakPos,PeakHght,'r'), hold off        % (Debug)
    
    % Restrict peaks to those inside of the Search-range
    DeletePeaks = PeakPos < SearchForPeak_LeftPt_Pts | PeakPos > SearchForPeak_RightPt_Pts;
    PeakPos(DeletePeaks) = []; PeakHght(DeletePeaks) = []; 
    
    % Use the one that is closest to SearchForPeak_Center_Pts
    tmp = min(abs(PeakPos - SearchForPeak_Center_Pts)) == abs(PeakPos - SearchForPeak_Center_Pts);
    PeakPos = PeakPos(tmp); PeakHght = PeakHght(tmp);
    
    % Circshift RefSpectrum
    RefSpec2 = circshift(RefSpec2,SearchForPeak_Center_Pts-PeakPos,Settings.ApplyAlongDim);
    
    % Debug:
%     figure; plot(abs(squeeze(RefSpec2))); hold on, scatter(PeakPos,PeakHght,'r'), hold off        % (Debug)
    clear tmp PeakPos PeakHght DeletePeaks
end



%% 2. Create Matrix of Differently Shifted RefSpecs
% Instead of shifting the individual spectra of all voxels N times, computing the DotProduct and finding the maximum,
% do instead: Shift the RefSpec once N times, save that, and calculate the DotProducts for all shifts of the RefSpec.
% In the fist case, we would need to do NumbOfVox x N shifts, whereas in the latter we only need N shifts.
if(~OnlyApplyShiftMap_flag)
	
	AdditionalOut.ShiftMap = zeros([size(MRStruct.Data,1) size(MRStruct.Data,2) size(MRStruct.Data,3)]);
    AdditionalOut.DotProdMap = AdditionalOut.ShiftMap;

    MaxShifts_Pts = round(Settings.MaxShifts_ppm * (MRStruct.RecoPar.LarmorFreq/1E6)* MRStruct.RecoPar.vecSize/(1E9/ MRStruct.RecoPar.Dwelltimes(1)));
    
	ReferenceSpecMat_Spec = squeeze(RefSpec2(1,1,1,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts)); 
    ReferenceSpecMat_Spec = ReferenceSpecMat_Spec/norm(ReferenceSpecMat_Spec);
	ReferenceSpecMat = zeros(size(ReferenceSpecMat_Spec,1));

	CircShiftVec = SearchForPeak_LeftPt_Pts-SearchForPeak_Center_Pts :1: SearchForPeak_RightPt_Pts-SearchForPeak_Center_Pts;
	for i = 1:abs(SearchForPeak_RightPt_Pts-SearchForPeak_LeftPt_Pts+1)
		ReferenceSpecMat(i,:) = circshift(ReferenceSpecMat_Spec,CircShiftVec(i));
	end
	
end


%% 3. Calculate & Apply ShiftMap

% MRStruct.Data = fftshift(fft(MRStruct.Data,[],Settings.ApplyAlongDim),Settings.ApplyAlongDim);

for x = 1:size(MRStruct.Data,1)
	for y = 1:size(MRStruct.Data,2)
		for z = 1:size(MRStruct.Data,3)
			
            if(mask(x,y,z) == 0 || (RefSpecWasProvided && x==RefVox(1) && y == RefVox(2) && z == RefVox(3)))			% Dont process if mask=0 or reference voxel
				continue
            end
            
            %if(x==13 && y == 29)
            %    sadf = 1; 
            %end
			
			TmpSpec = squeeze(SearchArray(x,y,z,SearchForPeak_LeftPt_Pts:SearchForPeak_RightPt_Pts));
%             TmpSpec = TmpSpec - fftshift(fft(exp(-transpose(0:numel(TmpSpec)-1)/1) .* ifft(ifftshift(TmpSpec))));
%             TmpSpec = abs(TmpSpec);
            TmpSpec = conj(TmpSpec) / norm(TmpSpec);
            DotProd = abs(ReferenceSpecMat * TmpSpec);
            AdditionalOut.DotProdMap(x,y,z) = max(DotProd);
            if(~OnlyApplyShiftMap_flag && max(abs(TmpSpec)) > Settings.PeakDetection.MinSNR*std(TmpSpec))
				% Calculate ShiftMap
                if(strcmpi(Settings.FindShiftAlgorithm,'DotProduct'))
                    MaxDotProd =  max(DotProd);
                    MaxDotProdMatch = find(DotProd == MaxDotProd); MaxDotProdMatch = MaxDotProdMatch(1);
                    ShiftPoint = -( CircShiftVec(MaxDotProdMatch) / Settings.ZerofillingFactor);
                else
                    [~, MaxPos] = max(abs(abs(TmpSpec))); 
                    ShiftPoint = (ceil(numel(TmpSpec)/2) - MaxPos )/ Settings.ZerofillingFactor; % Max should be always in center of TmpSpec!
                    MaxDotProd = 1E9;    % So that we have no condition on DotProd
                end
                if(MaxDotProd > Settings.PeakDetection.MinDotProd && ShiftPoint > min(MaxShifts_Pts) && ShiftPoint < max(MaxShifts_Pts))
                    AdditionalOut.ShiftMap(x,y,z) = ShiftPoint; % - because we shifted the reference, but below we want to shift the other spectra
                else
                    AdditionalOut.ShiftMap(x,y,z) = NaN; continue;
                end
            end
            
			% Apply ShiftMap
% 			MRStruct.Data(x,y,z,:) = circshift(squeeze(MRStruct.Data(x,y,z,:)),[AdditionalOut.ShiftMap(x,y,z) 1]);
			
		end
	end
end

% MRStruct.Data = ifft(fftshift(MRStruct.Data,Settings.ApplyAlongDim),[],Settings.ApplyAlongDim);
CurMap = AdditionalOut.ShiftMap; CurMap(isnan(CurMap)) = 0;         % Can't multiply with NaN, otherwise data gets NaNed
t = (0:(MRStruct.RecoPar.vecSize-1))/MRStruct.RecoPar.vecSize;
MRStruct.Data = MRStruct.Data .* exp(1i*2*pi*CurMap.*reshape(t,[ones([1 Settings.ApplyAlongDim-1]) numel(t) ones([1 ndims(MRStruct.Data)-Settings.ApplyAlongDim])]));

if(nargout > 1)
    HzPerPt = 1E9/MRStruct.RecoPar.Dwelltimes(1)/MRStruct.RecoPar.vecSize;
    AdditionalOut.B0Map_Hz = AdditionalOut.ShiftMap * HzPerPt;
end





%% 5. Postparations

supp_UpdateRecoSteps(MRStruct,Settings);





