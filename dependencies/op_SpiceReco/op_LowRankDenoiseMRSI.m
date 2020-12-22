function [InData, AdditionalOut] = op_LowRankDenoiseMRSI(InData,ExternalNoise,Settings)
%
% read_csi_dat Read in raw data from Siemens
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadInDataSets          ...     
%
% Output:
% -         kSpace                      ...     Output data in k-space. In case of SVS this is zero. size: channel x ROW x COL x SLC x Samples x Averages
% -         Info                        ...     
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: This function uses FFTs to get from k- to image-space. This is mathematically wrong, but Siemens seems to do the same when
% creating DICOMS. The only difference is that the images are flipped left/right and up/down.

% This function expects the input to be of form
% [nCha, nAngInt 


%% 0. Preparations


if(isempty(ExternalNoise))
    clear ExternalNoise
end
if(~isfield(Settings,'Debug_flag'))
    Settings.Debug_flag = false;
end
if(~isfield(Settings,'MaskShape'))
   Settings.MaskShape = 3; 
end
if(~isfield(InData,'RecoPar'))
    InData.RecoPar = InData.Par;
end

%%


% Estimate Noise
if(~exist('ExternalNoise','var'))
    NoiseVec = permute(InData.Data,[5 1 2 3 6 4]);
else
    NoiseVec = permute(ExternalNoise,[5 1 2 3 6 4]);
end
NoiseVec = GatherNoiseFromCSI(NoiseVec,Settings.MaskShape,200,10);      % Gather noise from end of FID


% Calculate Std
NoiseVecStd = complex(std(real(NoiseVec),[],2), std(imag(NoiseVec),[],2));

% Create noise with same size as InData.Data and std as defined by NoiseVecStd
Noise = complex(  myrepmat(real(NoiseVecStd),size(InData.Data)).*randn(size(InData.Data)),...
myrepmat(imag(NoiseVecStd),size(InData.Data)).*randn(size(InData.Data))  );


%% 


% Rearrange Data to be  [ S(r1,t1)  ... S(r1,tN) ]
%                       [ S(r2,t1)  ... S(r2,tN) ]
%                                   ...
%                       [ S(rM,t1)  ... S(rM,tN) ]
SVs = reshape(InData.Data.* InData.BrainMask,[prod(InData.RecoPar.DataSize(1:3)) InData.RecoPar.DataSize(4) prod(InData.RecoPar.DataSize(5:end))]); 

Min = min(size(SVs,1),size(SVs,2));
AdditionalOut.U = zeros([size(SVs,1) Min InData.RecoPar.DataSize(5:end)]);
AdditionalOut.Sigma = zeros([Min Min InData.RecoPar.DataSize(5:end)]);
AdditionalOut.V = zeros([size(SVs,2) Min InData.RecoPar.DataSize(5:end)]);

for AdditionalDimInd = 1:size(SVs,3)
    
    % Do an SVD
    [AdditionalOut.U(:,:,AdditionalDimInd), AdditionalOut.Sigma(:,:,AdditionalDimInd), AdditionalOut.V(:,:,AdditionalDimInd)] ...
    = svd(SVs(:,:,AdditionalDimInd),'econ');
    

    % Do SVD of noise
    CurNoise = Noise(:,:,:,:,AdditionalDimInd);
    CurNoise = CurNoise .* myrepmat(InData.BrainMask,size(CurNoise));
    NoiseSVs = reshape( CurNoise,[prod(InData.RecoPar.DataSize(1:3)) InData.RecoPar.DataSize(4)]);
    NoiseSVs = svd(NoiseSVs,'econ');


    % Plot Singular Values
    if(Settings.Debug_flag)
%         figure; imagesc( ( std(real(InData.Data(:,:,:,end-300:end)),0,4) + std(imag(InData.Data(:,:,:,end-300:end)),0,4) ) / 2)
        figure; plot(diag(AdditionalOut.Sigma(:,:,AdditionalDimInd)),'b'); hold on; plot(NoiseSVs,'r'); hold off;
        title('Singular Values B0Corr'); legend('SVs(B0Corr)','SVs(Noise)')
    end


    RankL_Suggested = find(diag(AdditionalOut.Sigma(:,:,AdditionalDimInd)) < NoiseSVs(1),1);
    RankL_Suggested(isempty(RankL_Suggested)) = Inf %#ok
    RankL = RankL_Suggested;
    if(isempty(RankL) || RankL > Settings.MaxRankL)
        fprintf('\nRestricting rank from %d to %d',RankL, Settings.MaxRankL);
        RankL = Settings.MaxRankL;
    end
    if(isfield(Settings,'MinRankL') && RankL <= Settings.MinRankL)
        RankL = Settings.MinRankL;
    end
    AdditionalOut.RankL(AdditionalDimInd) = RankL;
    
    
    % Resynthesize InData
    InData.Data(:,:,:,:,AdditionalDimInd) = reshape(AdditionalOut.U(:,1:RankL,AdditionalDimInd) * ...
    AdditionalOut.Sigma(1:RankL,1:RankL,AdditionalDimInd) * AdditionalOut.V(:,1:RankL,AdditionalDimInd)',InData.RecoPar.DataSize(1:4));
    

end



