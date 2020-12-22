function [MRStruct, AdditionalOut] = op_CorrSpectralB0(MRStruct,B0,Mask,Settings)
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
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadMRStructSets)
%
% Input: 
% -         file                    ...     Path of MRS(I) file.
% -         DesiredSize             ...     If you want to perform zerofilling or use only a part of the kspace, set
%                                           DesiredSize to your wanted [kx,ky,kz], e.g. [64 64 1].
% -         ReadMRStructSets          ...     
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

if(isempty(Mask))
    clear Mask
end
if(~exist('Settings','var'))
    Settings = [];
end
if(~isfield(Settings,'Debug_flag'))
    Settings.Debug_flag = false;
end
if(~isfield(Settings,'Overdiscrete_flag'))
    Settings.Overdiscrete_flag = false;
end
if(Settings.Overdiscrete_flag && ~isfield(Settings,'OverdiscreteSize'))
    Settings.OverdiscreteSize = size_MultiDims(B0.B0Map,1:3);
end
if(~Settings.Overdiscrete_flag || ~isfield(Settings,'Downsample_flag'))
    Settings.Downsample_flag = false;
end
if(~isfield(Settings,'ReverseB0_flag'))
    Settings.ReverseB0_flag = false;
end
if(~isfield(Settings,'RoundB0ToIntVecPts'))
    Settings.RoundB0ToIntVecPts = false;
end
if(~isfield(MRStruct.Par,'DataSize'))
    MRStruct.Par.DataSize = size(MRStruct.Data);
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end
if(~isfield(MRStruct.RecoPar,'DataSize'))
    MRStruct.RecoPar.DataSize = size(MRStruct.Data);
end


%% Calculate B0-Map 
% From Residual water and lipids in case B0-map is not provided

if(~exist('B0','var') || isempty(B0))
    TestIn.csi = MRStruct.Data;
    TestIn.mask = MRStruct.Mask;
    Sett = struct('PeakSearchPPM',[4.7],'PolyfitRegion',[4.4 5.0],'PeakSearchRangePPM',0.5); 
    Sett.LarmorFreq = MRStruct.RecoPar.LarmorFreq; Sett.Dwelltime = MRStruct.RecoPar.Dwelltimes(1); Sett.vecsize = MRStruct.RecoPar.vecSize; 
    TestIn.csi = MRStruct.Data; TestIn.mask = MRStruct.Mask;
    [TestOut,ShiftMap] = FrequencyAlignment(TestIn,Sett,4,2); 
    MRStruct.Data = TestOut;
end


%% Reverse B0

if(Settings.ReverseB0_flag)
    B0.B0Map = -1*B0.B0Map;
end


%% DEBUG: Plot Spectra Before B0Corr

if(Settings.Debug_flag)
    bla = MRStruct.Data(:,:,1,:,1);
    bla = bla .* myrepmat(MRStruct.Mask,size(bla));
    bla = abs(fftshift(fft(bla,[],4),4));
    bla = permute(bla,[4 1 2 3]);
    chemy = compute_chemshift_vector_1_2(MRStruct.RecoPar.LarmorFreq,MRStruct.RecoPar.Dwelltimes(1)/10^9,MRStruct.RecoPar.vecSize);
    figure; plot(chemy,bla(:,:))
    figure; plot(chemy,sum(bla(:,:),2))
    clear bla chemy
end

%% Get to Same Resolution for B0-map and MRStruct

if(Settings.Overdiscrete_flag)
    MRStruct.Data = ZerofillOrCutkSpace(MRStruct.Data,[Settings.OverdiscreteSize MRStruct.RecoPar.DataSize(4:end)],1);
    if(isfield(MRStruct,'NoiseData') && numel(MRStruct.NoiseData) > 1)
        MRStruct.NoiseData = ZerofillOrCutkSpace(MRStruct.NoiseData,[Settings.OverdiscreteSize MRStruct.RecoPar.DataSize(4:end)],1);       
    end
end

if(ndims(B0.B0Map) == 3)
    CurB0 = imresize3(B0.B0Map,size_MultiDims(MRStruct.Data,1:3));
else   
    CurB0 = imresize(B0.B0Map,size_MultiDims(MRStruct.Data,1:2));
end
if(exist('Mask','var'))
    if(ndims(Mask) == 3)
        Mask = imresize3(double(Mask),size(CurB0),'nearest');
    else
        Mask = imresize(double(Mask),size(CurB0),'nearest');        
    end
    CurB0 = CurB0 .* Mask;
end

% CAREFUL: FROM HERE ON MRStruct.RecoPar.DataSize MIGHT HAVE NOT THE REAL SIZE OF THE MRStruct.Data !!!


%% Apply B0Mat to D1Cart
% Round to shift only integer number of points
HzPerPt = 10^9/MRStruct.RecoPar.Dwelltimes(1) / MRStruct.RecoPar.vecSize;
if(Settings.RoundB0ToIntVecPts)
    CurB0 = round(CurB0/HzPerPt)*HzPerPt;
end

% Remove NaN's
CurB0(isnan(CurB0)) = 0;

AdditionalOut.B0Map = CurB0;
time   = (0:MRStruct.RecoPar.DataSize(4)-1)*MRStruct.RecoPar.Dwelltimes(1)/10^9;
AdditionalOut.TimeVec = time;
B0CorrMat_Spec = exp(myrepmat(2*pi*1i*CurB0,size(MRStruct.Data)) .* myrepmat(time,size(MRStruct.Data)));

MRStruct.Data = MRStruct.Data .* B0CorrMat_Spec;
if(isfield(MRStruct,'NoiseData') && numel(MRStruct.NoiseData) > 1)
    MRStruct.NoiseData = MRStruct.NoiseData .* B0CorrMat_Spec;    
end

%% Save B0CorrMat_Spec

if(nargout > 1)
    AdditionalOut.B0CorrMat_Spec = B0CorrMat_Spec;
end


%% Downsample Data
if(Settings.Overdiscrete_flag)
    if(Settings.Downsample_flag)
        MRStruct.Data = ZerofillOrCutkSpace(MRStruct.Data,MRStruct.RecoPar.DataSize,1);
        if(isfield(MRStruct,'NoiseData') && numel(MRStruct.NoiseData) > 1)
            MRStruct.NoiseData = ZerofillOrCutkSpace(MRStruct.NoiseData,MRStruct.RecoPar.DataSize,1);       
        end
    else
        MRStruct.RecoPar.DataSize(1:3) = size_MultiDims(B0.B0Map,1:3);
    end
end

%% DEBUG: Plot Spectra After B0Corr

if(Settings.Debug_flag)
    plot_SpecOfAllVoxels(MRStruct);
end


%% Postparations

MRStruct = supp_UpdateRecoSteps(MRStruct,Settings);

