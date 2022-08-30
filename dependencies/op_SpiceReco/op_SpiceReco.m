function [SpiceOut, AdditionalOut] = op_SpiceReco(D1,D2,B0,ModelFunction,SpiceOperators,Settings)
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

if(~isfield(Settings,'UseExternalPhi'))
    Settings.UseExternalPhi.Flag = false;
end
if(~isfield(Settings,'FreqAlign2_D1'))
    Settings.FreqAlign2_D1.Flag = true;
    Settings.FreqAlign2_D1 = struct('ApplyAlongDim',4,'ZerofillingFactor',2,'PeakSearchPPM',2.94,'PolyfitRegion',[2.3 3.4],'PeakSearchRangePPM',0.2);
end
if(~isfield(Settings,'ExpFilter_D1'))
    Settings.ExpFilter_D1.Flag = false;
end
if(~isfield(Settings,'NuisRem_D1'))
    Settings.NuisRem_D1{1}.Flag = false;
end
if(Settings.UseExternalPhi.Flag)
    Settings.FreqAlign2_D1.Flag = false;
    Settings.ExpFilter_D1.Flag = false;
    Settings.CorrSpectralB0_D1.Flag = false;
    for NuisLoop = 1:numel(Settings.NuisRem_D1); Settings.NuisRem_D1{NuisLoop}.Flag = false; end    
end

%% Remove NuisanceSignal from D2

% The functions used for nsRmSparse_bstrchanged have same name, but are different. This is a temporary hack
rmpath('./dependencies/op_NuisanceRemoval_HSVD')
addpath('./dependencies/op_NuisanceRemoval_HSVD_Sparse')



if(Settings.NuisRem_D2.NuisRem_flag)

    params.tD1 = (0:D1.RecoPar.Dwelltimes(1):D1.RecoPar.Dwelltimes(1)*(D1.RecoPar.vecSize-1))/10^9;
    
    % Works only for single temporal interleaf for now!!!
    params.tD2 = zeros(1,D2.Par.DataSize(1)*D2.Par.DataSize(end-1));
    for CurVecPt = 1:D2.Par.DataSize(end-1)
        params.tD2(  (CurVecPt-1)*D2.Par.DataSize(1)+1  : CurVecPt*D2.Par.DataSize(1)  ) = (CurVecPt-1)*D2.Par.TrajTotPts + (0:(D2.Par.DataSize(1)-1));
    end
    params.tD2 = params.tD2*(D2.Par.ADC_dt/10^9);
%     params.tD2 = (0:D2.Par.Dwelltimes(1):D2.Par.Dwelltimes(1)*(D2.Par.vecSize-1))/10^9;
    params.ND1 = D1.RecoPar.DataSize([1 2 4]);
    params.ND2 = SpiceOperators.OutDataSize(1:2);
    params.ND2In = SpiceOperators.InDataSize(1:2);
    params.ImMask = B0.Mask;
    params.waterMask = B0.Mask;    
    params.B0Map = -B0.B0Map;
    params.NosD1 = 2;
    params.dt = D1.RecoPar.Dwelltimes(1)/10^9;
    
    params.lambdaDfWat = 0; params.lambdaDfLip = 0; params.lambdaDfMet = 0;

    
%     params.csiInfo.Ny = params.ND1(2); params.csiInfo.Nx = params.ND1(1); params.csiInfo.Nt = params.ND1(3);
%     params.csiInfo.dt = D1.RecoPar.Dwelltimes(1)/10^9; 
    params.csiInfo.NtEcho = 0;
    
    params = nsRmSparseInitParams(params);
    
    params.xtD1W = squeeze(D1.Data);
    params.xtD1WrL = params.xtD1W;
    params.xtD1WrLr = params.xtD1W;
    params.waterMaskD2 = double(params.waterMaskD2);
    params.waterMaskD2(params.waterMaskD2 == 0) = 0.01;
    params.lipMaskD2(params.lipMaskD2 == 0) = 0.01;
    
    
    Factor = D1.Par.LarmorFreq * 1e-6;
    ConjSign = ...
    (-1)^(double(D1.RecoSteps.Step6_op_IterReconstructNonCartMRData.ConjInkSpace_flag) + double(D1.RecoSteps.Step6_op_IterReconstructNonCartMRData.ConjIniSpace_flag));
    params.optWat.fSel = [80 -80];
    params.optWat.fSel2 = ConjSign*(Settings.NuisRem_D1{1}.WaterPPMs([1 2 4]) - 4.65)*Factor;
    params.optWat.maxT2 = [10 1000000];
    params.optLip.fSel2 = ConjSign*(Settings.NuisRem_D1{1}.LipidPPMs - 4.65)*Factor;
    params.optLip.maxT2 = Settings.NuisRem_D1{1}.LipidT2s;
    params.optMeta.fSel2 = ConjSign*(Settings.NuisRem_D1{1}.MetaboPPMs - 4.65)*Factor;
%     params.optMeta.maxT2 = 25;

    Sizze = SpiceOperators.InDataSize;
    params.sampleMask = reshape(eye(Sizze(1)),[Sizze(1) 1 1 1 Sizze(1)]);
    params.sampleMask = repmat(params.sampleMask,[1 Sizze(2:4) 1 Sizze(5)]);
    params.sampleMask = reshape(params.sampleMask,[Sizze(1:4) Sizze(1)*Sizze(5)]);
%     params.sampleIndex = find(params.sampleMask);
    params.sampleIndex = []; 
    for vecPt = 1:Sizze(5)
        for AIPts = 1:Sizze(2)
            TrajPt = 1:Sizze(1);
            params.sampleIndex = cat(2,params.sampleIndex,((Sizze(1)*(vecPt-1)+TrajPt)-1)*Sizze(1)*Sizze(2)+TrajPt + (AIPts-1)*Sizze(1));
        end
    end

    params.sampleIndexShift = params.sampleIndex;
    
    params.Operators = SpiceOperators;
    params.Operators.SensMap = 1;
    params.Operators.InDataSize(6) = 1;
    
    params.verbose = 0;
    
    params.RW = 15; 
    params.RL = 0;
    params.RM = 3;
    
    params.MaxIterCG = 10;
    
    
    
    
    for CurCha = 1:size(D2.Data,6)
        fprintf('\n\nDecontaminating Channel %d\n',CurCha)
        d = (D2.Data(:,:,:,:,:,CurCha));
        [ dWL, D2.Data(:,:,:,:,:,CurCha), output, params ] = nsRmSparse_bstrchanged( d, params );
%         D2.Data(:,:,:,:,:,CurCha) = D2.Data(:,:,:,:,:,1);
    end
end


rmpath('./dependencies/op_NuisanceRemoval_HSVD_Sparse')
addpath('./dependencies/op_NuisanceRemoval_HSVD')



%% Correct Spectral B0 for D1

if(Settings.CorrSpectralB0_D1.Flag)
    D1 = op_CorrSpectralB0(D1,B0,D1.Mask,Settings.CorrSpectralB0_D1);
    % Alternative: Use FrequencyAlignment on the residual water. Interesting observation: If I nicely frequency-align the spectra according to the residual water,
    % my metabolites don't really get properly frequency aligned. I guess this has to do with the water suppression, which can make the water-maximum appear at 
    % slightly different places... So I don't recomend to do this alternative. Instead, first apply the normal B0-map, then remove lipids and water, and then 
    % align the Cr-peaks (NAA or Cho might also work).
%     TestIn.csi = D1.Data;
%     TestIn.mask = D1.Mask;
%     Sett = struct('PeakSearchPPM',[4.7],'PolyfitRegion',[4.4 5.0],'PeakSearchRangePPM',0.5); Sett.LarmorFreq = D1.RecoPar.LarmorFreq; Sett.Dwelltime = D1.RecoPar.Dwelltimes(1); Sett.vecsize = D1.RecoPar.vecSize; TestIn.csi = D1.Data; TestIn.mask = D1.Mask;
%     [TestOut,ShiftMap] = FrequencyAlignment(TestIn,Sett,4,2); 
%     D1.Data = TestOut;
    
end


%% DEBUG: Check Alignment

% plot_SpecOfAllVoxels(D1)



%% Remove NuisanceSignal of D1

mask = D1.BrainMask(:,:,:,1);
% mask = true(size(mask2));

for NuisLoop = 1:numel(Settings.NuisRem_D1)
    if(Settings.NuisRem_D1{1}.Flag)
        D1 = op_NuisanceRemoval_HSVD(D1,Settings.NuisRem_D1{NuisLoop},mask,0);
    end
end
AdditionalOut.D1 = D1;

clear mask2 mask;


%% DEBUG: Check NuisRem

% D1Temp = D1;
% D1Temp.Data = D1Temp.Data .* D1.BrainMask;
% plot_SpecOfAllVoxels(D1Temp); clear D1Temp



%% Frequency Align D1-Data

if(Settings.FreqAlign2_D1.Flag)
    [D1,ShiftMap] = op_FrequencyAlign2(D1,Settings.FreqAlign2_D1);
end



%% Exponentially Filter D1-Data

if(Settings.ExpFilter_D1.Flag)
    D1 = op_ExponentialFilter(D1,Settings.ExpFilter_D1);
end


%% Extract Phi (Spectral Basis) from D1

if(Settings.UseExternalPhi.Flag)
    USV.RankL = size(D1.Data,1);
    USV.RankL = min(USV.RankL,Settings.Denoise_D1.MaxRankL);
    USV.V = D1.Data(1:USV.RankL,:,:)';
    USV.U = ones(USV.RankL);        % Should not matter, I hope
    USV.Sigma = ones(USV.RankL);    % Should not matter, I hope
else
    [bla0,USV] = op_LowRankDenoiseMRSI(D1,[],Settings.Denoise_D1);
end

%% Interpolate Phi to FID-raster of Target Result
% if(D1.RecoPar.Dwelltimes(1) ~= (D2.Par.Dwelltimes(1)/D2.Par.TimeUndersamplFactor) || D1.RecoPar.DataSize(4) ~= D2.Par.DataSize(5))

if(D1.RecoPar.Dwelltimes(1) ~= D2.Par.Dwelltimes(1) || D1.RecoPar.DataSize(4) ~= D2.Par.DataSize(5))
    tD1 = (0:D1.RecoPar.vecSize-1)*D1.RecoPar.Dwelltimes(1)/10^9;
    tD2 = (0:D2.Par.vecSize-1)*D2.Par.Dwelltimes(1)/10^9/D2.Par.TimeUndersamplFactor;                                                            % No longer need TimeUndersamplFactor, bc
%     tD2 = (0:D2.Par.vecSize*D2.Par.TimeUndersamplFactor-1)*D2.Par.Dwelltimes(1)/10^9/D2.Par.TimeUndersamplFactor;  % we no longer remove data, but zero them 
    for CurCha = 1:size(USV.V,3)
        Ranky = USV.RankL(CurCha);
        CurU = USV.U(:,1:Ranky,CurCha);
        CurS = USV.Sigma(1:Ranky,1:Ranky,CurCha);
        CurV = USV.V(:,1:Ranky,CurCha)';
        AdditionalOut.Phi{CurCha} = InterpFID(CurV,CurS,CurU,tD1,tD2,false);  
    end
    clear CurU CurS CurV Ranky CurCha;
end


%% Create Phi if Not Existing
% (Happens when no interpolation is done in previous step)

if(~isfield(AdditionalOut,'Phi'))
    for CurCha = 1:size(USV.V,3)
        AdditionalOut.Phi{CurCha} = USV.V(:,1:USV.RankL(CurCha),CurCha)';
    end
end
AdditionalOut.RankL = USV.RankL;
clear USV;


%% Define SamplingOperator for D2

if(~isfield(SpiceOperators,'SamplingOperator'))
    SpiceOperators.SamplingOperator = zeros([D2.Par.DataSize(1:4) SpSpice.D1.Par.DataSize(5)]);
    SpiceOperators.SamplingOperator(:,:,:,:,1:D2.Par.TimeUndersamplFactor:end) = 1;

    % Randomly undersample in kx and ky
    % test = zeros([80 80]);
    % test(1:40,1:80) = 1;
    % test_access = randperm(80*80);
    % test2 = repmat(reshape(test(test_access),[80 80]),[1 1 1 SpSpice.D1.Par.DataSize(end)]);
    % D2.SamplingOperator(logical(test2)) = 1;
end




%%

m1 = 62; m2 = 62;
% delta = 1/600;
eta = 0.05;     

% Im_HighRes = imresize(B0.Mag,[m1 m2]);
% Im_HighRes( isnan( Im_HighRes(:) ) ) = 0;
% [Gx, Gy]    = gradient(abs(Im_HighRes));
% Im_HighRes_Grd  = sqrt( Gx.^2 + Gy.^2 );
% lambdaTV = 4.7E-3; % 10
% MaxwTV = 0.5;
% maskk = MaskShrinkOrGrow(reshape(SpiceOperators.Mask,[m1 m2]),0,1,1);
% maxx = Im_HighRes_Grd .* maskk; maxx(maxx == 0) = NaN; maxx = nanmax(maxx(:));
% Im_HighRes_Grd = Im_HighRes_Grd/maxx * 6*lambdaTV;
% % Formula 1 for calculating wTV
% wTV = (lambdaTV.^2+MaxwTV*lambdaTV*Im_HighRes_Grd)./sqrt(Im_HighRes_Grd.^2+lambdaTV^2); 
% % Formula 2 for calculating wTV
% wTV = (lambdaTV.^2)./sqrt(Im_HighRes_Grd.^2+lambdaTV^2); 
% % Formula 3 for calculating wTV (Xianqis original way)
% wTV = lambdaTV./sqrt(Im_HighRes_Grd.^2+lambdaTV^2); 
% % Additional thresholding, works well together with Formula 3
% wTV(wTV > 2*lambdaTV) = 2*lambdaTV;
% wTV(wTV < 0.3*lambdaTV) = 0.99*lambdaTV;
% % Outside of mask: Set to lambdaTV, to avoid massive noise 
% wTV(~maskk) = lambdaTV;

% % % Xianqis option to normalize the image
% % Im_HighRes_hdr = load_nifti( Im_HighRes_filename );
% % Im_HighRes     = Im_HighRes_hdr.vol;
% % Im_HighRes( isnan( Im_HighRes(:) ) ) = 0;
% % norm_Im_HighRes = sqrt(sum(Im_HighRes .^ 2, 3));   %// Calculate Euclidean length
% % norm_Im_HighRes(norm_Im_HighRes < eps) = 1;       %// Avoid division by zero
% % Im_HighRes = bsxfun(@rdivide, Im_HighRes, norm_Im_HighRes); %// Normalize
% % [Gx, Gy, Gz]    = gradient(Im_HighRes);
% % Im_HighRes_Grd  = sqrt( Gx.^2 + Gy.^2 + Gz.^2 );
% % eta = 1e-2;
% % wTV = eta./sqrt(Im_HighRes_Grd.^2 + eta^2);



% % TV without weighting
% wTV = 7E-3;
wTV = 1E-3;

% wTV_0 = 1.0e-1; %1e-2;
% wP = 10*wTV_0;


%% Reconstruct Spatial Map UTS (= U Times Sigma)


tic_IterReco = tic;
fprintf('\n\nStart Iterative Reco . . .')




%% Reconstruct Spatial Map UTS (= U Times Sigma)


% % Xianqis Methods
% 
% 
% % Define Input Data
% % The ratio of the FoV-overgrids is necessary, because the data is scaled by the matrix size when going Cartesian --> Spiral, but
% % the Cartesian data is zerofilled in image domain, and those zeros dont contribute any signal...
% % CartSize = SpSpice.GroundTruth_ForSpSpice.Par.DataSize(1);
% % sft2_Operator_ForSpSpice = sft2_Operator(transpose(squeeze(SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM)*CartSize),transpose(nsamp),1);
% % 
% % if(Ctrl.RecoPar.Correct4B0_flag)
% %     sft2_Operator_ForSpSpice = sft2_Operator_ForSpSpice .* SpSpice.Reco.B0CorrMat_Spatial;
% % end
% % 
% % % B0CorrMat = ones([prod(SpSpice.Reco.Par.DataSize(1:end-1)) SpSpice.Reco.Par.DataSize(end)]);
% % B0CorrMat = conj(reshape(SpSpice.Reco.B0CorrMat_Spec,[prod(SpSpice.Reco.Par.DataSize(1:end-1)) SpSpice.Reco.Par.DataSize(end)]));
% 
% SpiceOperators.SamplingIndex = find(repmat(SpiceOperators.SamplingOperator,[1 1 1 1 1 SpiceOperators.InDataSize(end)]));
% 
% for CurAdditInd = 1:numel(AdditionalOut.Phi)
%     
% %     CurD2 = D2.Data(:,:,:,:,:,CurAdditInd);
%     CurD2 = D2.Data;
%     CurD2 = reshape(CurD2, [numel(CurD2) 1]);
% 
% %     % For pcg etc.
% %     CurD2 = ModelFunction('Transj',CurD2,AdditionalOut.Phi{CurAdditInd},SpiceOperators);
% %     funh = @(x) ModelFunction('Transj', ModelFunction('NoTransj',x,...
% %     AdditionalOut.Phi{CurAdditInd},SpiceOperators),      ...
% %     AdditionalOut.Phi{CurAdditInd},SpiceOperators);
% %     InitGuess = [];
% 
% 
% % %     lsqr
% %     funh = @(x,transp_flag) ModelFunction(transp_flag, x, AdditionalOut.Phi{CurAdditInd},SpiceOperators);
% %     [AdditionalOut.UTS{CurAdditInd}] = lsqr(funh,CurD2,Settings.IterReco.Tolerance,Settings.IterReco.Iterations,[],[],InitGuess);
% 
% %     Xianqis Least Square
%     funh_T = @(x) ModelFunction('transp', x, AdditionalOut.Phi{CurAdditInd},SpiceOperators);
%     funh = @(x) ModelFunction('notransp', x, AdditionalOut.Phi{CurAdditInd},SpiceOperators);
%     eta = 1.0E0;
% %     eta = 1E-2; % w/o DCF
%     AdditionalOut.UTS{CurAdditInd} = LS_method(CurD2,funh,funh_T,Settings.IterReco.Iterations,eta);
%     
% 
% %     % Xianqis weighted TV method 
% %     funh_T = @(x) ModelFunction('transp', x, AdditionalOut.Phi{CurAdditInd},SpiceOperators);
% %     funh = @(x) ModelFunction('notransp', x, AdditionalOut.Phi{CurAdditInd},SpiceOperators);
% %     AdditionalOut.UTS{CurAdditInd} = wTV_method(CurD2,funh,funh_T,Settings.IterReco.Iterations,wTV,eta,m1,m2,AdditionalOut.RankL);
% 
%     
% 
%     AdditionalOut.UTS{CurAdditInd} = reshape(AdditionalOut.UTS{CurAdditInd},[prod(SpiceOperators.OutDataSize(1:2)) AdditionalOut.RankL(CurAdditInd)]);
% 
% %     D2.SamplingOperator = reshape(D2.SamplingOperator,D2.Par.DataSize);
% 
% end

% Antoines Reco

% Define Parameters
mrsiReconParams.modelOrder = AdditionalOut.RankL;
mrsiReconParams.maxit = Settings.IterReco.Iterations;
mrsiReconParams.minit = 15; mrsiReconParams.minit(mrsiReconParams.minit > mrsiReconParams.maxit) = ceil(mrsiReconParams.maxit/2);
mrsiReconParams.Threshold=1E-8;

mrsiReconParams.LRTGVModelParams.check_it = 1;
mrsiReconParams.LRTGVModelParams.Plot_it = 5;
mrsiReconParams.LRTGVModelParams.Orthogonalize_it = 9999;%=20 %Orthogonalize Sepctral component  every *this* step
mrsiReconParams.LRTGVModelParams.SpecItFact = 0;%2  % Will run *this* more spectral gradient descent than spatial convergence steps
mrsiReconParams.LRTGVModelParams.reduction = 100;% 100 or 1000 if diverge. Starting with  alpha0/reduction, alpha1/reduction at the begining then going back to alpha0, alpha1 original values
mrsiReconParams.LRTGVModelParams.min_SpectStep = 1/128;%=1/128 Minimum Spectral step size
mrsiReconParams.LRTGVModelParams.max_SpectStep = 1/4;%=1/4 % decrease if diverges
mrsiReconParams.LRTGVModelParams.min_taup = 1/32;%=1/128 Minimum Spatial step size   (Used: 1/32)
mrsiReconParams.LRTGVModelParams.max_taup = 8/16;%=1/16 % decrease if diverges       (Used: 8/16)

mrsiReconParams.mrProt.samplerate = 10^9/D2.Par.Dwelltimes(1);
mrsiReconParams.mrProt.VSize = D2.Par.vecSize;
mrsiReconParams.mu_tv = 2E-4;

mrsiReconParams.CSOperators = SpiceOperators;


RecoOper = @(X) ModelFunction('NoTransj', X, mrsiReconParams.CSOperators);
RecoOper_Adj = @(X) ModelFunction('Transj', X, mrsiReconParams.CSOperators);

for CurAdditInd = 1:numel(AdditionalOut.Phi)
    
    
    LowRankRecoSett = struct('WriteFiles_flag',false,'WriteFiles_Path','./DebugOutput/LowRankTGVRecon','FixedV_flag',true,'V',AdditionalOut.Phi{CurAdditInd}');
    if(iscell(D2.Data))
        CurD2 = cat(1,D2.Data{:});
    else
        CurD2 = D2.Data; 
    end

    [ResultData_trr,U_rrc,V_tc,S] = LowRankTGV_RecurModel_bstr(CurD2,mrsiReconParams,RecoOper,RecoOper_Adj,LowRankRecoSett);
%     [ResultData_trr,U_rrc,V_tc,S] = LowRankTGV_RecurModel_bstr(CurD2,mrsiReconParams,RecoOper,RecoOper_Adj,struct('FixedV_flag',false));
    AdditionalOut.UTS{CurAdditInd} = reshape(U_rrc,[numel(U_rrc)/size(S,1) size(S,1)]) * S;

end

%% Resynthesize Data Using Estimated UTS

SpiceOut.Data = zeros(Settings.ResynthesizeSVDData.OutSize);
for CurAdditInd = 1:numel(AdditionalOut.Phi)

    Temp = AdditionalOut.UTS{CurAdditInd} * AdditionalOut.Phi{CurAdditInd};
    SpiceOut.Data(:,:,:,:,CurAdditInd) = reshape(Temp,Settings.ResynthesizeSVDData.OutSize(1:4));

end

fprintf('\ntook %f s',toc(tic_IterReco))





