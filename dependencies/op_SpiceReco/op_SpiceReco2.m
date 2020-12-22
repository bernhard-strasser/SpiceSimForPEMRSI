function [SpiceOut, AdditionalOut] = op_SpiceReco2(D1,D2,B0,ModelFunction,OptimizationFunction,SpiceOperators,Settings)
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






    %% Correct Spectral B0 for D1

    if(Settings.CorrSpectralB0_D1.Flag)
        D1 = op_CorrSpectralB0(D1,B0,D1.Mask,Settings.CorrSpectralB0_D1);
    end


    %% DEBUG: Check Alignment
    % bla = D1.Data(:,:,1,:,1);
    % bla = bla .* myrepmat(D1.Mask,size(bla));
    % bla = abs(fftshift(fft(bla,[],4),4));
    % bla = permute(bla,[4 1 2 3]);
    % chemy = compute_chemshift_vector_1_2(D1.Par.LarmorFreq,D1.Par.Dwelltimes(1)/10^9,D1.Par.vecSize);
    % figure; plot(chemy,bla(:,9,5))
    % figure; plot(chemy,bla(:,6,3))
    % figure; plot(chemy,bla(:,:))
    % figure; plot(chemy,sum(bla(:,:),2))





    %% Extract Phi from D1

    % for NuisLoop = 1:numel(Settings.EstSpecBases)
    %     if(Settings.EstSpecBases{1}.Flag)
    %         D1 = op_NuisanceRemoval_HSVD(D1,Settings.EstSpecBases{NuisLoop},mask,0);
    %     end
    % end

    
    AdditionalOut.USigmaPhi.WaterD1 = [];
    AdditionalOut.USigmaPhi.LipidD1 = [];
    AdditionalOut.USigmaPhi.MetaboD1 = [];
    
    % Water
    if(Settings.EstSpecBases.Water{1}.Flag)
        selParams.dt = D1.RecoPar.Dwelltimes(1)/10^9; selParams.NtEcho = 0; selParams.signalType = 'water';
        optWat.fSel = Settings.EstSpecBases.Water{1}.WaterfSel; 
        optWat.fSel2 = Settings.EstSpecBases.Water{1}.WaterfSel2;
        optWat.maxT2 = Settings.EstSpecBases.Water{1}.WaterMaxT2s;

        optLip.fSel2 = Settings.EstSpecBases.Water{1}.LipidfSel2; 
        optLip.maxT2 = Settings.EstSpecBases.Water{1}.LipidMaxT2s;    

        optMeta.fSel2 = Settings.EstSpecBases.Water{1}.MetabofSel2; 
        optMeta.minT2 = Settings.EstSpecBases.Water{1}.MetaboMinT2s;  

        temp = squeeze(D1.Data);
        temp                   = nsRm(temp, Settings.EstSpecBases.Water{1}.Mask, [], selParams, optWat, optLip, optMeta);
        temp = reshape(temp,[prod(size_MultiDims(temp,1:2)),size(temp,3)]);
        [UD1,SingValsD1,VD1] = svd(temp,'econ');

        R = Settings.EstSpecBases.Water{1}.NoOfSingVals;
        AdditionalOut.USigmaPhi.WaterD1.SingVals = SingValsD1(1:R,1:R);
        AdditionalOut.USigmaPhi.WaterD1.U        = UD1(:,1:R)*AdditionalOut.USigmaPhi.WaterD1.SingVals;
        AdditionalOut.USigmaPhi.WaterD1.Phi        = VD1(:,1:R)';
    end

    % Lipids
    if(Settings.EstSpecBases.Lipid{1}.Flag)
        selParams.dt = D1.RecoPar.Dwelltimes(1)/10^9; selParams.NtEcho = 0; selParams.signalType = 'lipid';
        optWat.fSel = Settings.EstSpecBases.Lipid{1}.WaterfSel; 
        optWat.fSel2 = Settings.EstSpecBases.Lipid{1}.WaterfSel2;
        optWat.maxT2 = Settings.EstSpecBases.Lipid{1}.WaterMaxT2s;

        optLip.fSel2 = Settings.EstSpecBases.Lipid{1}.LipidfSel2; 
        optLip.maxT2 = Settings.EstSpecBases.Lipid{1}.LipidMaxT2s;    

        optMeta.fSel2 = Settings.EstSpecBases.Lipid{1}.MetabofSel2; 
        optMeta.minT2 = Settings.EstSpecBases.Lipid{1}.MetaboMinT2s;  

        temp = squeeze(D1.Data);
        temp                   = nsRm(temp, Settings.EstSpecBases.Lipid{1}.Mask, [], selParams, optWat, optLip, optMeta);
        temp = reshape(temp,[prod(size_MultiDims(temp,1:2)),size(temp,3)]);
        [UD1,SingValsD1,VD1] = svd(temp,'econ');

        R = Settings.EstSpecBases.Lipid{1}.NoOfSingVals;
        AdditionalOut.USigmaPhi.LipidD1.SingVals = SingValsD1(1:R,1:R);
        AdditionalOut.USigmaPhi.LipidD1.U        = UD1(:,1:R)*AdditionalOut.USigmaPhi.LipidD1.SingVals;
        AdditionalOut.USigmaPhi.LipidD1.Phi        = VD1(:,1:R)';
    end

    % Metabo
    if(Settings.EstSpecBases.Metabo{1}.Flag)
        selParams.dt = D1.RecoPar.Dwelltimes(1)/10^9; selParams.NtEcho = 0; selParams.signalType = 'meta';
        optWat.fSel = Settings.EstSpecBases.Metabo{1}.WaterfSel; 
        optWat.fSel2 = Settings.EstSpecBases.Metabo{1}.WaterfSel2;
        optWat.maxT2 = Settings.EstSpecBases.Metabo{1}.WaterMaxT2s;

        optLip.fSel2 = Settings.EstSpecBases.Metabo{1}.LipidfSel2; 
        optLip.maxT2 = Settings.EstSpecBases.Metabo{1}.LipidMaxT2s;    

        optMeta.fSel2 = Settings.EstSpecBases.Metabo{1}.MetabofSel2; 
        optMeta.minT2 = Settings.EstSpecBases.Metabo{1}.MetaboMinT2s;  

        temp = squeeze(D1.Data);
        temp                   = nsRm(temp, Settings.EstSpecBases.Metabo{1}.Mask, [], selParams, optWat, optLip, optMeta);
        temp = reshape(temp,[prod(size_MultiDims(temp,1:2)),size(temp,3)]);
        [UD1,SingValsD1,VD1] = svd(temp,'econ');

        R = Settings.EstSpecBases.Metabo{1}.NoOfSingVals;
        AdditionalOut.USigmaPhi.MetaboD1.SingVals = SingValsD1(1:R,1:R);
        AdditionalOut.USigmaPhi.MetaboD1.U        = UD1(:,1:R)*AdditionalOut.USigmaPhi.MetaboD1.SingVals;
        AdditionalOut.USigmaPhi.MetaboD1.Phi        = VD1(:,1:R)';
    end

    % Interpolate Phis
    AdditionalOut.USigmaPhi.WaterD1 = interpV_bstr(AdditionalOut.USigmaPhi.WaterD1,Settings.EstSpecBases.TimeVecD1,Settings.EstSpecBases.TimeVecD2,false);
    AdditionalOut.USigmaPhi.LipidD1 = interpV_bstr(AdditionalOut.USigmaPhi.LipidD1,Settings.EstSpecBases.TimeVecD1,Settings.EstSpecBases.TimeVecD2,false);
    AdditionalOut.USigmaPhi.MetaboD1 = interpV_bstr(AdditionalOut.USigmaPhi.MetaboD1,Settings.EstSpecBases.TimeVecD1,Settings.EstSpecBases.TimeVecD2,false);
    

    %% Mask D1 and D2 Spectrally to Metabo Region
    % D1.Data = fftshift(fft(D1.Data,[],4),4);
    % D1.Data(:,:,:,[1:184,286:318],:) = 0;
    % D1.Data = ifft(ifftshift(D1.Data,4),[],4);
    % % D1.NoiseData = fftshift(fft(D1.NoiseData,[],4),4);
    % % D1.NoiseData(:,:,:,[1:184,286:318],:) = 0;
    % % D1.NoiseData = ifft(ifftshift(D1.NoiseData,4),[],4);
    % 
    % D2.Data = fftshift(fft(D2.Data,[],5),5);
    % D2.Data(:,:,:,:,[1:184,286:318],:) = 0;
    % D2.Data = ifft(ifftshift(D2.Data,5),[],5);
    % D2.NoiseData = fftshift(fft(D2.NoiseData,[],5),5);
    % D2.NoiseData(:,:,:,:,[1:184,286:318],:) = 0;
    % D2.NoiseData = ifft(ifftshift(D2.NoiseData,5),[],5);


    %% DEBUG: Check NuisRem
    % bla = D1.Data(:,:,1,:,1);
    % bla = bla .* myrepmat(D1.Mask,size(bla));
    % bla = abs(fftshift(fft(bla,[],4),4));
    % bla = permute(bla,[4 1 2 3]);
    % chemy = compute_chemshift_vector_1_2(D1.Par.LarmorFreq,D1.Par.Dwelltimes(1)/10^9,D1.Par.vecSize);
    % figure; plot(chemy,bla(:,:))
    % figure; plot(chemy,sum(bla(:,:),2))



    %% Extract Phi (Spectral Basis) from D1


    % [bla0,bla] = op_LowRankDenoiseMRSI(D1,[],Settings.Denoise_D1);
    % for CurCha = 1:size(bla.V,3)
    %     AdditionalOut.Phi{CurCha} = bla.V(:,1:bla.RankL(CurCha),CurCha)';
    % %         AdditionalOut.Phi{CurCha} = bla.V(:,1:25,CurCha)';    % REMOVE ME LATER
    % %     AdditionalOut.Phi{CurCha} = SpiceOperators.Phi(:,1:25,CurCha)';        % REMOVE ME LATER
    % %     AdditionalOut.RankL(CurCha) = 25;
    % end
    % AdditionalOut.RankL = bla.RankL;
    % clear bla



    %% Remove NuisanceSignal of D2

%     if(Settings.NuisRem_D2.NuisRem_flag)
%        IDontKnowYetHowToDoThis = true; 
%     end



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


    %% Reconstruct Spatial Map UTS (= U Times Sigma)


    tic_PCG = tic;
    fprintf('\n\nStart pcg . . .')

    % Define Input Data
    % The ratio of the FoV-overgrids is necessary, because the data is scaled by the matrix size when going Cartesian --> Spiral, but
    % the Cartesian data is zerofilled in image domain, and those zeros dont contribute any signal...
    % CartSize = SpSpice.GroundTruth_ForSpSpice.Par.DataSize(1);
    % sft2_Operator_ForSpSpice = sft2_Operator(transpose(squeeze(SpSpice.GroundTruth_ForSpSpice.TrajIdeal.GM)*CartSize),transpose(nsamp),1);
    % 
    % if(Ctrl.RecoPar.Correct4B0_flag)
    %     sft2_Operator_ForSpSpice = sft2_Operator_ForSpSpice .* SpSpice.Reco.B0CorrMat_Spatial;
    % end
    % 
    % % B0CorrMat = ones([prod(SpSpice.Reco.Par.DataSize(1:end-1)) SpSpice.Reco.Par.DataSize(end)]);
    % B0CorrMat = conj(reshape(SpSpice.Reco.B0CorrMat_Spec,[prod(SpSpice.Reco.Par.DataSize(1:end-1)) SpSpice.Reco.Par.DataSize(end)]));


    for CurAdditInd = 1:size(D2.Data,6)

        CurD2 = D2.Data(:,:,:,:,:,CurAdditInd);
        
        % Zerofill D2-data
%         CurD2(SpiceOperators.SamplingInd) = D2.Data(:,:,:,:,:,CurAdditInd);

        CurD2 = reshape(CurD2, [numel(CurD2) 1]);
        CurD2 = ModelFunction('Transj',CurD2,AdditionalOut.USigmaPhi,SpiceOperators);

%         funh = @(x) ModelFunction('Transj', ModelFunction('NoTransj',x,...
%         AdditionalOut.USigmaPhi,SpiceOperators),      ...
%         AdditionalOut.USigmaPhi,SpiceOperators);
    
        OptFun = @(x) OptimizationFunction(x, ModelFunction, AdditionalOut.USigmaPhi,SpiceOperators);

        InitGuess = [];

        AdditionalOut.UTS{CurAdditInd} = pcg(OptFun,CurD2,Settings.pcg.Tolerance,Settings.pcg.Iterations,[],[],InitGuess);
        AdditionalOut.UTS{CurAdditInd} = reshape(AdditionalOut.UTS{CurAdditInd},[size(SpiceOperators.B0CorrMat_Spec,1) AdditionalOut.RankL(CurAdditInd)]);

    %     D2.SamplingOperator = reshape(D2.SamplingOperator,D2.Par.DataSize);

    end



    %% Resynthesize Data Using Estimated UTS

    SpiceOut.Data = zeros(Settings.ResynthesizeSVDData.OutSize);
    for CurAdditInd = 1:numel(AdditionalOut.Phi)

        Temp = AdditionalOut.UTS{CurAdditInd} * AdditionalOut.Phi{CurAdditInd};
        SpiceOut.Data(:,:,:,:,CurAdditInd) = reshape(Temp,[Settings.ResynthesizeSVDData.OutSize(1:3) Settings.ResynthesizeSVDData.OutSize(4)]);

    end

    fprintf('\ntook %f s',toc(tic_PCG))

end

%% Subfunction: Interpolation
function SVDMatOut             = interpV_bstr(SVDMatIn,tD1,tD2,svd_on)
    if (~isempty(SVDMatIn) && ~isempty(SVDMatIn.Phi))
        if (length(tD1)~=length(tD2)) || norm(tD1-tD2) ~= 0 % if tD1 ~= tD2, do interpolation
            R                      = size(SVDMatIn.Phi,1);
            SVDMatOut.Phi                    = zeros(R,length(tD2));       
            for ind = 1:R
                SVDMatOut.Phi(ind,:)         = complex(interp1(tD1,real(squeeze(SVDMatIn.Phi(ind,:))),tD2,'spline'),...
                                             interp1(tD1,imag(squeeze(SVDMatIn.Phi(ind,:))),tD2,'spline'));
            end      
            if svd_on % do svd again
                temp               = SVDMatIn.U*SVDMatOut.Phi;
                [SVDMatOut.U,SVDMatOut.SingVals,tVD2]     = svd(temp,'econ');
                SVDMatOut.Phi                = tVD2;
                SVDMatOut.U                = SVDMatOut.U*SVDMatOut.SingVals;
            else
                SVDMatOut.SingVals                = SVDMatIn.SingVals;
                SVDMatOut.U                = SVDMatIn.U;
            end
        else % otherwise, no extra operation
            SVDMatOut.Phi                    = SVDMatIn.Phi;
            SVDMatOut.U                    = SVDMatIn.U;
            SVDMatOut.SingVals                    = SVDMatIn.SingVals;
        end
    else
        SVDMatOut.Phi                        = [];
        SVDMatOut.U                        = [];
        SVDMatOut.SingVals                        = [];
    end       
end    



