%%

clearvars
addpath(genpath('./dependencies'))





%% Definitions

SaveFigs.Flag = false;
SaveFigs.Path = './Figs/';
if(SaveFigs.Flag && ~exist(SaveFigs.Path,'dir'))
    mkdir(SaveFigs.Path)
end

Par.LarmorFreq = 297.223*10^6;  % 7T
Par.SBW = 926*3;
Par.Dwelltimes = 3.6e5;
Par.Nx = 64; Par.Ny = 64; Par.Nz = 1; Par.vecSize = 1024;
Par.DataSize = [Par.Nx Par.Ny Par.Nz Par.vecSize];

ParD1 = Par;
ParD1.Nx = 22;  ParD1.Ny = 22; ParD1.Nz = 1; ParD1.vecSize = 1024;
ParD1.DataSize = [ParD1.Nx ParD1.Ny ParD1.Nz ParD1.vecSize];

TimeUndersamplFactor = 2;

CalcMapsSetts.NAA.SumFromToPPM = [2.01-0.1 2.01+0.1];
CalcMapsSetts.Cho.SumFromToPPM = [3.21-0.1 3.21+0.1];
CalcMapsSetts.Cr.SumFromToPPM = [3.03-0.1 3.03+0.1];

PlotSpecs = {[33 33 1], [18 36 1], [39 47 1]};

DenoisingExample_flag = true;
SpiceExample_flag = true;


%% Define Spectral Bases, Metabolic Maps, InData.Maps.Masks etc

load('InData.mat');

% Take only some Metabolites
TakeOnlyMetabos = [20 21 23];

% Misuse metabolic maps for the following metabolites to use as water residual maps of different components:
% Glc, Asp, GSH, Gln, Ins, MM_mea, Scyllo-Inositol, Glu, Ins+Gly, Glx
% Remark: Using random water maps doesnt work properly, because when I calculate the low-res D1 dataset, it will smear all those water-components, and the
% smearing will be very strong when using noise (the low-res representation of noise is not accurate in comparison to the original high-res noise).
% This will cause problems when removing the water. That's why I changed now to using metabolite maps for water-maps, where the low-res representation
% is much more truthful to the original.
WaterMaps = reshape(InData.Maps.Metabos(:,:,1,[1 2 6 7 8 9 14 17 22 24]),[Par.Nx Par.Ny 1 10]);
% WaterMaps = randn([size(InData.Maps.Mask) 1 10]) .* reshape(1:-0.1:0.1,[1 1 1 10]);

InData.Maps.Metabos = InData.Maps.Metabos(:,:,1,TakeOnlyMetabos);

% Convert B0-map ppm --> Hz
InData.Maps.B0 = InData.Maps.B0 * Par.LarmorFreq / 10^6;
% Also overwrite two B0-values which are outliers
InData.Maps.B0(42,18) = 9.808; InData.Maps.B0(43,18) = 11.591;

Tmp = Simulate_FID_Spectra(3.21,4.65,0,0,0.03,0.1*9,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Cho 
InData.SpecComp(1,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(2.01,4.65,0,0,0.03,0.1*3,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % NAA
InData.SpecComp(2,:) = Tmp(2,:);
[Tmp,ppmvec] = Simulate_FID_Spectra(3.03,4.65,0,0,0.03,0.1*3,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);   % Cr 
InData.SpecComp(3,:) = Tmp(2,:);

% Residual Waters 1-10
Tmp = Simulate_FID_Spectra(4.65,4.65,10,0,0.03,4.0,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 1
InData.SpecComp(4,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.69,4.65,14,0,0.05,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 2
InData.SpecComp(5,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.61,4.65,-20,0,0.036,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 3
InData.SpecComp(6,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.63,4.65,-46,0,0.01,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 4
InData.SpecComp(7,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.88,4.65,90,0,0.02,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 5
InData.SpecComp(8,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.4,4.65,-10,0,0.015,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 6
InData.SpecComp(9,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.59,4.65,16,0,0.03,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 7
InData.SpecComp(10,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.73,4.65,21,0,0.03,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 8
InData.SpecComp(11,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.75,4.65,4,0,0.03,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 9
InData.SpecComp(12,:) = Tmp(2,:);
Tmp = Simulate_FID_Spectra(4.57,4.65,-3,0,0.03,4,0,1/Par.SBW,Par.vecSize,Par.LarmorFreq);            % Residual Water 10
InData.SpecComp(13,:) = Tmp(2,:);

InData.Maps.Metabos = cat(4,InData.Maps.Metabos,WaterMaps); clear WaterMaps


ppmvec = ppmvec(1,:);

clear AllMaps Tmp;

B0Scale = 1; % How strong the B0 should be simulated. If 1, it's the (measured) LCModel-shift map


if(DenoisingExample_flag)
    %% Rank of Gaussian Noise
    
    
    S = randn([5*5 Par.vecSize]) + 1i*randn([5*5 Par.vecSize]);

    fprintf('\nGaussian noise has rank: %d',rank(S))
    Results.Sigma = svd(S);
    GaussSingValFig = figure; plot(Results.Sigma,'LineWidth',2), title('Singular Values')
    Gauss_MetaboSpecFig = figure;
    subplot(1,2,1),imagesc(real(reshape(S(:,1),[5 5]))); 
    subplot(1,2,2),plot(real(S(1,:))); 


    if(SaveFigs.Flag)
        saveas(GaussSingValFig,[SaveFigs.Path '/GaussianNoise_SingVals.png']);
        saveas(Gauss_MetaboSpecFig,[SaveFigs.Path '/GaussianNoise_MetMapsSpecs.png']);

    end



    %% Only NAA Signal



    PickComponents = 2;
    NoOfMetabos = numel(PickComponents);

    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize);


    Results.NAAMap = CalcMetaboMap(S,InData.Maps.Mask,ppmvec,CalcMapsSetts.NAA);


    % Plot rank and singular values
    fprintf('\nMRSI data with only NAA has rank: %d',rank(reshape(S,[Par.Nx*Par.Ny Par.vecSize])))
    Results.Sigma = svd(reshape(S,[Par.Nx*Par.Ny Par.vecSize]));
    NAASingValsFig = figure; plot(Results.Sigma(1:30),'LineWidth',2), title('Singular Values')


    % Plot metabolic maps and spectra
    NAAMetMapSpecFig = figure; 
    subplot(NoOfMetabos,2,1)
    imagesc(Results.NAAMap); xticks([]); yticks([]); axis square; title('NAA, RMSE = 0')
    subplot(NoOfMetabos,2,2)
    plot(ppmvec,real(squeeze(fftshift(fft(S(PlotSpecs{1}(1),PlotSpecs{1}(2),1,:),[],4),4))),'LineWidth',2)


    if(SaveFigs.Flag)
        saveas(NAASingValsFig,[SaveFigs.Path '/NAA_SingVals.png']);
        saveas(NAAMetMapSpecFig,[SaveFigs.Path '/NAA_MetMapsSpecs.png']);

    end


    clear PickComponents NoOfMetabos CurMap S Results



    %% NAA, Cr, Cho


    PickComponents = 1:3;
    NoOfMetabos = numel(PickComponents);

    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize);


    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec);
    SaveFigs.Name = 'NAACrCho_GroundTruth';
    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs)
    Results_GroundTruth = Results;

    
   %% NAA, Cr, Cho, Enforce rank < 3

    PickComponents = 1:3;
    NoOfMetabos = numel(PickComponents);
    AssumedRank = 1;

    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize);

    % Denoise
    [U,Sigma,V] = svd(reshape(S,[Par.Nx*Par.Ny Par.vecSize]),'econ');
    S = U(:,1:AssumedRank) * Sigma(1:AssumedRank,1:AssumedRank) * V(:,1:AssumedRank)';
    S = reshape(S,Par.DataSize); clear CurMap

    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = 'NAACrCho_EnforceRank2';
    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs)

    

    %% Projection onto singular vectors vs spectral components
    % Investigate why the SVD takes for all singular vectors spectral peaks from the whole spectral range, and not 
    % SingularVector1st = NAA, SingularVector2nd = Cho, SingularVector3rd = Cr (or similar). The singular value decomposition tries
    % to decompose the given signals in such singular vectors, that explains in decreasing order most of the ground truth signal.
    % So lets compare how many singular (spectral) vectors we need to explain our data when projecting onto the singular vectors vs 
    % when projecting onto a basis consisting of (NAA,Cr,Cho).
    % Let's use a rank 1 model (all metabolite maps are that of NAA) to test that.
    % Remark: For some reason, my projections on the individual components are giving weird results when projecting on all three metabolites.
    % I guess I would need to account for the covariances between Cr and Cho (the projection on Cr also gives a bit of Cho and vice versa bc they are close)
    
    PickComponents = [1 2 3];
    NoOfMetabos = numel(PickComponents);

    CurMap = InData.Maps.Metabos(:,:,:,ones([1 NoOfMetabos])); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize);

    
    % Compare how well we can reproduce one voxel
    [U,Sigma,V] = svd(reshape(S,[Par.Nx*Par.Ny Par.vecSize]),'econ');
    V_DiffPeaks = InData.SpecComp(PickComponents,:)';
    V_DiffPeaks = V_DiffPeaks ./ transpose(sqrt(diag(V_DiffPeaks' * V_DiffPeaks)));
    
    
    % Calculate projections:
    S_ProjSingVec = reshape( reshape(S,[64*64 1024]) * V(:,1) * (V(:,1))',[64 64 1 1024]);
    S_ProjPeaks1 = reshape( reshape(S,[64*64 1024]) * V_DiffPeaks(:,1) * (V_DiffPeaks(:,1))',[64 64 1 1024]);
    S_ProjPeaks12 = reshape( reshape(S,[64*64 1024]) * V_DiffPeaks(:,1:2) * (V_DiffPeaks(:,1:2))',[64 64 1 1024]);
    S_ProjPeaks123 = reshape( reshape(S,[64*64 1024]) * V_DiffPeaks(:,1:3) * (V_DiffPeaks(:,1:3))',[64 64 1 1024]);
    
    % Add a bit of noise to S
    S = S + 0.01*randn(size(S)) +1i*0.01*randn(size(S));
    
    % Plot voxel 32 32 1, together with the projections
    figure; 
    subplot(2,2,1); plot(ppmvec,real(squeeze(fftshift(fft(S(32,32,1,:)))))); hold on;
    plot(ppmvec,real(squeeze(fftshift(fft(S_ProjSingVec(32,32,1,:))))),'r'); hold off;                   % Projection onto singular vectors
    title('Proj on SingVector1')
    
    
    subplot(2,2,2); plot(ppmvec,real(squeeze(fftshift(fft(S(32,32,1,:)))))); hold on;
    plot(ppmvec,real(squeeze(fftshift(fft(S_ProjPeaks1(32,32,1,:))))),'r'); hold off;                   % Projection onto singular vectors
    title('Proj on SpecComp1')
    
    subplot(2,2,3); plot(ppmvec,real(squeeze(fftshift(fft(S(32,32,1,:)))))); hold on;
    plot(ppmvec,real(squeeze(fftshift(fft(S_ProjPeaks12(32,32,1,:))))),'r'); hold off;                   % Projection onto singular vectors
    title('Proj on SpecComp1+2')
    
    subplot(2,2,4); plot(ppmvec,real(squeeze(fftshift(fft(S(32,32,1,:)))))); hold on;
    plot(ppmvec,real(squeeze(fftshift(fft(S_ProjPeaks123(32,32,1,:))))),'r'); hold off;                   % Projection onto singular vectors
    title('Proj on SpecComp1+2+3')
    
    
    %% Another explanation 
    % why the SVD does not just give the NAA, Cr, Cho signals:
    % This should simply show how you can decompose a rank-1 3x3 matrix by either using the singular vector ([1 2 3]) to explain all rows,
    % or by using three "basis vectors" (which could be e.g. the signals of Cr, Cho, NAA). When we chose the latter, we need three basis vectors,
    % while for the singular vector we need only 1. So by chosing the singular vector as basis, we are more efficient and need a lower number of vectors
    % to explain the same thing.
   
    
    Test = [1 2 3; 2 4 6; 3 6 9];
    BasisVec1 = [1 0 0];
    BasisVec2 = [0 2 0];
    BasisVec3 = [0 0 3];
    BasisVec4 = [1 2 3];
    
    Test
    Test_Rep2 = [BasisVec4; 2*BasisVec4; 3*BasisVec4]
    Test_Rep1 = [BasisVec1+BasisVec2+BasisVec3; 2*(BasisVec1+BasisVec2+BasisVec3); 3*(BasisVec1+BasisVec2+BasisVec3)]
    

    %% NAA, Cr, Cho + B0-Inhomogeneity
    % What happens to the rank if we add B0-inhomogeneity?

    PickComponents = 1:3;
    NoOfMetabos = numel(PickComponents);

    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize); clear CurMap

    % B0-Effect
    time   = (0:Par.vecSize-1)*(1/Par.SBW);
    B0CorrMat_Spec = exp(2*pi*1i*InData.Maps.B0*B0Scale .* reshape(time,[1 1 1 numel(time)]));
    S = S .* B0CorrMat_Spec;
    clear B0CorrMat_Spec time;


    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = 'NAACrCho_B0';
    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs)

    figure; imagesc(InData.Maps.B0*B0Scale);
    
    
    %% NAA, Cr, Cho + B0-Inhomogeneity. Enforce rank = 3
    % What happens if we enforce rank = 3 in case of B0-inhomogeneity?

    PickComponents = 1:3;
    NoOfMetabos = numel(PickComponents);
    AssumedRank = NoOfMetabos;
    B0Scale = 1.5; % How strong the B0 should be simulated. If 1, it's the (measured) LCMode-shift map
    
    
    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize); clear CurMap

    % B0-Effect
    time   = (0:Par.vecSize-1)*(1/Par.SBW);
    B0CorrMat_Spec = exp(2*pi*1i*InData.Maps.B0*B0Scale .* reshape(time,[1 1 1 numel(time)]));
    S = S .* B0CorrMat_Spec;
    clear B0CorrMat_Spec time;

    
    % Denoise
    [U,Sigma,V] = svd(reshape(S,[Par.Nx*Par.Ny Par.vecSize]),'econ');
    S = U(:,1:AssumedRank) * Sigma(1:AssumedRank,1:AssumedRank) * V(:,1:AssumedRank)';
    S = reshape(S,Par.DataSize); clear CurMap
    

    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = 'NAACrCho_B0_EnforceRank3';
    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs)

    figure; imagesc(InData.Maps.B0*B0Scale); colorbar, title('B0 map')




    %% NAA, Cr, Cho + residual water
    % What happens if we add some water residuals? We use 10 water components for that, with noise-like metabolic maps
    
    PickComponents = 1:13;
    NoOfMetabos = numel(PickComponents);

    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize); clear CurMap

    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = 'NAACrCho_ResWater';
    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs)

    
    %% NAA, Cr, Cho + residual water. Enforce rank = 3
    % What happens if we enforce rank 3 on this dataset with residual water? The metabolites will be mostly gone, bc the water signal dominates the svd.

    PickComponents = 1:13;
    NoOfMetabos = numel(PickComponents);
    AssumedRank = 3;

    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize); clear CurMap

    % Denoise
    [U,Sigma,V] = svd(reshape(S,[Par.Nx*Par.Ny Par.vecSize]),'econ');
    S = U(:,1:AssumedRank) * Sigma(1:AssumedRank,1:AssumedRank) * V(:,1:AssumedRank)';
    S = reshape(S,Par.DataSize); clear CurMap
    
    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = 'NAACrCho_ResWater_EnforceRank3';
    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs)
    


    %% NAA, Cr, Cho + Noise
    % What happens with the rank if we add noise to our data?
    
    PickComponents = 1:3;
    NoOfMetabos = numel(PickComponents);

    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]); 
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize); clear CurMap
    S_GroundTruth = S;


    % Add Noise
    NoiseScale = 0.25; 
    S = S + NoiseScale*randn(size(S)) + NoiseScale*1i*randn(size(S));
    S_NoisyFullRank = S;

    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = ['NAACrCho_Noise_Scale' num2str(NoiseScale)];
    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs,'Noisy Data')



    % NAA, Cr, Cho + Noise, Denoised
    % Use data from above

    AssumedRank = NoOfMetabos;


    % Denoise
    [U,Sigma,V] = svd(reshape(S,[Par.Nx*Par.Ny Par.vecSize]),'econ');
    S = U(:,1:AssumedRank) * Sigma(1:AssumedRank,1:AssumedRank) * V(:,1:AssumedRank)';
    S = reshape(S,Par.DataSize); clear CurMap

    Results = CalcResults(S,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = ['NAACrCho_Noise_Scale' num2str(NoiseScale) '_Denoised'];

    ShowResults(Results,PlotSpecs,ppmvec,S,SaveFigs,'Denoised Data')

    
    
    %% How Does the Spectral Residuals Look Like?
    % Is the noise evenly distributed in the spectrum?
    
    ShowResults(Results,PlotSpecs,ppmvec,S-S_GroundTruth,SaveFigs,'Denoised - GroundTruth')
    ShowResults(Results,PlotSpecs,ppmvec,S-S_NoisyFullRank,SaveFigs,'Denoised - Noisy')    
    

end



if(SpiceExample_flag)

    %% Calculate GroundTruth & Definitions
    
    B0Scale = 1;
    
    PickComponents = 1:3;
    NoOfMetabos = numel(PickComponents);
    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S_GroundTruth = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize);
    Results_GroundTruth = CalcResults(S_GroundTruth,InData.Maps.Mask,Par,CalcMapsSetts,ppmvec);
    SaveFigs.Name = 'NAACrCho_GroundTruth';
    
    
    %% Spice Reco: Create D1 & D2 Data
    PickComponents = 1:13;
    NoOfMetabos = numel(PickComponents);
    CurMap = InData.Maps.Metabos(:,:,:,PickComponents); CurMap = reshape(CurMap,[numel(CurMap)/NoOfMetabos NoOfMetabos]);
    S = reshape(CurMap * InData.SpecComp(PickComponents,:),Par.DataSize); clear CurMap

    % B0-Effect
    time   = (0:Par.vecSize-1)*(1/Par.SBW);
    B0CorrMat_Spec = exp(2*pi*1i*InData.Maps.B0*B0Scale .* reshape(time,[1 1 1 numel(time)]));
    S = S .* B0CorrMat_Spec;
    clear B0CorrMat_Spec time;

    % Add Noise
    NoiseScale = 0.25;
    S = S + NoiseScale*randn(size(S)) + NoiseScale*1i*randn(size(S));


    % Undersample data to size of D1
    D1.Data = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(conj(S),1),2),[],1),[],2),1),2);
    D1.Par = ParD1;
    D1.Data = ZerofillOrCutkSpace(D1.Data,D1.Par.DataSize,0);
    D1.Mask = imresize(InData.Maps.Mask,D1.Par.DataSize(1:2),'nearest');
    D1.BrainMask = D1.Mask;
    D1.LipMask = zeros(size(D1.Mask));

    % Undersample data to size of D2
    D2.Data = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(conj(S),1),2),[],1),[],2),1),2);
    D2.Par = Par;
    D2.Mask = InData.Maps.Mask;
    D2.BrainMask = D2.Mask;
    D2.LipMask = zeros(size(D2.Mask));
    SamplingOperator = zeros(D2.Par.DataSize); SamplingOperator(:,:,:,1:TimeUndersamplFactor:end) = 1;
    D2.Data = D2.Data .* SamplingOperator;

    D2_BefNuisRem = D2;
    
    % Show D1
    D1Reco = D1;
    D1Reco.Data = FFTOfMRIData(D1Reco.Data,1,[1 2],0,1,0);
    
    % D1Reco_zf
    D1Reco_zf = D1;
    D1Reco_zf.Data = ZerofillOrCutkSpace(D1.Data,D2.Par.DataSize,0);
    D1Reco_zf.Par = D2.Par;
    D1Reco_zf.Data = FFTOfMRIData(D1Reco_zf.Data,1,[1 2],0,1,0);
    


    %% Spice Reco: Remove nuisance signal from D2
    % Substeps: 
    %       Extract Water & Metabo Phis from D1
    %       Estimate gamma (spatial components) using above phi & undersampled D2 data
    %       Remove nuisance signal by calculating D2NuisRem = D2 - gamma_water*phi_water

    
    % A bit messy all those preparatory stuff, but that's how it is ;)

    
     % The functions used for nsRmSparse_bstrchanged have same name, but are different. This is a hack
    rmpath('./dependencies/op_NuisanceRemoval_HSVD')
    addpath('./dependencies/op_NuisanceRemoval_HSVD_Sparse')


    
%     SpiceOperators.B0 = InData.Maps.B0; 
%     SpiceOperators.SamplingOperator = SamplingOperator;

    MetaboR = 5;

    params.RM = MetaboR;
    params.RW = 15; 
    params.RL = 0;
    
    params.tD1 = (0:D1Reco.Par.Dwelltimes(1):D1Reco.Par.Dwelltimes(1)*(D1Reco.Par.vecSize-1))/10^9;
    params.tD2 = params.tD1;
                                                                   
                                                            
    params.ND1 = D1Reco.Par.DataSize([1 2 4]);
    params.ND2 = D2.Par.DataSize([1 2 4]);
    params.ImMask = D2.Mask;         % This needs to be brain + lipids mask, bc in nsRmSparseInitParams the lipmask is calculated lipmask = ImMask - waterMask;
    params.waterMask = D2.Mask; % This needs to be the brain-mask (NOT brain + lipids!) 
    params.B0Map = InData.Maps.B0*B0Scale;  
    params.NosD1 = 1;
    params.dt = D1Reco.Par.Dwelltimes(1)/10^9;
    params.lambdaDfWat = 0; params.lambdaDfLip = 0; params.lambdaDfMet = 0;
    params.csiInfo.NtEcho = 0;
    params = nsRmSparseInitParams(params);
    params.debug_on = 0;
    
    
    params.xtD1W = squeeze(D1Reco.Data);
    params.xtD1WrL = params.xtD1W;
    params.xtD1WrLr = params.xtD1W;
    params.waterMaskD2 = double(params.waterMaskD2);
    params.waterMaskD2(params.waterMaskD2 == 0) = 0.01;
    params.lipMaskD2(params.lipMaskD2 == 0) = 0.01;
    ConjFlag = 1;

    Factor = D1Reco.Par.LarmorFreq * 1e-6;
    ConjSign = (-1)^(double(ConjFlag));
    
    params.optWat.fSel = [80 -80];
    params.optWat.fSel2 = sort(ConjSign*([4.2,  5.2] - 4.65)*Factor);
    params.optWat.maxT2 = [1e6];
%     params.optLip.fSel2 = sort(ConjSign*(SettingsTemp.NuisRem_D1{1}.LipidPPMs([1:2:numel(SettingsTemp.NuisRem_D1{1}.LipidPPMs) numel(SettingsTemp.NuisRem_D1{1}.LipidPPMs)]) - 4.65)*Factor);
%     params.optLip.maxT2 = SettingsTemp.NuisRem_D1{1}.LipidT2s;
    params.optMeta.fSel2 = sort(ConjSign*([3.7575, 1.9] - 4.65)*Factor);
    %     params.optMeta.maxT2 = 25;


    params.sampleMask = false(D2.Par.DataSize);
    params.sampleMask(:,:,:,1:TimeUndersamplFactor:end) = true;
    if(numel(params.sampleMask) > 2^32-1)
        UseUintFun = @uint64; 
    else
        UseUintFun = @uint32;
    end
    params.sampleIndex = find(params.sampleMask);
    params.sampleIndexShift = params.sampleIndex;

    params.Operators.SensMap = 1;

    params.Operators.SamplingOperator = params.sampleMask;
    params.Operators.InDataSize = [D2.Par.DataSize(1:3) 1 ceil(D2.Par.DataSize(4)/TimeUndersamplFactor)];
    params.Operators.OutDataSize = D2.Par.DataSize;
    
    params.verbose = 0;


    params.MaxIterCG = 100;
    params.D1RecoPar = D1Reco.Par;
    params.ND2In = D2.Par.DataSize(1:2);
  
    params.ShowSpectra_flag = 1;
    params.ShowVoxels = {[32 32 1], [32 33 1], [33 32 1],[33 33 1]};
    params.tolCG = 1E-10;
    
    
    for CurCha = 1:1
        fprintf('\n\nDecontaminating Channel %d\n',CurCha)
        d = reshape(D2.Data(:,:,:,1:TimeUndersamplFactor:end,CurCha),params.Operators.InDataSize);
        [ dWL, dOut, nsRmSparseOutput, params ] = nsRmSparse_bstrchanged( d, params );
        D2.Data(:,:,:,1:TimeUndersamplFactor:end,CurCha) = dOut;
    end   
	clear params dWL d ConjFlag ConjSign CurCha Factor UseUintFun dOut

    
    rmpath('./dependencies/op_NuisanceRemoval_HSVD_Sparse')
    addpath('./dependencies/op_NuisanceRemoval_HSVD')
   
    
    
    %% Spice Reco: Perform Spice Reco
    
    
    % A bit messy all those preparatory stuff, but that's how it is ;)
    addpath('./dependencies/op_SpiceReco')

    % Reshape Data
    D2_bak = D2;
    SpiceOperators.InDataSize = [D2.Par.DataSize(1:3) 1 D2.Par.DataSize(4)];
    D2.Data = reshape(D2.Data,SpiceOperators.InDataSize);
    D2.Par.DataSize = SpiceOperators.InDataSize;
    
%     SpiceOperators.B0 = InData.Maps.B0; 
    SpiceOperators.SamplingOperator = reshape(SamplingOperator,SpiceOperators.InDataSize);
    SpiceOperators.SensMap = 1;
    SpiceOperators.FoVShift = 1;
    SpiceOperators.OutDataSize = D2.Par.DataSize([1:3 5]);
    SpiceOperators.Mask = D2.Mask;
    
    time   = (0:SpiceOperators.OutDataSize(4)-1)*D2.Par.Dwelltimes(1)/10^9; time = reshape(time,[1 1 1 numel(time)]);
    B0CorrMat_Spec = exp(-2*pi*1i*InData.Maps.B0 *B0Scale .* time);
    SpiceOperators.B0CorrMat_Spec = B0CorrMat_Spec;
    clear time B0CorrMat_Spec;
    

    ModelFunction = @CompSens_ModelData; % For Antoines function, because he does the projections on V manually



    % % % Set Up Settings for Spice % % %

    % CorrSpectralB0
    Settings.Spice.CorrSpectralB0_D1.Flag = true;
    Settings.Spice.CorrSpectralB0_D1.RoundB0ToIntVecPts = false;
    Settings.Spice.CorrSpectralB0_D1.Debug_flag = false;

    % Remove Nuisance Signal
    Settings.Spice.NuisRem_D1{1}.Flag = true;     
    Settings.Spice.NuisRem_D1{1}.NoOfSingVals = 15;
    Settings.Spice.NuisRem_D1{1}.WaterPPMs = [5.9, 4.2]; % Originally had: [8.0, 3.8386; 3.8386, 3.6764]. Results in same as original-nsrm: [5.463, 5.3; 5.3, 3.024]
    Settings.Spice.NuisRem_D1{1}.WaterT2s = [1e6]; % [1e6;20]
    Settings.Spice.NuisRem_D1{1}.LipidPPMs = [8.7069, 8.7069];
    Settings.Spice.NuisRem_D1{1}.LipidT2s = [-2000];
    Settings.Spice.NuisRem_D1{1}.MetaboPPMs = [4.2, 1.9];
    Settings.Spice.NuisRem_D1{1}.OtherPPMs = []; % 9.0,-1
    Settings.Spice.NuisRem_D1{1}.OtherT2s = []; %1e6
    Settings.Spice.NuisRem_D1{1}.Debug_flag = false;
    Settings.Spice.NuisRem_D1{1}.SignalType = 'water';
 
    
    Settings.Spice.FreqAlign2_D1 = struct('Flag',false,'ApplyAlongDim',4,'ZerofillingFactor',2,'PeakSearchPPM',2.01,'PeakSearchRangePPM',0.2,'AlignRefPeak_flag',1,'UseSVDForRef_flag',1);
 
    % Exp Filter
    Settings.Spice.ExpFilter_D1 = struct('Flag',false,'ApplyAlongDim',4,'ExpFilterStrength_Hz',-10);
    
    % Rank for estimated phi
    Settings.Spice.Denoise_D1.MaxRankL = MetaboR;
    Settings.Spice.Denoise_D1.MinRankL = MetaboR;
    Settings.Spice.Denoise_D1.Debug_flag = false;

    % NuisRemance of D2
    Settings.Spice.NuisRem_D2.NuisRem_flag = false;

    % IterReco
    Settings.Spice.IterReco.Iterations = 100;
    Settings.Spice.IterReco.Tolerance = 10^-6;

    % Resynthesize
    Settings.Spice.ResynthesizeSVDData.OutSize = SpiceOperators.OutDataSize;

    % Zerofill Input-Data       
    D2.Par.TimeUndersamplFactor = TimeUndersamplFactor;
    D2.Data = D2.Data .* SpiceOperators.SamplingOperator;

    % D1 Data
    D1ForSpice = D1Reco;


    % RECO: Perform SPICE Reconstruction
    TmpB0.B0Map = -InData.Maps.B0*B0Scale;
    [D2SpiceReco, SpiceAddOutput] = op_SpiceReco(D1ForSpice,D2,TmpB0,ModelFunction,SpiceOperators,Settings.Spice);

    D2SpiceReco.Mask = D2.Mask;
    D2SpiceReco.BrainMask = D2.BrainMask;
    D2SpiceReco.Par = D2.Par;
    D2SpiceReco.Par.DataSize = size(D2SpiceReco.Data);
    

    
    D2 = D2_bak;
    
    
    %% Spice Reco: Show results:
    
    % GroundTruth
    ShowResults(Results_GroundTruth,PlotSpecs,ppmvec,S_GroundTruth,false,'GroundTruth D2')

    
    Results = CalcResults(D1Reco_zf.Data,D2.Mask,D1Reco_zf.Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = ['NAACrCho_D1'];
    ShowResults(Results,PlotSpecs,ppmvec,D1Reco_zf.Data,SaveFigs,'GroundTruth D1')
    % Show D2BefNuisRem
    D2Reco = D2;
    D2Reco.Data = FFTOfMRIData(D2_BefNuisRem.Data,1,[1 2],0,1,0);    
    Results = CalcResults(D2Reco.Data,D2Reco.Mask,D2Reco.Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = ['NAACrCho_D2'];
    ShowResults(Results,PlotSpecs,ppmvec,D2Reco.Data,SaveFigs,'D2 BeforeReco')    

    clear D2Reco;
    
    
    % Metabo & Water phis (spectral components) for D2NuisRem
    H2OPhi1Fig = figure;
    for ii = 1:size(nsRmSparseOutput.VWD1,1)
        subplot(4,4,ii); plot(ppmvec,abs(fftshift(fft(nsRmSparseOutput.VWD1(ii,:))))); title(['H2OSpecComp' num2str(ii)])        
    end
    MetPhi1Fig = figure;
    SubPlotSize = ceil(sqrt(size(nsRmSparseOutput.VMD1,1)));
    for ii = 1:size(nsRmSparseOutput.VMD1,1)
        subplot(SubPlotSize,SubPlotSize,ii); plot(ppmvec,abs(fftshift(fft(nsRmSparseOutput.VMD1(ii,:))))); title(['MetSpecComp' num2str(ii)])        
    end

    % Metabo & Water gammas (spatial components)
    Tmp = reshape(nsRmSparseOutput.xSol,[D2.Par.DataSize(1) D2.Par.DataSize(1) size(nsRmSparseOutput.VWD1,1)+size(nsRmSparseOutput.VMD1,1)]) .* D2.Mask;
    H2OGamma1Fig = figure;
    SubPlotSize = ceil(sqrt(size(nsRmSparseOutput.VWD1,1)));    
    for ii = 1:size(nsRmSparseOutput.VWD1,1)
        subplot(SubPlotSize,SubPlotSize,ii); imagesc(abs(Tmp(:,:,ii))); title(['H2OSpatComp' num2str(ii)])
    end
    MetGamma1Fig = figure;
    SubPlotSize = ceil(sqrt(size(nsRmSparseOutput.VMD1,1)));
    for ii = (size(nsRmSparseOutput.VWD1,1)+1):size(nsRmSparseOutput.VWD1,1)+size(nsRmSparseOutput.VMD1,1)
        subplot(SubPlotSize,SubPlotSize,ii-size(nsRmSparseOutput.VWD1,1)); imagesc(abs(Tmp(:,:,ii))); title(['MetSpatComp' num2str(ii-size(nsRmSparseOutput.VWD1,1))])
    end    
  
    
    % Show D2 NuisRem
    D2NuisRemReco = D2; 
    D2NuisRemReco.Data = FFTOfMRIData(conj(D2.Data),1,[1 2],0,1,0);    
    Results = CalcResults(D2NuisRemReco.Data,D2NuisRemReco.Mask,D2NuisRemReco.Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = ['NAACrCho_D2NuisRem'];
    ShowResults(Results,PlotSpecs,ppmvec,D2NuisRemReco.Data,SaveFigs,'D2 NuisRem')    

    
    % Show D1NuisRem
    D1NuisRemReco = SpiceAddOutput.D1; 
    D1NuisRemReco.Data = ZerofillOrCutkSpace(D1NuisRemReco.Data,D2.Par.DataSize,1);
    D1NuisRemReco.Mask = D2.Mask; D1NuisRemReco.Par = D2.Par;
    Results = CalcResults(D1NuisRemReco.Data,D1NuisRemReco.Mask,D1NuisRemReco.Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = ['NAACrCho_D1NuisRem'];
    ShowResults(Results,PlotSpecs,ppmvec,D1NuisRemReco.Data,SaveFigs,'D1 NuisRem')    
            
    
    % Show Metabo phis from D1NuisRem
    MetPhi2Fig = figure;
    SubplotNo = ceil(sqrt(size(SpiceAddOutput.Phi{1},1)));
    for ii = 1:size(SpiceAddOutput.Phi{1},1)
        subplot(SubplotNo,SubplotNo,ii); plot(ppmvec,abs(fftshift(fft(SpiceAddOutput.Phi{1}(ii,:))))); title(['MetSpecComp' num2str(ii)])        
    end    
    
    
    % Show Metabo gamma
    MetGamma2Fig = figure;
    for ii = 1:size(SpiceAddOutput.UTS{1},2)
        subplot(SubplotNo,SubplotNo,ii); imagesc(reshape(abs(SpiceAddOutput.UTS{1}(:,ii)),D2.Par.DataSize(1:2))); title(['MetSpatComp' num2str(ii)])
    end        
    
  
    % Show SpiceReco
%     D2NuisRemReco.Data = FFTOfMRIData(D2.Data,1,[1 2],0,1,0);    
    Results = CalcResults(D2SpiceReco.Data,D2SpiceReco.Mask,D2SpiceReco.Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    SaveFigs.Name = ['NAACrCho_SpiceReco'];
    ShowResults(Results,PlotSpecs,ppmvec,D2SpiceReco.Data,SaveFigs,'Spice Reco')    
            
    %
    Results = CalcResults(D2SpiceReco.Data - S_GroundTruth,D2SpiceReco.Mask,D2SpiceReco.Par,CalcMapsSetts,ppmvec,Results_GroundTruth);
    ShowResults(Results,PlotSpecs,ppmvec,D2SpiceReco.Data,SaveFigs,'Spice Reco - GT')    
    
    
    if(SaveFigs.Flag)
        saveas(H2OPhi1Fig,[SaveFigs.Path '/H2O_Phi1.png']);
        saveas(MetPhi1Fig,[SaveFigs.Path '/Met_Phi1.png']);
        saveas(H2OGamma1Fig,[SaveFigs.Path '/H2O_Gamma1.png']);
        saveas(MetGamma1Fig,[SaveFigs.Path '/Met_Gamma1.png']);
        saveas(MetPhi2Fig,[SaveFigs.Path '/Met_Phi2.png']);
        saveas(MetGamma2Fig,[SaveFigs.Path '/Met_Gamma2.png']);

        
    end
    

end


%% Function for calculating metabolic maps from signal
function MetaboMap = CalcMetaboMap(Signal,Mask,ppmvec,Settings)

    Signal_Spec = fftshift(fft(Signal,[],4),4);
    Tmp = FindClosestIndex(ppmvec,Settings.SumFromToPPM(1)); SumFromTo(1) = Tmp{1};
    Tmp = FindClosestIndex(ppmvec,Settings.SumFromToPPM(2)); SumFromTo(2) = Tmp{1};
    SumFromTo = sort(SumFromTo);
    
    MetaboMap = sum(real(Signal_Spec(:,:,:,SumFromTo(1):SumFromTo(2))),4) .* Mask;

end



%% Function for Plotting Results
function ShowResults(Results,PlotSpecs, ppmvec,S,SaveFigs,SpecTitle)

    if(~isstruct(SaveFigs))
        Tmp = SaveFigs; clear SaveFigs; SaveFigs.Flag = Tmp; clear Tmp;
        SaveFigs.Name = 'NotSpecified';
    end
    % Plot rank and singular values
    SingValFig = figure; plot(Results.Sigma(1:30),'LineWidth',2), title('Singular Values'), ylim([0 310]);


    % Plot metabolic maps and spectra
    MetMapSpecFig = figure; 
    subplot(3,2,1)
    imagesc(Results.NAAMap); title(['NAAMap, RMSE = ' num2str(Results.RMSEs.NAA)]), xticks([]); yticks([]); %, axis square
    subplot(3,2,3)
    imagesc(Results.CrMap); title(['CrMap, RMSE = ' num2str(Results.RMSEs.Cr)]), xticks([]); yticks([]) %, axis square
    subplot(3,2,5)
    imagesc(Results.ChoMap); title(['ChoMap, RMSE = ' num2str(Results.RMSEs.Cho)]), xticks([]); yticks([]) %, axis square
    subplot(3,2,2)
    plot(ppmvec,real(squeeze(fftshift(fft(S(PlotSpecs{1}(1),PlotSpecs{1}(2),1,:),[],4),4))),'LineWidth',1.5), yticks([]), xlim([1 5.5])
    if(exist('SpecTitle','var') && ~isempty(SpecTitle))
        title(SpecTitle) 
    end
    subplot(3,2,4)
    plot(ppmvec,real(squeeze(fftshift(fft(S(PlotSpecs{2}(1),PlotSpecs{2}(2),1,:),[],4),4))),'LineWidth',1.5), yticks([]), xlim([1 5.5])
    subplot(3,2,6)
    plot(ppmvec,real(squeeze(fftshift(fft(S(PlotSpecs{3}(1),PlotSpecs{3}(2),1,:),[],4),4))),'LineWidth',1.5), yticks([]), xlim([1 5.5])
    
    fprintf('\nDataset %s has rank: %d',SaveFigs.Name,Results.Rank)
    
    if(SaveFigs.Flag)
        saveas(SingValFig,[SaveFigs.Path '/' SaveFigs.Name '_SingVals.png']);
        saveas(MetMapSpecFig,[SaveFigs.Path '/' SaveFigs.Name '_MetMapsSpecs.png']);
%         imwrite((SingValFig),[SaveFigs.Path '/' SaveFigs.Name '_SingVals.png'],'BitDepth',8);
%         imwrite((MetMapSpecFig),[SaveFigs.Path '/' SaveFigs.Name '_MetMapsSpecs.png'],'BitDepth',8);
    end
    

end



%% Function for Calculating Results
function Results = CalcResults(S,Mask,Par,CalcMapsSetts,ppmvec,Results_GroundTruth)

    Results.NAAMap = CalcMetaboMap(S,Mask,ppmvec,CalcMapsSetts.NAA);
    Results.CrMap = CalcMetaboMap(S,Mask,ppmvec,CalcMapsSetts.Cr);
    Results.ChoMap = CalcMetaboMap(S,Mask,ppmvec,CalcMapsSetts.Cho);
    Results.Sigma = svd(reshape(S,[Par.Nx*Par.Ny Par.vecSize]));
    Results.Rank = rank(reshape(S,[Par.Nx*Par.Ny Par.vecSize]));
    
    if(exist('Results_GroundTruth','var'))
        Results.RMSEs.NAA = sqrt(mean(abs(Results_GroundTruth.NAAMap(~Mask==0) - Results.NAAMap(~Mask==0)).^2));
        Results.RMSEs.Cr = sqrt(mean(abs(Results_GroundTruth.CrMap(~Mask==0) - Results.CrMap(~Mask==0)).^2));
        Results.RMSEs.Cho = sqrt(mean(abs(Results_GroundTruth.ChoMap(~Mask==0) - Results.ChoMap(~Mask==0)).^2));
    else
        Results.RMSEs.NAA = 0;
        Results.RMSEs.Cr = 0;
        Results.RMSEs.Cho = 0;
    end

end
