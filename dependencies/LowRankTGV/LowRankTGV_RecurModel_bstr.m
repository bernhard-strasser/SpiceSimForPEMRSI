function [ResultData_trr,U_rrrc,V_tc,S] = LowRankTGV_RecurModel_bstr(mrsiData_kkkstc,mrsiReconParams,A,Ah,Settings)
% *************************************************************************
% What is needed in mrsiReconParams:
% mrsiReconParams.modelOrder: Low rank order
% mrsiReconParams.BrainMask2D: 2D Brain Mask
% mrsiReconParams.WaterFreqMap: B0 field map in Herz
% mrsiReconParams.mrProt.samplerate: SampleRate (Acqu. BW) in Hz
% mrsiReconParams.SENSE: coil sensitivity profiles
% mrsiReconParams.mu_tv: TGV regularization paramter (5E-4 - 1E-3 for 2D)
% mrsiReconParams.kmask:  Mask of the aquired K-space points
% *************************************************************************
% Compute SVD of Adjoint Solution
%*************************************************************************

if(~exist('Settings','var'))
    Settings = struct();
end
if(~isfield(Settings,'WriteFiles_flag'))
    Settings.WriteFiles_flag = false;
end
if(~isfield(Settings,'WriteFiles_Path'))
    Settings.WriteFiles_Path = './LowRankTGVRecon_DebugOutput'; 
end
if(~isfield(Settings,'FixedV_flag'))
    Settings.FixedV_flag = false;
end
if(Settings.FixedV_flag)
	mrsiReconParams.LRTGVModelParams.SpecItFact = 0;
end

fprintf('Compute initial SVD Adjoint solution...\n');
SiOri=size(mrsiData_kkkstc);
nDimsOri= ndims(mrsiData_kkkstc);

%Assign useful variables
% NbT=SiOri(4);
% FreqMap=mrsiReconParams.WaterFreqMap; % B0map ;
% Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;

% Time_rrt=permute(repmat(([0 :(NbT-1)]'/Fs),[1 SiOri(2) SiOri(3)]),[2,3,1]);               % bstr: I precalculated the exp(2pi*Time*Freq), no need to do it again
% Freqshift_rrt=exp(-2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT])); %-2pi to go from k -> r
% Freqshift_rrt = squeeze(FreqMap);
% Freqshift_crrt= permute(repmat(Freqshift_rrt,[1 1 1 SiOri(1)]),[4 1 2 3]);
% clear  Time_rrt Freqshift_rrt

% Brainmask_crrt=permute(repmat(mrsiReconParams.BrainMask2D,[1 1 size(mrsiData_kkkstc,1) size(mrsiData_kkkstc,4)]),[3 1 2 4]);
% SENSE_crrt=repmat(mrsiReconParams.SENSE,[1 1 1 NbT]);     % bstr: Did this repmat alrdy before

AdjointData_rrrt = Ah(mrsiData_kkkstc);
SiOut = size(AdjointData_rrrt);
% AdjointData_rrrt =conj(SENSE_crrt).*ifft(ifft(mrsiData_kkkstc,[],2),[],3).*Freqshift_crrt.*Brainmask_crrt; % coil-r-r-t
% AdjointData_rrrt = squeeze(sum(AdjointData_rrrt,1));% r-r-t

clear Brainmask_crrt Freqshift_crrt SENSE_crrt

[Uorig,Sorig,Vorig] = svd(reshape(AdjointData_rrrt,[],SiOut(4)),0);
if(Settings.FixedV_flag)
    V_tc = Settings.V;
else
    V_tc=Vorig(:,1:mrsiReconParams.modelOrder);
end
S=Sorig(1:mrsiReconParams.modelOrder,1:mrsiReconParams.modelOrder);
U_rrrc=reshape(Uorig(:,1:mrsiReconParams.modelOrder), SiOut(1),SiOut(2),SiOut(3),[]);
Init_U_rrrc = U_rrrc;

% Print out the SVD spectra in a eps file
if(Settings.WriteFiles_flag)
    if ~exist(Settings.WriteFiles_Path,'dir')
        mkdir(Settings.WriteFiles_Path)
    end
    VisualizeSpectral( V_tc,S, [Settings.WriteFiles_Path '/Initial_Spectral_Components']);
end

clear  Uorig Vorig 

%%%%%%%%%%%
%Test Random Initialization (Doesn't Converge so far)
%{
U_rrrc=rand(size(U_rrrc));
V_tc=2*(rand(size(V_tc))+1j*rand(size(V_tc))) -1-1j;
for c=1:size(V_tc,2)
	V_tc=V_tc/norm(V_tc(:,c));
	U_rrrc=U_rrrc/norm(U_rrrc(:,:,c));
end
S=eye(size(S))*S(1,1)/10;
%}

SizeVol = size(U_rrrc);


% DOES NOT EXIST!!!
% fprintf(makeSectionDisplayHeader('Start the recursive reconstruction...\n'));

%% TGV Parameters

alpha = mrsiReconParams.mu_tv;
maxit = mrsiReconParams.maxit;  %Maximum nb of iterations. 1200 should be more than safe 
minit = mrsiReconParams.minit;  %Minimum. It won't let stop before that
Threshold = mrsiReconParams.Threshold;  % Will stop in threshold is reached and if the nb of iteration > minit


[U_rrrc,V_tc,S, costFunVal] = tgv2_l2_2D_multiCoil_LowRank_CombinedConv_bstr(mrsiData_kkkstc,A,Ah,U_rrrc,V_tc,S, 2*alpha, alpha, maxit,minit,mrsiReconParams,Threshold,Settings);

% Print Finale Components in eps files
if(Settings.WriteFiles_flag)
    VisualizeTGV( Init_U_rrrc ,U_rrrc,[Settings.WriteFiles_Path '/LowRank-TGV_Resulting_Spatial_Components']);
    VisualizeSpectral( V_tc,S, [Settings.WriteFiles_Path '/LowRank-TGV_Resulting_Spectral_Components'])
    s=[ './LowRank-TGV_Final_FrequencyMap.ps'];
    % figs=figure('visible','off');imagesc(DynamicFreqMap,[-30 30]);
    % colorbar; title('Final Adaptative Frequency Map');print(figs, '-dpsc2',s);
end
    
fprintf('Recombination of the Data...\n');
ResultData_trr=formTensorProduct(U_rrrc,V_tc*S,nDimsOri-3);

end
