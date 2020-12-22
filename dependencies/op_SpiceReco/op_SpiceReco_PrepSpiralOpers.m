function [SpiceOperators] = op_SpiceReco_PrepSpiralOpers(InputOperators,B0Info)
%
% LowRankTGV_RecurModel_PrepOper Prepare Operators for Antoines LowRankTGV Compressed Sensing Reco 
%
% This function was written by Bernhard Strasser, March 2019.
%
%
% The function can read in MRS(I) data in the Siemens raw file format ".DAT" and performs
% some easy Postprocessing steps like zerofilling, Hadamard decoding, Noise Decorrelation etc.
%
%
% [kSpace, Info] = read_csi_dat(file, DesiredSize,ReadInDataSets)
%
% Input: 
% -         Input                    ...     The Input. Must have fields 'Data', and 'Par'.
% -         CoilWeightMap             ...     The Coil-weighting-map. Must have field 'Data'
% -         Settings          ...            Structure to specify how the coil combination should be performed.
%
% Output:
% -         Input                      ...     The Output (Input is overwritten). Will have field 'RecoPar' additionally.
% -         AdditionalOut                        ...  Variable for storing additional output (if asked for by user).   
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: For now, the coil index in CoilWeightMap is the first, but in Input it's the last. Change that later!
% Also, the function cannot handle 3D / multislice data for now. And the channel index is fixed for now.


% This function expects the input to be of form
% ??? 


%% B0-Correction Operator

VSize = B0Info.RecoPar.vecSize;
DummyData.Par = B0Info.RecoPar;
DummyData.Data = ones(DummyData.Par.DataSize);
[Delme2,Dummy] = op_CorrSpectralB0(DummyData,B0Info.B0,B0Info.Mask,B0Info.Settings);
SpiceOperators.B0CorrMat_Spec = reshape(Dummy.B0CorrMat_Spec,[size(InputOperators.sft2_Oper,2) VSize]);
clear Delme2 DummyData Dummy


%% Mask

SpiceOperators.Mask = imresize(B0Info.Mask,B0Info.RecoPar.DataSize(1:2),'nearest');
SpiceOperators.Mask = SpiceOperators.Mask(:);


%% All Other Operators

% Just Make other operators to correct size
% Later, really calculate them from scratch
SpiceOperators.SamplingOperator = logical(squeeze(InputOperators.SamplingOperator(:,:,:,:,1:VSize,:)));
SpiceOperators.SamplingOperator = reshape(SpiceOperators.SamplingOperator,[numel(SpiceOperators.SamplingOperator(:,:,1)) VSize]);
SpiceOperators.TiltTrajMat = InputOperators.TiltTrajMat;
SpiceOperators.DCFPreG = InputOperators.DCFPreG(:) / (norm(InputOperators.DCFPreG(:))/sqrt(numel(InputOperators.DCFPreG)));

FoVMask = EllipticalFilter(ones(B0Info.RecoPar.DataSize(1:2)),[1 2],[1 1 1 B0Info.RecoPar.DataSize(1)/2-1],1); 
FoVMask = FoVMask(:);
SpiceOperators.sft2_Oper = InputOperators.sft2_Oper;
SpiceOperators.sft2_Oper(:,~logical(FoVMask(:))) = 0;

% SpiceOperators.SENSE = reshape(ones(size(InputOperators.SensMap(:,:,1,1,1))),[size(SpiceOperators.sft2_Oper,2) 1]); % Fake SensMap, if coilcomb was alrdy done

% Scale = B0Info.Mask ./ sqrt(sum(abs(InputOperators.SensMap(:,:,:,1,:)).^2,5));
Scale = 1 ./ sqrt(sum(abs(InputOperators.SensMap(:,:,:,1,:)).^2,5));
Scale(isinf(Scale) | isnan(Scale) | Scale == 0 ) = 1;
SpiceOperators.SENSE = reshape(InputOperators.SensMap(:,:,:,1,:),[size(SpiceOperators.sft2_Oper,2) 1 size(InputOperators.SensMap,5)]).*Scale(:).*SpiceOperators.Mask;
SpiceOperators.SENSE(isinf(SpiceOperators.SENSE) | isnan(SpiceOperators.SENSE)) = 0;
% SpiceOperators.SENSE(SpiceOperators.SENSE == 0) = 1;
% SpiceOperators.SENSE = SpiceOperators.SENSE / (norm(SpiceOperators.SENSE(:))/sqrt(numel(SpiceOperators.SENSE))) / 1;    % Normalize?!


% SpiceOperators.B0CorrMat_Spec = squeeze(InputOperators.B0CorrMat_Spec);


% / (norm(InputOperators.DCFPreG(:))/sqrt(numel(InputOperators.DCFPreG))): The rationale behind this is to rescale the operator having the same "energy" as an operator
% full of ones, which does nothing to the data when applied element-wise. Such an 1-operator has energy sqrt(numel(x))


%% old:
% SpiceOperators.DCFPreG = InputOperators.DCFPreG(:) / (norm(InputOperators.DCFPreG(:)));
% SpiceOperators.SENSE = reshape(InputOperators.SensMap(:,:,:,1,:),[size(SpiceOperators.sft2_Oper,2) 1 size(InputOperators.SensMap,5)]);


