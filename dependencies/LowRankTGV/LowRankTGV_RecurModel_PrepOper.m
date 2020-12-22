function [CSOperators] = LowRankTGV_RecurModel_PrepOper(InputOperators,B0Info)
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

DummyData.Par = B0Info.RecoPar;
DummyData.Data = ones(DummyData.Par.DataSize);
[Delme2,Dummy] = op_CorrSpectralB0(DummyData,B0Info.B0,B0Info.Mask,B0Info.Settings);
CSOperators.B0CorrMat_Spec = squeeze(Dummy.B0CorrMat_Spec);
clear Delme2 DummyData Dummy


%% All Other Operators

% Just Make other operators to correct size
% Later, really calculate them from scratch
VSize = B0Info.RecoPar.vecSize;
CSOperators.SamplingOperator = InputOperators.SamplingOperator(:,:,1,1,1,1);
InSize = size(CSOperators.SamplingOperator); 
CSOperators.TiltTrajMat = reshape(InputOperators.TiltTrajMat,[InSize(1:2) size(InputOperators.TiltTrajMat,2)]); %[InSize(1:2) size(InputOperators.TiltTrajMat,2)] actually...
CSOperators.DCFPreG = InputOperators.DCFPreG / norm(InputOperators.DCFPreG(:)/sqrt(numel(InputOperators.DCFPreG)));
CSOperators.sft2_Oper = InputOperators.sft2_Oper;

Scale = 1 ./ sqrt(sum(abs(InputOperators.SensMap(:,:,:,1,:)).^2,5));
Scale(isinf(Scale) | isnan(Scale) | Scale == 0 ) = 1;
CSOperators.SENSE = InputOperators.SensMap(:,:,:,1,:).*Scale.*B0Info.Mask;
CSOperators.SENSE(isinf(CSOperators.SENSE) | isnan(CSOperators.SENSE)) = 0;
CSOperators.SENSE = squeeze_single_dim(CSOperators.SENSE(:,:,:,1,:),3);


% CSOperators.SENSE = squeeze_single_dim(InputOperators.SensMap(:,:,:,1,:),3);
% CSOperators.B0CorrMat_Spec = squeeze(InputOperators.B0CorrMat_Spec);

CSOperators.Mask = imresize(B0Info.Mask,B0Info.RecoPar.DataSize(1:2),'nearest');
% CSOperators.Mask = InputOperators.Mask;




