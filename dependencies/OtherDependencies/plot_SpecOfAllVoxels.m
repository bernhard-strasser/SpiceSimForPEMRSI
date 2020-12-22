function plot_SpecOfAllVoxels(MRStruct, Settings,Mask)
%
% op_PermuteMRData Permute MR Data
%
% This function was written by Bernhard Strasser, Oct 2019.
%
%
% The function just runs the "permute" function on MR data as specified by Settings.PermuteVec.
% It also takes care of things like updating the .Par.DataSize, and permuting the NoiseData, if available.
%
%
% [MRStruct] = op_PermuteMRData(MRStruct,Settings)
%
% Input: 
% -         MRStruct               ...      The structure containing the MR data. Can have different fields. For this function, it
%                                           it has to have field
%                                           * .Data
% -         Settings               ...      Structure with settings how data should be processed. Should have field:
%                                           * .PermuteVec: Data will be permuted to according to this vector.
%
% Output:
% -         MRStruct               ...      The structure containing the MR data. Can have different fields.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: memused_linux,Analyze_csi_mdh_1_3, read_ascconv_1_2, hadamard_encoding_7.m

% Further remarks: 

%% Preparations

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'UseThisInStructMask'))
    Settings.UseThisInStructMask = 'BrainMask';
end

if(~isstruct(MRStruct))
    bak = MRStruct; clear MRStruct; MRStruct.Data = bak; clear bak; 
end
if(isfield(Settings,'UseThisInStructMask') && ~exist('Mask','var') && isfield(MRStruct,(Settings.UseThisInStructMask)))
    Mask = MRStruct.(Settings.UseThisInStructMask);
end
if(~exist('Mask','var'))
    Mask = ones(size_MultiDims(MRStruct.Data,1:3));
end



%% 

MRStruct.Data = MRStruct.Data(:,:,:,:,1);
MRStruct.Data = MRStruct.Data .* Mask;
MRStruct.Data = abs(fftshift(fft(MRStruct.Data,[],4),4));
MRStruct.Data = permute(MRStruct.Data,[4 1 2 3]);


if(isfield(MRStruct,'RecoPar'))
    chemy = compute_chemshift_vector(MRStruct.RecoPar);
elseif(isfield(MRStruct,'Par'))
    chemy = compute_chemshift_vector(MRStruct.Par);
else
    chemy = 1:size(MRStruct.Data,1); 
end
    
    
% convert from PPM to freq pts
if(isfield(Settings,'PlotPPMRange'))
    Settings.f_first = FindClosestIndex(chemy,max(Settings.PlotPPMRange)); Settings.f_first = Settings.f_first{1};
    Settings.f_end = FindClosestIndex(chemy,min(Settings.PlotPPMRange)); Settings.f_end = Settings.f_end{1};
    MRStruct.Data = MRStruct.Data(Settings.f_first:Settings.f_end,:,:,:,:);
    chemy = chemy(Settings.f_first:Settings.f_end);
end


figure; plot(chemy,sum(MRStruct.Data(:,:),2))
figure; plot(chemy,MRStruct.Data(:,:))



