function [MRStruct,exp_filter_mat] = op_ExponentialFilter(MRStruct,Settings)
%
% ExponentialFilter Apply an Exponential filter to time-domain signals (e.g. FIDs)
%
% This function was written by Wolfgang Bogner 2013, revised by Bernhard Strasser, October 2013.
%
%
% The function computes an exponential filter in Hertz
%
%
% [MRStruct,exp_filter_funct] = ExponentialFilter(InArray,Settings)
%
% Input: 
% -         InArray                     ...    Input array to which the filter should be applied
% -         dwelltime                   ...    The dwelltime in [ns], i.e. the time between two consecutive time points.
% -         Settings.ApplyAlongDim               ...    Along this dimension the filter is applied. 
% -         Settings.ExpFilterStrength_Hz               ...    The filter strength in Hz
%
% Output:
% -         MRStruct                    ...     The filtered output array
% -         exp_filter_funct            ...     The filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy:

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% Define variables if not defined


if(~isfield(MRStruct,'Par'))
    MRStruct.Par = struct;
    if(isfield(MRStruct,'Data_file'))
        MRStruct.Par = read_ascconv(MRStruct.Data_file); 
    end
end
if(~isfield(MRStruct.Par,'DataSize'))
    MRStruct.Par.DataSize = size(MRStruct.Data);
end
if(~isfield(MRStruct,'RecoPar'))
    MRStruct.RecoPar = MRStruct.Par;
end

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'ExpFilterStrength_Hz'))
    Settings.ExpFilterStrength_Hz = 3;
end
if(~isfield(Settings,'ApplyAlongDim'))
    Settings.ApplyAlongDim = numel(size(MRStruct));
end


% Define vecSize
vecSize = size(MRStruct.Data,Settings.ApplyAlongDim);

% Define Dwelltime
dwelltime = MRStruct.RecoPar.Dwelltimes(1);



%% 1. Compute exponential Time-Domain Filter

dwelltime_in_s = dwelltime/1000000000;

t= 0:dwelltime_in_s:dwelltime_in_s*(vecSize-1);
exp_filter_funct = exp(-Settings.ExpFilterStrength_Hz*t);     %exp(-t/a) wobei "1/a" = "exp_filter" Linebroadening in Hz




%% 2. Replicate Filter to size of InArray

exp_filter_mat = myrepmat(exp_filter_funct,size(MRStruct.Data));





%% 3. Apply Hamming Filter


MRStruct.Data = exp_filter_mat.*MRStruct.Data;

if(isfield(MRStruct,'NoiseData'))
    MRStruct.NoiseData = exp_filter_mat.*MRStruct.NoiseData;
end



