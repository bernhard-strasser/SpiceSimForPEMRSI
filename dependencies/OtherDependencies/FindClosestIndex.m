function MatchIndex = FindClosestIndex(InVec,MatchValue,Epsilon)
%
% FindClosestIndex Find Index of Closest Match of Vector to Provided Values Up to Epsilon
%
% This function was written by Bernhard Strasser, April 2018.
%
%
% The function gives the index of . 
%
%
% MatchIndex = FindClosestIndex(InVec,MatchValue)
%
% Input: 
% -         InVec                     ...     Vector for which the indices should be found.
% -         MatchValue                ...     Value(s) for which the closest match in InVec should be found.
% -         Epsilon                   ...     Both indices i,j are treated equally as matches, if abs(InVec(i) - InVec(j)) <= Epsilon and i is a match 
%
% Output:
% -         MatchIndex                ...     Index or Indices of the matched values. For each entry in "MatchValue", one cell with all the
%                                             found matches is provided. If several matches are found for one "MatchValue", all their positions
%                                             are provided
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None






%% 0. Preparations & Housekeeping

if(~exist('Epsilon','var'))
    Epsilon = eps;
end


%% 1. 

MatchIndex = cell([1 numel(MatchValue)]);
for CurValInd = 1:numel(MatchValue)
    MatchIndex{CurValInd} = find( abs(InVec - MatchValue(CurValInd)) <=  min(abs(InVec - MatchValue(CurValInd))) + Epsilon );
end


