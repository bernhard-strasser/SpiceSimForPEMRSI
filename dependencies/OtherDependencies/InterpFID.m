function [VD2,SD2,UD2] = InterpFID(VD1,SD1,UD1,tD1,tD2,svd_on)
%
% InterpFID This function interpolated FIDs.
%
% This function was written by Chao Ma, [month] [year].
%
%
% The function was taken from Chao Ma's NuisanceRemoval package for sparse data. It interpolates VD1 at time points tD1 to VD2 at time points tD2.
% VD1 can hold several FIDs, the time-dimension must be the second, the first dimension is the different FIDs. If svd_on is true, then also
% SD1 and UD1 must be given, which are the singular values matrix, and the left (spatial) singular vectors of an svd, i.e. [UD1,SD1,VD1] = svd(CsiData).
% In general, the functions is written for VD1s of such type, but I think all FIDs should be handled.
% 
%
%
% [A,B] = FunctionTemplateBstrasser(inputvar1,inputvar2)
%
% Input: 
% -         VD1                   ...    A matrix of size [N x NTime1], where each VD1(ii,:) hold and FID, and N is the number of FIDs (e.g. right (temp) sing val)
% -         SD1                   ...    Matrix of singular values from [UD1,SD1,VD1] = svd(CsiData). Only necessary if svd_on = true.
% -         UD1                   ...    Matrix of left (spatial) singular vectors from [UD1,SD1,VD1] = svd(CsiData). Only necessary if svd_on = true.
% -         tD1                   ...    Vector specifying the time points of VD1 in s. Must be of size NTime1.
% -         tD2                   ...    Vector specifying the time points to which the interpol should be performed. Size: 1xNTime2
% -         svd_on                ...    If true, another svd is performed after interpolation.
%
% Output:
% -         VD2                           ...     The FIDs interpolated to size [N x NTime2]
% -         SD2                           ...     If svd_on = true, another svd is performed on the interpolated data, and this is the resulting sing value matrix.
% -         UD2                           ...     If svd_on = true, another svd is performed on the interpolated data, and this is the resulting spat sing vector mat.
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: None

% Further remarks: 



%% 0. Preparations, Definitions

% % 0.1 Preparations
% if(~exist('inputvar1','var'))
%     inputvar1 = 0;
% end
% if(~exist('inputvar2','var'))
%     inputvar2 = 0;
% end

% 0.3 Definitions
    
%% Function: Interpolation

if ~isempty(VD1)
    if (length(tD1)~=length(tD2)) || norm(tD1-tD2) ~= 0 % if tD1 ~= tD2, do interpolation
        R                      = size(VD1,1);
        VD2                    = zeros(R,length(tD2));       
        for ind = 1:R
            VD2(ind,:)         = complex(interp1(tD1,real(squeeze(VD1(ind,:))),tD2,'spline'),...
                                         interp1(tD1,imag(squeeze(VD1(ind,:))),tD2,'spline'));
        end      
        if svd_on % do svd again
            temp               = UD1*VD2;
            [UD2,SD2,tVD2]     = svd(temp,'econ');
            VD2                = tVD2;
            UD2                = UD2*SD2;
        else
            SD2                = SD1;
            UD2                = UD1;
        end
    else % otherwise, no extra operation
        VD2                    = VD1;
        UD2                    = UD1;
        SD2                    = SD1;
    end
else
    VD2                        = [];
    UD2                        = [];
    SD2                        = [];
end       
      



%% 3. Postparations





