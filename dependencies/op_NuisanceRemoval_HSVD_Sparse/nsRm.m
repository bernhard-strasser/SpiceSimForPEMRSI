function [xtns,xtnsr] = nsRm(xt, ImMask, B0Map, params, optWat, optLip, optMeta, optOther)
% water/fat removal
% Input:
% 	xt			: csi/epsi  data
% 	ImMask		: mask
% 	optWatFat   : parameters of water/fat removal algorithm for regions inside the mask
% 	optOther	: parameters of water/fat removal algorithm for regions out the mask
%   B0Map       :
%   NtEcho      :
%
% Output:
% 	xtns 	    : water/fat signal
% 	xtnsr 		: signal after water/fat removal
%
% Created by Chao Ma
% Last modified : 04-20-2014
% History:
% 04-20-2014    : Add field inhomogeneity
%

    % preparation
    Ny                                  = size(xt,1);
    Nx                                  = size(xt,2);
    Nt                                  = size(xt,3);
    xtns                                = zeros(Ny,Nx,Nt);
    xtnsr                               = zeros(Ny,Nx,Nt);
    if ~exist('B0Map','var') || isempty(B0Map)
        B0Map                           = zeros(Ny,Nx);
    end
    if sum(size(B0Map) ~= [Ny,Nx]) > 0
        B0Map                           = imresize(B0Map,[Ny,Nx]);
    end
    if ~exist('optOther','var')
        optOther                        = [];    
    end
    
    f                                   = linspace(-1/params.dt/2,1/params.dt/2,Nt);
    % nuisance signal removal point by point
    for indy = 1:Ny
        for indx = 1:Nx
            s_xt                        = xt(indy,indx,:);
            df                          = B0Map(indy,indx);
            if ImMask(indy,indx)        == 1
                s_xt                    = reshape(s_xt,1,[]);
                params.df               = df;
                s_xtWaterFat            = sigselhsvd(s_xt, params, optWat, optLip, optMeta);
                s_xtWaterFatRemoved     = s_xt - s_xtWaterFat;
            elseif ~isempty(optOther)
                s_xt                    = reshape(s_xt,1,[]);
                params.df               = df;
                s_xtWaterFat            = sigselhsvd(s_xt, params, optOther, optLip, optMeta);
                s_xtWaterFatRemoved     = s_xt - s_xtWaterFat;
            else
                s_xtWaterFat            = 0;
                s_xtWaterFatRemoved     = s_xt;

            end
            xtns(indy,indx,:)           = s_xtWaterFat;
            xtnsr(indy,indx,:)          = s_xtWaterFatRemoved;
        end
    end
end