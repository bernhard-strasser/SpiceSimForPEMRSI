function [xtnsr,xtns] = nsRm_bstrchanged(xt, ImMask, B0Map, params, optWat, optLip, optMeta, optOther, poolObj)
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
    xtns                                = zeros(Ny*Nx,Nt);
    xtnsr                               = zeros(Ny*Nx,Nt);
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

    [indy, indx]                        = ind2sub2d(Ny,1:Ny*Nx);
    if nargin < 9
        poolObj                         = [];
    end
    % nuisance signal removal point by point
    if isempty(poolObj)
        for ind = 1:Ny*Nx
            s_xt                        = reshape(xt(indy(ind),indx(ind),:),1,[]);
            if(all(s_xt == 0))
                continue
            end
            df                          = B0Map(indy(ind),indx(ind));
            if ImMask(indy(ind),indx(ind))        == 1
%                 [x,y] = ind2sub([Ny Nx],ind); fprintf('\nVoxel: %d, [%d, %d]',ind,x,y)    % DEBUG ONLY
                sigselhsvd_params       = params;
                sigselhsvd_params.n = sigselhsvd_params.n(indy(ind),indx(ind));
                sigselhsvd_params.df    = df;
                s_xtWaterFat            = sigselhsvd(s_xt, sigselhsvd_params, optWat, optLip, optMeta);
                s_xtWaterFatRemoved     = s_xt - s_xtWaterFat;
            elseif ~isempty(optOther)
                sigselhsvd_params       = params;
                sigselhsvd_params.df    = df;
                s_xtWaterFat            = sigselhsvd(s_xt, sigselhsvd_params, optOther, optLip, optMeta);
                s_xtWaterFatRemoved     = s_xt - s_xtWaterFat;
            else
                s_xtWaterFat            = 0;
                s_xtWaterFatRemoved     = s_xt;
            end
            xtns(ind,:)                 = s_xtWaterFat;
            xtnsr(ind,:)                = s_xtWaterFatRemoved;
        end
    else
        parfor ind = 1:Ny*Nx
            s_xt                        = reshape(xt(indy(ind),indx(ind),:),1,[]);
            df                          = B0Map(indy(ind),indx(ind));
            if ImMask(indy(ind),indx(ind))        == 1
                sigselhsvd_params       = params;
                sigselhsvd_params.df    = df;
                s_xtWaterFat            = sigselhsvd(s_xt, sigselhsvd_params, optWat, optLip, optMeta);
                s_xtWaterFatRemoved     = s_xt - s_xtWaterFat;
            elseif ~isempty(optOther)
                sigselhsvd_params       = params;
                sigselhsvd_params.df    = df;
                s_xtWaterFat            = sigselhsvd(s_xt, sigselhsvd_params, optOther, optLip, optMeta);
                s_xtWaterFatRemoved     = s_xt - s_xtWaterFat;
            else
                s_xtWaterFat            = zeros(1,Nt);
                s_xtWaterFatRemoved     = s_xt;
            end
            xtns(ind,:)                 = s_xtWaterFat;
            xtnsr(ind,:)                = s_xtWaterFatRemoved;
        end
    end

    xtns                                = reshape(xtns,[Ny,Nx,Nt]);
    xtnsr                               = reshape(xtnsr,[Ny,Nx,Nt]);
end