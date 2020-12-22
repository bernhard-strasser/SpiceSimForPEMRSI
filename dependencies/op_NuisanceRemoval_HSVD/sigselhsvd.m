function xtSel            = sigselhsvd(xt, params, optWat, optLip, optMeta)
% signal selection using harmonics retrieval
% Chao Ma
% Last modified: 04-20-2014
% History:
% 04-20-2014: 1) Add B0 map information


    % data preparation
    temp                  = size(xt);
    Nr                    = temp(1);
    Nt                    = temp(2);    
    dt                    = params.dt;

    % set default values
    [params, optWat, optLip, optMeta] ...
                          = initParams(params, optWat, optLip, optMeta);

    % select signals based on HSVD
    z                     = cell(Nr,1);
    a                     = cell(Nr,1);
    df                    = cell(Nr,1);
    for indr = 1:Nr
        % Apply HSVD to RHS signal 
        xt_R              = reshape(xt(indr,params.NtEcho+1:end),[],1);
        % estimate noise
        temp              = mfft1d(xt_R,[],1);
        noise_std         = 2*std(temp(end-50:end));
        NR                = length(xt_R);
        L                 = round(NR/2);
%         [tz2,tfk2,tZ,ta]  = hsvd(xt_R, L, NR, params.n, dt);
        [tz2,tfk2]  = hsvd(xt_R, L, NR, params.n, dt);
%         while sum(isnan(ta))>0
%             params.n      = round(params.n/2);
%             [tz2,tfk2,tZ,ta]= hsvd(xt_R, L, NR, params.n, dt);
%         end

        % enforce nonnegative T=2
        if params.realT2 == 1
            tind          = find(abs(tz2)<1);
        else
            tind          = 1:length(tz2);
        end
        tz_R              = tz2(tind);
        tfk_R             = tfk2(tind);
        tT2_R             = -params.dt./log(abs(tz_R(:)))*1e3;
%         tZ_R              = tZ(:,tind);
%         ta_R              = ta(tind);
        tZ_R = 0; ta_R = 0;

        % find water, lip and meta signal
        indWat_R          = selectSignal(tfk_R, tz_R, params, optWat, tZ_R, ta_R, noise_std, optMeta);
        indLip_R          = selectSignal(tfk_R, tz_R, params, optLip, tZ_R, ta_R, noise_std, optMeta);
        indMeta_R         = selectSignal(tfk_R, tz_R, params, optMeta, tZ_R, ta_R, noise_std, optMeta);

        % forward matrix for RHS
        AR_Lip            = getA(NR,tz_R(indLip_R));
        AR_Wat            = getA(NR,tz_R(indWat_R));
        AR_Meta           = getA(NR,tz_R(indMeta_R));
        AR                = [AR_Lip,AR_Wat,AR_Meta];
        xtSel_R           = prepOut(params,AR,xt_R,indWat_R,indLip_R,indMeta_R);

        if params.NtEcho == 0 % if FID 
            xtSel(indr,:) = xtSel_R;           
        else % if echo signal
            % Apply HSVD to LHS signal 
            xt_L          = reshape(xt(indr,1:params.NtEcho),[],1);
            txt_L         = conj(xt_L(end:-1:1));
            NL            = length(xt_L);
            [tz2,tfk2]    = hsvd(txt_L, round(NL/2), NL, params.n, dt);

            % enforce nonnegative T=2
            if params.realT2 == 1
                tind      = find(abs(tz2)<1);
            else
                tind      = 1:length(tz2);
            end
            tz_L          = tz2(tind);
            tfk_L         = tfk2(tind);
            tT2_L         = -params.dt./log(abs(tz_L(:)))*1e3;

            % find water, lip and meta signal
            if ~isempty(optLip)
                optLip.maxT2 ...
                          = optLip.maxT2 + 10;
            end
            if ~isempty(optMeta)
                optMeta.minT2 ...
                          = optMeta.minT2 + 10;
            end
            indWat_L1     = selectSignal(tfk_L, tz_L, params, optWat);
            indLip_L1     = selectSignal(tfk_L, tz_L, params, optLip);
            indMeta_L1    = selectSignal(tfk_L, tz_L, params, optMeta);
            deltaf        = 5;
            indWat_L2     = freqSearch(tfk_L, tfk_R(indWat_R), deltaf);
            indLip_L2     = freqSearch(tfk_L, tfk_R(indLip_R), deltaf);
            indMeta_L2    = freqSearch(tfk_L, tfk_R(indMeta_R), deltaf);
            indWat_L      = union(indWat_L1,indWat_L2);
            indLip_L      = union(indLip_L1,indLip_L2);
            indMeta_L     = union(indMeta_L1,indMeta_L2);

            % forward matrix for LHS
            AL_Lip        = getA(NL,tz_L(indLip_L));
            AL_Wat        = getA(NL,tz_L(indWat_L));
            AL_Meta       = getA(NL,tz_L(indMeta_L));
            AL            = [AL_Lip,AL_Wat,AL_Meta];
            txtSel_L      = prepOut(params,AL,txt_L,indWat_L,indLip_L,indMeta_L);
            xtSel_L       = conj(txtSel_L(end:-1:1));

            xtSel(indr,:) = [xtSel_L(:);xtSel_R(:)];
        end
    end
end

% set default values
function [params, optWat, optLip, optMeta] = initParams(params, optWat, optLip, optMeta) 
    % default parameters of params
    params                = setDefaultField(params,'realT2',0);
    params                = setDefaultField(params,'InherentRef',0);
    params                = setDefaultField(params,'LSQ2',0);
    params                = setDefaultField(params,'NtEcho',0);
    params                = setDefaultField(params,'n',15);
    params                = setDefaultField(params,'threshold',0);
    params                = setDefaultField(params,'T2pMin',2);
    params                = setDefaultField(params,'T2pMax',100);

    % default optWat
    if ~isempty(optWat)
        optWat            = setDefaultField(optWat,'maxT2',[]);
        optWat            = setDefaultField(optWat,'minT2',[]);
    end

    % default optLip
    if ~isempty(optLip)
        optLip            = setDefaultField(optLip,'maxT2',[]);
        optLip            = setDefaultField(optLip,'minT2',[]);
    end

    % default optMeta
    if ~isempty(optMeta)
        optMeta           = setDefaultField(optMeta,'maxT2',[]);
        optMeta           = setDefaultField(optMeta,'minT2',[]);
    end
end

% select signals based on frequency and T2
function indSel           = selectSignal(tfk, tz, params, opt, tZ, ta, noise_std, opt2)
    if ~isempty(opt)
        % select the signal of interest
        % set frequency range
        if params.InherentRef == 1
            temp          = find((min(opt.fSel(:))<tfk) & (tfk<=max(opt.fSel(:))));
            [~,ind]       = max(abs(ta(temp)));
            temp2         = tfk(temp);
            if ~isempty(ind)
                df        = temp2(ind);
                tfSel     = df + opt.fSel2;
            else
                df        = 0;
                tfSel     = zeros(length(opt.maxT2)+length(opt.minT2)+1,1);
            end
        else
            if isfield(params,'df') && ~isempty(params.df)
                df        = params.df;
                tfSel     = df + opt.fSel2;
            else
                df        = 0;
                tfSel     = opt.fSel2;
            end
        end

        % find the signal of interest based on frequency and T2 values
        tT2               = -params.dt./log(abs(tz(:)))*1e3;
        indSel            = [];
        for ind           = 1:(length(opt.maxT2)+length(opt.minT2))
            if isempty(opt.minT2)
                temp      = find((min(tfSel(ind,:))<tfk) & (tfk<=max(tfSel(ind,:))) ...                 % bstr: added "real(tT2(:))>0", bc
                             & ((real(tT2(:))<= opt.maxT2(ind)) & real(tT2(:))>0 | imag(tT2(:))~=0));   % negative T2's make no much sense...?
            else
                temp      = find((min(tfSel(ind,:))<tfk) & (tfk<=max(tfSel(ind,:))) ...
                             & ((real(tT2(:))>= opt.minT2(ind)) | imag(tT2(:))~=0));
            end
            indSel        = [indSel;temp(:)];
        end
    else
        indSel            = [];
    end
    if params.threshold > 0
        tindSel           = [];
        tfSel2             = df + opt2.fSel2;
        for ind = 1:length(indSel)
            if ( (tfk(ind)>min(tfSel2(:))) & (tfk(ind)<max(tfSel2(:))) )
                temp        = tZ(:,indSel(ind))*ta(indSel(ind));
                tempf       = mfft1d(temp,[],1);
                tempf_max   = max(abs(tempf(:)));
                if tempf_max > params.threshold*noise_std
                    tindSel = [tindSel;indSel(ind)];
                end
            else
                tindSel   = [tindSel;indSel(ind)];
            end
        end
        indSel            = tindSel;
    end
end

% frequency search
function indSel           = freqSearch(tfk, f, deltaf)
    if ~exist('deltaf','var')
        deltaf            = 5;
    end
    n                     = length(f);
    indSel                = [];
    for ind = 1:n
        tind              = find(abs(tfk-f(ind))<deltaf);
        if ~isempty(tind)
            temp          = abs(tfk(tind)-f(ind));
            [~,temp2]     = min(temp(:));
            indSel        = [indSel;tind(temp2(1))];
        end
    end
end
% build forward matrix for right-hand-side
function A                = getA(N,tz)
    if ~isempty(tz)
        A                 = zeros(N,length(tz));
        tz = double(tz);    % bstr: So that "inf"-values are less likely in A, if abs(tz(p)) > 1.
        for m = 1:N
            for p = 1:length(tz)
                A(m,p)    = tz(p)^(m-1);
            end
        end
    else
        A                 = [];
    end
end

% prepare output
function bSel            = prepOut(params,A,b,indWat,indLip,indMeta)
    % output
        switch params.signalType
            case 'lipid'
                indSel    = 1:length(indLip);
            case 'water'
                indSel    = length(indLip)+1:length(indLip)+length(indWat);
            case 'nuisance'
                indSel    = 1:length(indLip)+length(indWat);
            case 'meta'
                indSel    = length(indLip)+length(indWat)+1:length(indLip)+length(indWat)+length(indMeta);
            otherwise
                indSel    = indWat;
        end
        if ~isempty(indSel)
            if params.LSQ2
                ta        = pinv(A(:,indSel))*b;     
            else
                ta        = pinv(A)*b;
                ta        = ta(indSel);
            end
            bSel          = A(:,indSel)*ta;      
        else
            bSel          = zeros(length(b),1);
        end
end
