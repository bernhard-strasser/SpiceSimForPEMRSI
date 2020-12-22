function params                = nsRmSparseInitParams(params)
    if ~isfield(params,'D1shift')
        params.D1shift         = [0,0];
    end
    if ~isfield(params,'D2shift')
        params.D2shift         = [0,0];
    end

    % mask for D1
    tImMask                    = logical(imresize(logical(params.ImMask),...
                                      round([params.ND1(1),params.ND1(2)]*params.NosD1)));
    twaterMask                 = logical(imresize(logical(params.waterMask),...
                                      round([params.ND1(1),params.ND1(2)]*params.NosD1)));
    tImMask                    = circshift(tImMask,params.D1shift);
    twaterMask                 = circshift(twaterMask,params.D1shift);
    tImMask                    = bwmorph(tImMask,'dilate');
    twaterMask                 = bwmorph(twaterMask,'erode');
    tlipMask                   = tImMask - twaterMask;
%     params.waterMaskD1         = tImMask;         % bstr changed: Somehow I get so much lipids in the lipid layer, which is included in tImMask. And there are signals
    params.waterMaskD1         = twaterMask;        % at ~5 ppm, which are treated as water often times. Try therefore smaller mask.
    params.metaMaskD1          = twaterMask;
    params.lipMaskD1           = tlipMask;

    % mask for D2
    tImMask                    = logical(imresize(logical(params.ImMask),...
                                      round([params.ND2(1),params.ND2(2)])));
    twaterMask                 = logical(imresize(logical(params.waterMask),...
                                      round([params.ND2(1),params.ND2(2)])));
    tImMask                    = circshift(tImMask,params.D2shift);
    twaterMask                 = circshift(twaterMask,params.D2shift);
    tImMask                    = bwmorph(tImMask,'dilate');
    twaterMask                 = bwmorph(twaterMask,'erode');
    tlipMask                   = tImMask - twaterMask;
    params.metaMaskD2          = twaterMask;
    params.lipMaskD2           = tlipMask;
    params.waterMaskD2         = tImMask;

    % field map
    params                     = setDefaultField(params,'B0Corr',1);
    if params.B0Corr == 1
        B0MapRsCsi             = imresize(params.B0Map,round([params.ND1(1),params.ND1(2)]*params.NosD1));
        B0MapRsCsi             = circshift(B0MapRsCsi,params.D1shift);
        params.B0MapD1Wat      = B0MapRsCsi;
        params.B0MapD1Met      = B0MapRsCsi;
        params.B0MapD1Lip      = B0MapRsCsi;

        B0MapRsEpsi            = imresize(params.B0Map,round([params.ND2(1),params.ND2(2)]));
        B0MapRsEpsi            = circshift(B0MapRsEpsi,params.D2shift);
        params.B0MapD2Wat      = B0MapRsEpsi;
        params.B0MapD2Met      = B0MapRsEpsi;
        params.B0MapD2Lip      = B0MapRsEpsi;
    else
        params.B0MapD1Wat      = [];
        params.B0MapD1Met      = [];
        params.B0MapD1Lip      = [];
        params.B0MapD2Wat      = [];
        params.B0MapD2Met      = [];
        params.B0MapD2Lip      = [];
    end
    
    % selParams
    selParams                  = [];
    selParams.dt               = params.dt;
    selParams.NtEcho           = params.csiInfo.NtEcho;
    params.selParams           = selParams;

    % spectral parameters
    [optWat,optLip,optMeta]    = setDefaultSpectInfo;
    params                     = setDefaultField(params,'optWat',optWat);
    params                     = setDefaultField(params,'optLip',optLip);
    params                     = setDefaultField(params,'optMeta',optMeta);

    % regularization
    params                     = setDefaultField(params,'specMask',...
                                 ones(params.ND2(1)*params.ND2(2),length(params.tD1)-1));
    params                     = setDefaultField(params,'lambdaWat',1e-4);
    params                     = setDefaultField(params,'lambdaMet',1e-4);
    params                     = setDefaultField(params,'lambdaLip',1e-4);
    params                     = setDefaultField(params,'lambdaDfWat',sqrt(5e-3));
    params                     = setDefaultField(params,'lambdaDfMet',0);
    params                     = setDefaultField(params,'lambdaDfLip',sqrt(5e-3));

    % algorithm related
    params                     = setDefaultField(params,'noMeta',0);
    params                     = setDefaultField(params,'NosD1',2);
    params                     = setDefaultField(params,'W', 0.01);
    params                     = setDefaultField(params,'NUWeights',1);
    params                     = setDefaultField(params,'InterpSvdOn',0);
    params                     = setDefaultField(params,'zeroInitial',1);
    params                     = setDefaultField(params,'tolCG',1e-3);
    params                     = setDefaultField(params,'MaxIterCG',200);
    params                     = setDefaultField(params,'Nc',1);
    params                     = setDefaultField(params,'verbose',1);
    params                     = setDefaultField(params,'debug_on',1);
end


