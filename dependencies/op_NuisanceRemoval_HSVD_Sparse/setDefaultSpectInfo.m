function [optWat,optLip,optMeta] = setDefaultSpectInfo
    % water signal
    optWat                       = [];
    optWat.fSel                  = [-80, 80];
    optWat.fSel2                 = [-100, -80, 200];
    optWat.maxT2                 = [10,1e6];

    % lipid signal
    optLip                       = [];
    optLip.fSel2                 = [-500, -450, -350, -300,-200];
    optLip.maxT2                 = [ 100,   20,    5,   10];

    % metabolite signal
    optMeta                      = [];
    optMeta.fSel2                = [-350, -80];
    optMeta.minT2                = 25;
end