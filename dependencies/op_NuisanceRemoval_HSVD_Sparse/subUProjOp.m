function y = subUProjOp( x, designPara)
% Adjoint of the forward operator for updating U with given V

% prepare parameters
NHr           = designPara.NHr;
NyHr          = NHr(1);
NxHr          = NHr(2);
RWat          = size(designPara.VtWat,1);
RLip          = size(designPara.VtLip,1);
RMet          = size(designPara.VtMet,1);

if isfield(designPara,'BigLambda')
    BigLambda = diag(designPara.BigLambda);
else
    BigLambda = 1;
end

% AcAx
Ax            = subUAOp( x, 'notransp_df', designPara);
AcAx          = subUAOp( Ax, 'transp_df', designPara);

% weighted l2 regularization
reg_wl2       = [];
x             = reshape(x,NyHr*NxHr,RWat+RLip+RMet);
if RWat       ~= 0
    UxWat     = x(:,1:RWat);
    temp      = designPara.lambdaWat*UxWat*BigLambda(1:RWat,1:RWat)*BigLambda(1:RWat,1:RWat)';
    reg_wl2   = [reg_wl2;temp(:)];
end
if RLip       ~= 0
    UxLip     = x(:,(RWat+1):(RWat+RLip));
    temp      = designPara.lambdaLip*UxLip*BigLambda((RWat+1):(RWat+RLip),(RWat+1):(RWat+RLip))...
                *BigLambda((RWat+1):(RWat+RLip),(RWat+1):(RWat+RLip))';
    reg_wl2   = [reg_wl2;temp(:)];
end
if RMet       ~= 0     
    UxMet     = x(:,(RWat+RLip+1):end);
    temp      = designPara.lambdaMet*UxMet*BigLambda((RWat+RLip+1):end,(RWat+RLip+1):end)...
                *BigLambda((RWat+RLip+1):end,(RWat+RLip+1):end)';
    reg_wl2   = [reg_wl2;temp(:)];
end

% cost function
y             = AcAx + reg_wl2(:);

end