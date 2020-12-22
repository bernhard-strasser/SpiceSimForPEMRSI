function [U_rrrc,V_tc,S,costFunVal] = tgv2_l2_2D_multiCoil_LowRank_CombinedConv_bstr(Data_kkkstc,A,Ah,U_rrrc,V_tc,S, alpha0, alpha1, maxits,minits,mrsiReconParams,Threshold,Settings)
% Primal dual TGV2 algorithm

if(~exist('Settings','var'))
    Settings = struct;
end
if(~isfield(Settings,'WriteFiles_flag'))
    Settings.WriteFiles_flag = false;
end
if(~isfield(Settings,'WriteFiles_Path'))
    Settings.WriteFiles_Path = './LowRankTGVRecon_DebugOutput'; 
end


[dxm,dym,dzm,dxp,dyp,dzp] = defineDiffOperators();
check_it = mrsiReconParams.LRTGVModelParams.check_it;%25; %Will show several converge values at evrey *this* steps
Plot_it=mrsiReconParams.LRTGVModelParams.Plot_it;% 50 %Will plot figures every *this* steps

Orthogonalize_it=mrsiReconParams.LRTGVModelParams.Orthogonalize_it;%=20 %Orthogonalize Sepctral component  every *this* step
SpecItFact=mrsiReconParams.LRTGVModelParams.SpecItFact;%2  % Will run *this* more spectral gradient descent than spatial convergence steps
reduction =mrsiReconParams.LRTGVModelParams.reduction;% 100 or 1000 if diverge. Starting with  alpha0/reduction, alpha1/reduction at the begining then going back to alpha0, alpha1 original values

min_SpectStep=mrsiReconParams.LRTGVModelParams.min_SpectStep;%=1/128 Minimum Spectral step size
max_SpectStep=mrsiReconParams.LRTGVModelParams.max_SpectStep;%=1/4 % decrease if diverges

min_taup=mrsiReconParams.LRTGVModelParams.min_taup;%=1/128 Minimum Spatial step size
max_taup=mrsiReconParams.LRTGVModelParams.max_taup;%=1/16 % decrease if diverges

CorrB0Map_count=1;

% Assign useful variables
Init_U_rrrc= U_rrrc;
[ M N O NbComp] = size(U_rrrc); % numSamplesOnSpoke, numSamplesOnSpoke, nCh
NbCoil = size(Data_kkkstc,6);
NbT= size(Data_kkkstc,5);
SizeVol = size(U_rrrc);
DimVol = ndims(U_rrrc)-1;
UIndcs  = repmat({':'}, [1, numel(SizeVol)]);
numSpatialPts=size(U_rrrc, 1)*size(U_rrrc,2);
% BMask_rrt=repmat(mrsiReconParams.CSOperators.BrainMask2D,[1 1 NbT]);
% kmask_ckkt=permute(repmat(mrsiReconParams.kmask,[1 1 NbCoil NbT]),[3,1,2,4]);
% kmask_ckkt = mrsiReconParams.CSOperators.SamplingOperator;

%Prepare operator for Forward Transform %Not Needed with A and Ah functions
%{ 
FreqMap=mrsiReconParams.WaterFreqMap ;
Fs=mrsiReconParams.mrProt.samplerate*NbT/mrsiReconParams.mrProt.VSize;
DelayT=mrsiReconParams.DelayT;
Time_rrt=permute(repmat(([0 :(NbT-1)]'/Fs+DelayT),[1 M N]),[2,3,1]);
Freqshift_rrt=exp(2*pi*1i*Time_rrt.*repmat(FreqMap,[1 1 NbT]));
RepData_crrt=0*Data_kkkstc;
CompData_kkkstc=0*Data_kkkstc;
kmask_ckkt=permute(repmat(mrsiReconParams.kmask,[1 1 NbCoil NbT]),[3,1,2,4]);
SENSE_crrt=repmat(mrsiReconParams.SENSE,[1 1 1 NbT]);
TempData=0*SENSE_crrt;
G=0*formTensorProduct(U_rrrc, V_tc*S,DimVol);
%}

v = 0*Data_kkkstc;%zeros([NbCoil SizeVol]); % coil-k-k-c
p = zeros([SizeVol,2]);
q = zeros([SizeVol,3]);
xi = zeros([SizeVol,2]);
xi_ = xi; % v in article

StepNormGrad = [];
StepDiff=Threshold;
StepDiffu = [];
StepDiffr = [];
StepDiffxi = [];
StepDiffq = [ ];
StepDiffp = [];
StepDiffww = [ ];
StepDiffdivp = [];

tau_p = min_taup;
tau_d = tau_p*2;
stepSize = min_SpectStep; 

PrevNormGradU=1E12;
NormGradWentDown=0;
NbDivergSteps=0;

RelDiffSq_V=1;
RelDiffSq_U=1;

t_old  = 2;
PrevNormGrad=0;
step_noConv=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign Iterative Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = U_rrrc;
uold = u; 
xiold = xi;
V_old = V_tc;
V_proj  = V_tc;

% First dataConsistencyCost calculation
dataConsistencyCost = norm(Data_kkkstc(:) - vectorizeArray(A(formTensorProduct(U_rrrc, V_tc*S,DimVol))));
costFunVal = dataConsistencyCost;

% Iteration Index
k=-1;

% Initial Reduced regularization values
alpha00 = alpha0/reduction;
alpha10 = alpha1/reduction;
% Middle-way Increased regularization values
alpha002 = alpha0*reduction;
alpha102 = alpha1*reduction;
alpha01 = alpha0;
alpha11 = alpha1;

while ( abs(StepDiff(end))>Threshold  & (k<maxits) ) | k<(minits); %| k<(minits+Orthogonalize_all_last_it);
    k=k+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPATIAL Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % update alpha's
   if k<=(minits/2)
  	 alpha0 = exp(k/(minits/2)*log(alpha01) + ((minits/2)-k)/(minits/2)*log(alpha00));
  	 alpha1 = exp(k/(minits/2)*log(alpha11) + ((minits/2)-k)/(minits/2)*log(alpha10));
  elseif k<=(minits)
   	 alpha0 = exp((k-minits/2)/(minits/2)*log(alpha01) + ((minits/2)-(k-minits/2))/(minits/2)*log(alpha002));
  	 alpha1 = exp((k-minits/2)/(minits/2)*log(alpha11) + ((minits/2)-(k-minits/2))/(minits/2)*log(alpha102));
  else   
	 alpha0 = alpha01;
  	 alpha1 = alpha11;
  end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SAVE VARIABLES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uold = u;
    xiold = xi;
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DUAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Compute gradients
    ux = dxp(U_rrrc);
    uy = dyp(U_rrrc);
    
    p(UIndcs{:},1) = p(UIndcs{:},1) - tau_d*(ux + xi_(UIndcs{:},1));
    p(UIndcs{:},2) = p(UIndcs{:},2) - tau_d*(uy + xi_(UIndcs{:},2));

    % projection    
    for comp=1:NbComp
        absp = sqrt(abs(p(UIndcs{1:(end-1)},comp,1)).^2 + abs(p(UIndcs{1:(end-1)},comp,2)).^2 );
        denom = max(1,absp/(alpha1*max(diag(S))/S(comp,comp)));
        p(UIndcs{1:(end-1)},comp,1) = p(UIndcs{1:(end-1)},comp,1)./denom;
        p(UIndcs{1:(end-1)},comp,2) = p(UIndcs{1:(end-1)},comp,2)./denom;
    end   
    
    % symmetrized gradient
    gradxi1 = dxm(xi_(UIndcs{:},1));
    gradxi2 = dym(xi_(UIndcs{:},2));
    gradxi3 = (dym(xi_(UIndcs{:},1)) + dxm(xi_(UIndcs{:},2)))/2;
    
    q(UIndcs{:},1) = q(UIndcs{:},1) - tau_d*gradxi1; % line
    q(UIndcs{:},2) = q(UIndcs{:},2) - tau_d*gradxi2;
    q(UIndcs{:},3) = q(UIndcs{:},3) - tau_d*gradxi3;
    
    
    for comp=1:NbComp
        absq = sqrt(abs(q(UIndcs{1:(end-1)},comp,1)).^2 + abs(q(UIndcs{1:(end-1)},comp,2)).^2 + 2*abs(q(UIndcs{1:(end-1)},comp,3)).^2);
        denom = max(1,absq/(alpha0*max(diag(S))/S(comp,comp)));
        q(UIndcs{1:(end-1)},comp,1) = q(UIndcs{1:(end-1)},comp,1)./denom;
        q(UIndcs{1:(end-1)},comp,2) = q(UIndcs{1:(end-1)},comp,2)./denom;
        q(UIndcs{1:(end-1)},comp,3) = q(UIndcs{1:(end-1)},comp,3)./denom;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRIMAL UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    G  = formTensorProduct(U_rrrc, V_tc*S,DimVol); %rrt	
    % Error v of the Forward transform to the raw data
    v =  (A(G) - Data_kkkstc ); % Data Fidelity gradient % Dimensions % coil-k-k-t
    
	%Adjoint Transform on the error v

    ww= formTensorProduct(Ah(v) ,(V_tc*inv(S))', DimVol ) ; % r - r - component
    
    % divergence
    divp = dxm(p(UIndcs{:},1)) + dym(p(UIndcs{:},2));
    
    u = u - tau_p*(ww + divp); %lines 8
      
    % divergence
    divq1 = dxp(q(UIndcs{:},1)) + dyp(q(UIndcs{:},3));
    divq2 = dxp(q(UIndcs{:},3)) + dyp(q(UIndcs{:},2));
    
    xi(UIndcs{:},1) = xi(UIndcs{:},1) - tau_p*(divq1 - p(UIndcs{:},1));%line 11
    xi(UIndcs{:},2) = xi(UIndcs{:},2) - tau_p*(divq2 - p(UIndcs{:},2));%line 11
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % AUXILIARY UPDATE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    U_rrrc =2*u - uold;
    xi_ = 2*xi - xiold;
    
	%Impose Component norm to be 1. Store the change in norm in S (originally the singular values)
    for c=1:NbComp
        NormUrrc=sqrt(sum(sum(abs(U_rrrc(:,:,c)).^2,1),2));
        S(c,c)=S(c,c)*NormUrrc;
        U_rrrc(:,:,c) = U_rrrc(:,:,c)/NormUrrc;
        u(:,:,c) = u(:,:,c)/NormUrrc;
    end

    RelDiffSq_U = norm(U_rrrc(:)-u(:))^2/norm(u(:))^2 + norm(xi_(:)-xi(:))^2/norm(xi(:))^2;
    norm_p=norm(p(:),1);
    norm_q=norm(q(:),1);
    RelDiffSq_old=RelDiffSq_U;
    
    NormGradU=norm(ww(:) + divp(:));
    
    %Adjust step size base on previous gradient value
    if k>5 & NormGradU<PrevNormGradU;
        NormGradWentDown=1;
        tau_p = tau_p*1.40;
        if tau_p>max_taup; tau_p=max_taup;end;
    elseif NormGradWentDown==1;
        tau_p=tau_p*0.55;
        if tau_p<min_taup;tau_p=min_taup;end;
    end
    tau_d=tau_p*2;
    PrevNormGradU=NormGradU;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPECTRAL Update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for SpecInIt=1:SpecItFact


         G  = formTensorProduct(U_rrrc, V_tc*S,DimVol); %rrt	
        % Error v of the Forward transform to the raw data
         v =  (A(G) - Data_kkkstc ); % Data Fidelity gradient % Dimensions % coil-k-k-t
        %Adjoint transform and scalar product over spatial component to obtain spectral components gradient 
        gradStep = ((reshape(U_rrrc, [], NbComp)*inv(S))' * reshape(Ah(v) , [numSpatialPts, NbT]))' ;  % t - component

        % Orthogonalize the gradient step (make each component direction orthogonal)        
        [ gradStep] = OrthogonalizeTimeComponents(gradStep);

        V_proj   = V_tc -stepSize * gradStep;
        
        %Spectral components norm must always be 1
        for c=1:size(V_tc,2)
              V_proj(:,c) = V_proj(:,c) ./ norm(V_proj(:,c));
        end
       
        if k>0 
            RelDiffSq_V=norm(V_proj(:)-V_old(:))/norm(V_proj(:));
        end
        %Adjust Gradient descent step size base on previous gradient value
        NormGrad=norm(gradStep(:));
        if (NormGrad<PrevNormGrad || k==0)
            step_noConv=0;
            stepSize = 1.2 * stepSize;
            if stepSize>max_SpectStep;stepSize=max_SpectStep;end
        else
            step_noConv=step_noConv+1;
            stepSize = 0.75 * stepSize;
            if stepSize<min_SpectStep;stepSize=min_SpectStep;end
        end
        
        PrevNormGrad=NormGrad;
        StepNormGrad = [StepNormGrad NormGrad];
        
        t_new    = (1 + sqrt(1 + 4*t_old.^2)) / 2;
        fact=((t_old - 1) / t_new);
        % V_update = V + ((t_old - 1) / t_new) * (V - V_old);
        V_tc = (1-fact)*V_old +(fact)* V_proj;
       
        for c=1:size(V_tc,2)
              V_tc(:,c) = V_tc(:,c) ./ norm(V_tc(:,c));
        end
        V_old    = V_proj;
        t_old    = t_new;
          
    end
    if  mod(k+1,Orthogonalize_it)==0
        V_tc = OrthogonalizeTimeComponents(V_tc);
        V_old= OrthogonalizeTimeComponents(V_old);
    end

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % keep other norms for monitoring
    
    norm_xi=norm(xi_(:));
    norm_u=norm(U_rrrc(:));%norm(ww(:) + divp(:));
    norm_r=norm(NormGradU(:));
    norm_divp=norm(divp(:));
    norm_ww=norm(ww(:));
    
  
    StepDiffu = [ StepDiffu norm_u];
    StepDiffr = [ StepDiffr norm_r];
    StepDiffxi = [ StepDiffxi norm_xi];
    StepDiffq = [ StepDiffq norm_q];
    StepDiffp = [ StepDiffp norm_p];
    StepDiffww = [ StepDiffww norm_ww ];
    StepDiffdivp = [ StepDiffdivp norm_divp ];
    
	% Display Convergence values 	
    if mod(k+1,check_it) == 0
        
        dataConsistencyCost = norm(Data_kkkstc(:) - vectorizeArray(A(formTensorProduct(U_rrrc, V_tc*S,DimVol))));
        StepDiff = [ StepDiff, (dataConsistencyCost-costFunVal(end))/(dataConsistencyCost*check_it)];
        costFunVal = [costFunVal, dataConsistencyCost];

        fprintf('\nTGV2-L2-2D: it = %g, costFunDiff = %g,tau_p = %g, NormGradU = %g ', k+1,StepDiff(end) ,tau_p, NormGradU);
        if(SpecItFact > 0)  % bstr: No need to print anything if we don't update V
            fprintf('\nSpectral Iteration: %d , spectral step_size = %d, RelDiffSq_V = %d, step_noConv = %d', k+1, stepSize,RelDiffSq_V,step_noConv);
        end
        CovV=V_tc'*V_tc;
        U_rc =reshape(U_rrrc,[],size(U_rrrc,3));
        CovU=(U_rc'*U_rc);
        fprintf('\nSpectral component independance: %d , Spatial component independance: %d\n', trace(abs(CovV))/sum(abs(CovV(:))),  trace(abs(CovU))/sum(abs(CovU(:))) );
        
    end
    % Print intermediate Components in eps files	
    if  (Settings.WriteFiles_flag && exist(Settings.WriteFiles_Path,'dir') && mod(k+1,Plot_it)==0)
        VisualizeTGV( Init_U_rrrc ,U_rrrc,[Settings.WriteFiles_Path '/LowRank-TGV_Spatial_Components_muTV', num2str(mrsiReconParams.mu_tv),'_step', num2str(k+1)]);
        VisualizeSpectral( V_tc,S, [ Settings.WriteFiles_Path '/LowRank-TGV_Spectral_Components_step', num2str(k+1)])
        close all;
    end
end

%Reorder Component following singular values:
[SDiag, DescOrder]=sort(diag(S),'descend');
U_rrrc=U_rrrc(:,:,:,DescOrder).* mrsiReconParams.CSOperators.Mask(:,:);
V_tc=V_tc(:,DescOrder);
S=diag(SDiag);
fprintf([ '\nTGV Recon done in ', num2str(k),' steps. Relative Step Diff U =',  num2str(RelDiffSq_U), ' steps. Relative Step Diff V =',  num2str(RelDiffSq_V),'\n']);

