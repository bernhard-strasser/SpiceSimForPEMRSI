function [ dWL, dWrLr, output, params ] = nsRmSparse_bstrchanged( d, params )
% Removal of nuisance signals from sparsely sampled data
% Created by Chao Ma
% Last modified on 02-20-2013
%

    %% some default paramters
    if ~isfield(params,'tD2')
        params.tD2                 = params.tD1;
    end
    if(~isfield(params,'ShowSpectra_flag'))
        params.ShowSpectra_flag = false; 
    end
    params.NtTr                    = length(params.tD2);
    NyD1                           = params.ND1(1);
    NxD1                           = params.ND1(2);
    NtD1                           = params.ND1(3);
    NosD1                          = params.NosD1;
    NyD2                           = params.ND2(1);
    NxD2                           = params.ND2(2);
    NyD2In                         = params.ND2In(1);
    NxD2In                         = params.ND2In(2);    
    NtD2                           = length(params.tD2);

    % generate a noise matrix to estimate sigular values from noise
    tnoise                         = params.xtD1WrLr(2,2,end-100:end-30);
    noiseStd                       = std(tnoise(:));
    noiseMat                       = noiseStd/sqrt(2)*(randn(size(params.xtD1WrLr))+1i*randn(size(params.xtD1WrLr)));
    noiseMat                       = reshape(noiseMat,[round(NyD1*NxD1),NtD1]);
    noiseS                         = svd(noiseMat,'econ');
    
    % output
    output                         = [];
    %% if not given, extract temporal basis functions from D1 data
    % water basis functions
    if ~isfield(params,'VWD1')
        if ~isempty(params.xtD1W) && (params.RW~=0)
            % extract water signal
            disp('Estimating temporal basis functions of water signals ...');
            tic;
            selParams              = params.selParams;
            selParams.signalType   = 'water';
            selParams.NtEcho       = params.csiInfo.NtEcho;
            [UWD1,SWD1,VWD1,params.RW] ...
                                   = estVt(params.xtD1W,params.B0MapD1Wat,params.tD1,params.waterMaskD1,...
                                     NosD1, params.ND1, params.optWat, params.optLip,params.optMeta, ...
                                     params.RW, selParams, params.verbose,noiseS);
            temp                   = toc;
            if params.debug_on >0
                disp(['Estimation finished in ',num2str(temp),' s']);
                if params.debug_on == 2
                    f              = linspace(-1/params.dt/2,1/params.dt/2,NtD1);
                    show_basis(fftshift(fft(VWD1,[],2),2),f,900);
                    set(gcf,'Name','Basis functions for water signal');
                    pause(1);
                end
            end
        else
            VWD1                   = [];
            SWD1                   = [];
            UWD1                   = [];
        end
    else
        VWD1                       = params.VWD1;
        if isfield(params,'SWD1')
            SWD1                   = params.SWD1;
        else
            SWD1                   = ones(size(VWD1,1),1);
        end
        if isfield(params,'UWD1')
            UWD1                   = params.UWD1;
        else
            UWD1                   = ones(NyD1*NxD1,size(VWD1,1));   
        end      
    end
    output.VWD1                    = VWD1;
    output.SWD1                    = SWD1;
    output.UWD1                    = UWD1;

    % lipid basis functions
    if ~isfield(params,'VLD1')
        if ~isempty(params.xtD1WrL) && (params.RL~=0)
            % extract lip signal
            disp('Estimating temporal basis functions of fat signals ...');
            tic;
            selParams              = params.selParams;
            selParams.signalType   = 'lipid';
            selParams.NtEcho       = params.csiInfo.NtEcho;
            

            
            [ULD1,SLD1,VLD1,params.RL] ...
                                   = estVt(params.xtD1WrL,params.B0MapD1Lip,params.tD1,params.lipMaskD1,...
                                     NosD1, params.ND1, params.optWat, params.optLip,params.optMeta, ...
                                     params.RL, selParams, params.verbose,noiseS);
            temp                   = toc;
            if params.debug_on >0
                disp(['Estimation finished in ',num2str(temp),' s']);
                if params.debug_on == 2
                    f              = linspace(-1/params.dt/2,1/params.dt/2,NtD1);
                    show_basis(fftshift(fft(VLD1,[],2),2),f,901);
                    set(gcf,'Name','Basis functions for fat signal');
                    pause(1);
                end
            end
        else
            VLD1                   = [];
            SLD1                   = [];
            ULD1                   = [];
        end
    else
        VLD1                       = params.VLD1;
        if isfield(params,'SLD1')
            SLD1                   = params.SLD1;
        else
            SLD1                   = ones(size(VLD1,1),1);
        end
        if isfield(params,'ULD1')
            ULD1                   = params.ULD1;
        else
            ULD1                   = ones(NyD1*NxD1,size(VLD1,1));   
        end        
    end
    output.VLD1                    = VLD1;
    output.SLD1                    = SLD1;
    output.ULD1                    = ULD1;

    % metabolite basis functions
    if ~isfield(params,'VMD1')
        if ~isempty(params.xtD1WrLr) && (params.RM~=0)
            disp('Estimating temporal basis functions of metabolite signals ...');
            tic;
            selParams              = params.selParams;
            selParams.signalType   = 'meta';
            selParams.NtEcho       = params.csiInfo.NtEcho;
            [UMD1,SMD1,VMD1,params.RM] ...
                                   = estVt(params.xtD1WrLr,params.B0MapD1Met,params.tD1,params.metaMaskD1,...
                                     NosD1, params.ND1, params.optWat, params.optLip, params.optMeta, ...
                                     params.RM, selParams, params.verbose,noiseS);
            temp                   = toc;
            if params.debug_on >0
                disp(['Estimation finished in ',num2str(temp),' s']);
                if params.debug_on == 2
                    f              = linspace(-1/params.dt/2,1/params.dt/2,NtD1);
                    show_basis(fftshift(fft(VMD1,[],2),2),f,902);
                    set(gcf,'Name','Basis functions for metabolite signal');
                    pause(1);
                end
            end
        else
            VMD1                   = [];
            SMD1                   = [];
            UMD1                   = [];
        end
    else
        VMD1                       = params.VMD1;
        if isfield(params,'SMD1')
            SMD1                   = params.SMD1;
        else
            SMD1                   = ones(size(VMD1,1),1);
        end
        if isfield(params,'UMD1')
            UMD1                   = params.UMD1;
        else
            UMD1                   = ones(NyD1*NxD1,size(VMD1,1));   
        end        
    end
    output.VMD1                    = VMD1;
    output.SMD1                    = SMD1;
    output.UMD1                    = UMD1;

    % only use water and fat basis functions
    if isfield(params,'noMeta') && (params.noMeta == 1)
        VMD1                       = [];
        SMD1                       = [];
        UMD1                       = [];
        params.RM                  = 0;
    end

    % interpolation of the temporal basis functions for D2 data, if necessary (tD2~=tD1)
    [VWD2,SWD2,UWD2]               = interpV(VWD1,SWD1,UWD1,params.tD1,params.tD2,params.InterpSvdOn);
    [VLD2,SLD2,ULD2]               = interpV(VLD1,SLD1,ULD1,params.tD1,params.tD2,params.InterpSvdOn);
    [VMD2,SMD2,UMD2]               = interpV(VMD1,SMD1,UMD1,params.tD1,params.tD2,params.InterpSvdOn);

    %% parameter preparation
    designPara                     = [];
    designPara.NHr                 = [NyD2,NxD2];
    designPara.lambdaWat           = params.lambdaWat;
    designPara.lambdaLip           = params.lambdaLip;
    designPara.lambdaMet           = params.lambdaMet;
    designPara.lambdaDfWat         = params.lambdaDfWat;
    designPara.lambdaDfLip         = params.lambdaDfLip;
    designPara.lambdaDfMet         = params.lambdaDfMet;
    designPara.specMask            = params.specMask;

    % shift masks to fit a forward model
    % bstrchanged: Do NOT shift, because my FFT-operator does not produce shifted images
    designPara.waterMask           = prepMask(params.waterMaskD2,VWD2);
    designPara.lipMask             = prepMask(params.lipMaskD2,VLD2);
    designPara.metaMask            = prepMask(params.metaMaskD2,VMD2);

    % phase term caused by B0 inhomogeneity for D2
    tD2Tr                          = reshape(params.tD2(1:params.NtTr),1,params.NtTr);
    designPara.PsiWat              = prepPsi(params.B0MapD2Wat,tD2Tr,0,VWD1);
    if(numel(designPara.PsiWat) == prod([params.Operators.OutDataSize(1:3) params.NtTr]))   % If no water/lipids/metabos are present, the resulting size can be 1x1
        designPara.PsiWat = reshape(designPara.PsiWat,[params.Operators.OutDataSize(1:3) params.NtTr]);
    end
    designPara.PsiLip              = prepPsi(params.B0MapD2Lip,tD2Tr,0,VLD1);
    if(numel(designPara.PsiLip) == prod([params.Operators.OutDataSize(1:3) params.NtTr])) 
        designPara.PsiLip              = reshape(designPara.PsiLip,[params.Operators.OutDataSize(1:3) params.NtTr]);        
    end
    designPara.PsiMet              = prepPsi(params.B0MapD2Met,tD2Tr,0,VMD1);
    if(numel(designPara.PsiMet) == prod([params.Operators.OutDataSize(1:3) params.NtTr])) 
        designPara.PsiMet              = reshape(designPara.PsiMet,[params.Operators.OutDataSize(1:3) params.NtTr]);
    end
    
    % for smoothness penalty
    designPara.PsiWat_1            = prepPsi(params.B0MapD2Wat,params.tD1,0,VWD1);
    designPara.PsiLip_1            = prepPsi(params.B0MapD2Lip,params.tD1,0,VLD1);
    designPara.PsiMet_1            = prepPsi(params.B0MapD2Met,params.tD1,0,VMD1);

    % BigLambda, x0, designPara.Vtx
    BigLambda                      = [];
    x0                             = [];
    [designPara.VtWat,BigLambda,x0]= prepDesignParaVbx(BigLambda,x0,VWD2,UWD2,SWD2,params.NtTr,NyD1,NxD1,NyD2,NxD2);
    [designPara.VtLip,BigLambda,x0]= prepDesignParaVbx(BigLambda,x0,VLD2,ULD2,SLD2,params.NtTr,NyD1,NxD1,NyD2,NxD2);
    [designPara.VtMet,BigLambda,x0]= prepDesignParaVbx(BigLambda,x0,VMD2,UMD2,SMD2,params.NtTr,NyD1,NxD1,NyD2,NxD2);
    % for smoothness penalty
    designPara.VtWat_1             = VWD1;
    designPara.VtLip_1             = VLD1;
    designPara.VtMet_1             = VMD1;
    
    if params.NUWeights
        designPara.BigLambda       = 1./sqrt(BigLambda./max(BigLambda(:)));
    else
        designPara.BigLambda       = 1;
    end

    % prepare sample index
%     sampleMask                     = zeros([NyD2In,NxD2In,length(params.tD2)]);
%     sampleMask(params.sampleIndex) = 1;
%     sampleMask                     = sampleMask(:,:,1:params.NtTr);
%     sampleIndex                    = find(sampleMask==1);
    sampleIndex = params.sampleIndex;
%     sampleMaskShift                = ifftshift(ifftshift(sampleMask,1),2);
%     sampleIndexShift               = find(sampleMaskShift==1);
    sampleIndexShift               = sampleIndex;           % bstrchanged: My FFT-operator does not produce shifted images.
    designPara.sampleIndex         = sampleIndexShift;
    designPara.Operators = params.Operators;

    % prepare data
%     temp                           = zeros([NyD2In,NxD2In,length(params.tD2)]);
    td = d(:);
%     temp(params.sampleIndex)       = d(:);
%     temp                           = temp(:,:,1:params.NtTr);
%     temp                           = ifftshift(ifftshift(temp,1),2);  % bstrchanged: My FFT-operator does not produce shifted images.
%     td                             = temp(sampleIndexShift);
    b                              = subUAOp(td(:), 'transp', designPara);

    %% Solve the resulting optimization problem
    if params.zeroInitial
        x0                         = zeros(NyD2*NxD2*(params.RW+params.RL+params.RM),1);
    end
    tic;
    [xSol, flag, RELRES, ITER, RESVEC] ...
                                   = pcg(@(x)subUProjOp( x, designPara ), b, params.tolCG, ...
                                         params.MaxIterCG,[],[],x0);
    tElapsed = toc;
    disp(['flag = ', num2str(flag), ', RELRES = ', num2str(RELRES), ', ITER =',num2str(ITER)]);
    disp(['Calculation time = ', num2str(tElapsed), ' s']);
    if params.debug_on > 0
        if params.debug_on == 2
            figure;
            plot(mag2db(RESVEC));
            title('Residual error vs # of iterations');
        end

        % calculate noise standard deviation
        temp                       = d(end-30:end);
        noiseStd                   = std(temp(:));
        noiseEnergy                = noiseStd^2*length(d(:));
        sigEnergy                  = norm(d(:))^2;
        disp(['Noise energy: ',num2str(noiseEnergy)]);
        disp(['Signal energy: ',num2str(sigEnergy)]);

        % calculate data consistency
        tdSol                      = subUAOp(xSol,'notransp',designPara);
        disp(['Data consistency error: ', num2str(norm(tdSol(:)-td(:))^2)]);

        % calculate regularization term
        % tBigLambda                 = diag(designPara.BigLambda);
        % txw                        = reshape(xSol,[],size(tBigLambda,1))*tBigLambda;
        % treg                       = designPara.lambda*norm(txw(:))^2;
        % disp(['Regularization term: ',num2str(treg)]);
    end

    %% Prepare output
    output.VLD2                    = VLD2;
    output.VMD2                    = VMD2;
    output.VWD2                    = VWD2;
    output.xSol                    = xSol;

    % re-prepare tD2, phase and basis functions
    designPara.PsiWat              = prepPsi(params.B0MapD2Wat,params.tD2,0,VWD2);
    designPara.PsiLip              = prepPsi(params.B0MapD2Lip,params.tD2,0,VLD2);
    designPara.PsiMet              = prepPsi(params.B0MapD2Met,params.tD2,0,VMD2);
    designPara.VtWat               = VWD2;
    designPara.VtLip               = VLD2;
    designPara.VtMet               = VMD2;
    designPara.sampleIndex         = params.sampleIndexShift;
    % calcuate dWL and dWrLr
    [dWL,dWrLr,output.xtWL]        = prepOutput(xSol,'reconWL',designPara,params,NyD2,NxD2,NtD2,d);
%     output.xtWL                    = fftshift(fftshift(output.xtWL,1),2);
    % calculate other output
    [output.dW,output.dWr]         = prepOutput(xSol,'reconW',designPara,params,NyD2,NxD2,NtD2,d);
    [output.dL,output.dLr]         = prepOutput(xSol,'reconL',designPara,params,NyD2,NxD2,NtD2,d);
    [output.dM,output.dMr]         = prepOutput(xSol,'reconM',designPara,params,NyD2,NxD2,NtD2,d);

    % reprepare tD1
    designPara.PsiWat              = prepPsi(params.B0MapD2Wat,params.tD1,0,VWD1);
    designPara.PsiLip              = prepPsi(params.B0MapD2Lip,params.tD1,0,VLD1);
    designPara.PsiMet              = prepPsi(params.B0MapD2Met,params.tD1,0,VMD1);
    designPara.VtWat               = VWD1;
    designPara.VtLip               = VLD1;
    designPara.VtMet               = VMD1;
    designPara.sampleIndex         = ones([NyD2,NxD2,NtD1]);
    % calculate other output
    [~,~,output.xtW]               = prepOutput(xSol,'reconW',designPara,params,NyD2,NxD2,NtD2,d);
    [~,~,output.xtL]               = prepOutput(xSol,'reconL',designPara,params,NyD2,NxD2,NtD2,d);
    [~,~,output.xtM]               = prepOutput(xSol,'reconM',designPara,params,NyD2,NxD2,NtD2,d);
%     output.xtW                     = fftshift(fftshift(output.xtW,1),2);
%     output.xtL                     = fftshift(fftshift(output.xtL,1),2);
%     output.xtM                     = fftshift(fftshift(output.xtM,1),2);



    if(params.ShowSpectra_flag)
        fabsorreal = @real;
        figure; 
        Rows = ceil(sqrt(numel(params.ShowVoxels)));
        Cols = ceil(numel(params.ShowVoxels)./Rows);
        for CurVox = 1:numel(params.ShowVoxels)
            subplot(Rows,Cols,CurVox)
            if(numel(params.ShowVoxels{CurVox}) < 3)
                params.ShowVoxels{CurVox}(3:4) = 1;
            end
            if(numel(params.ShowVoxels{CurVox}) < 4)
                params.ShowVoxels{CurVox}(4) = 1;
            end
%             UndersampledDt = params.D1RecoPar.Dwelltimes(1)/10^9*round(params.D1RecoPar.vecSize/size(d,5));
%             chemy = compute_chemshift_vector_1_2(params.D1RecoPar.LarmorFreq,UndersampledDt,size(d,5));
%             plot(chemy,squeeze(feval(fabsorreal,fftshift(fft(d(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),params.ShowVoxels{CurVox}(4),:,1))))))
%             hold on
%             plot(chemy,squeeze(feval(fabsorreal,fftshift(fft(dWL(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),params.ShowVoxels{CurVox}(4),:,1))))),'r')
%             plot(chemy,squeeze(feval(fabsorreal,fftshift(fft(d(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),params.ShowVoxels{CurVox}(4),:,1))) - fftshift(fft(dWL(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),params.ShowVoxels{CurVox}(4),:,1))))),'g')
%             hold off
%             title(sprintf('(%d,%d,%d,%d)', params.ShowVoxels{CurVox}(1), params.ShowVoxels{CurVox}(2), params.ShowVoxels{CurVox}(3), params.ShowVoxels{CurVox}(4)))
            
            tmpIn = zeros([size_MultiDims(d,1:3) size(output.xtWL,4)]);
            tmpIn(params.Operators.SamplingOperator(:)) = d;
            tmpOut = zeros([size_MultiDims(d,1:3) size(output.xtWL,4)]);
            tmpOut(params.Operators.SamplingOperator(:)) = dWrLr;
            UndersampledDt = params.D1RecoPar.Dwelltimes(1)/10^9;
            chemy = compute_chemshift_vector_1_2(params.D1RecoPar.LarmorFreq,UndersampledDt,size(tmpOut,4));
            plot(chemy,squeeze(feval(fabsorreal,fftshift(fft(tmpIn(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),:,params.ShowVoxels{CurVox}(4),1))))))
            hold on
            plot(chemy,squeeze(feval(fabsorreal,fftshift(fft(tmpIn(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),:,params.ShowVoxels{CurVox}(4),1))) - fftshift(fft(tmpOut(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),:,params.ShowVoxels{CurVox}(4),1))))),'r')
            plot(chemy,squeeze(feval(fabsorreal,fftshift(fft(tmpOut(params.ShowVoxels{CurVox}(1),params.ShowVoxels{CurVox}(2),params.ShowVoxels{CurVox}(3),:,params.ShowVoxels{CurVox}(4),1))))),'g')
            hold off
            title(sprintf('(%d,%d,%d,%d)', params.ShowVoxels{CurVox}(1), params.ShowVoxels{CurVox}(2), params.ShowVoxels{CurVox}(3), params.ShowVoxels{CurVox}(4)))
        
            
        end
    end




end

%% Subfunction: Calculate phase caused by field inhomogeneity
function Psi                       = prepPsi(B0Map,t,shift_on,VD)
    if ~isempty(VD)
        Nt                         = length(t);
        t                          = reshape(t, 1, Nt);
        if ~isempty(B0Map)
            if shift_on
                B0Map              = ifftshift(ifftshift(B0Map,1),2); 
            end
            temp                   = reshape(kron(B0Map(:),t),[],Nt);
            Psi                    = exp(1i*2*pi*temp);
        else
            Psi                    = 1;
        end
    else
        Psi                        = 1;
    end
end

%% Subfunction: Extract temporal basis functions
function [UD1,SD1,VD1,R]           = estVt(xtD1,B0Map,t,mask,NosD1,ND1,optWat,optLip,optMeta,...
                                           R,selParams,verbose, noiseS)
    NyD1                           = ND1(1);
    NxD1                           = ND1(2);
    NtD1                           = ND1(3);
    %  B0 correction
    PsiD1                          = prepPsi(B0Map,t,0,1);
    if length(PsiD1) > 1
        PsiD1                      = reshape(PsiD1,round([NyD1*NosD1,NxD1*NosD1,NtD1]));
    end
    xtD1Bc                         = mcprec2d(xtD1,PsiD1,NosD1,1);

    % extract signal of interest
    tmask                          = logical(imresize(logical(mask),round([NyD1,NxD1])));
    if ~isempty(optWat) || ~isempty(optLip) || ~isempty(optMeta)
        % HSVD
        xtD1Bc_1                   = nsRm(xtD1Bc, tmask, [], selParams, optWat, optLip, optMeta);
    else
        xtD1Bc_1                   = xtD1Bc.*repmat(tmask,[1,1,NtD1]);
    end
            
    % SVD
    temp                           = reshape(xtD1Bc_1,[round(NyD1*NxD1),NtD1]);
    [UD1,SD1,tVD1]                 = svd(temp,'econ');
    if verbose > 0    
        figure;
        plot(mag2db(diag(SD1)./noiseS(1)),'r'); hold on;
        plot(mag2db(noiseS./noiseS(1)),'k');
        xlim([0,100]);
        title('Normalized singular value');
        tR                         = input('Choose the number of basis functions:');
        if ~isempty(tR)
            R                      = tR;
        end
    end
    if(numel(find(SD1 ~= 0)) == 0)
        R = 0; SD1 = []; UD1 = []; VD1 = [];
    else
        SD1                            = SD1(1:R,1:R);
        UD1                            = UD1(:,1:R)*SD1;
        VD1                            = tVD1(:,1:R)';
    end
end

%% Subfunction: Interpolation
function [VD2,SD2,UD2]             = interpV(VD1,SD1,UD1,tD1,tD2,svd_on)
    if ~isempty(VD1)
        if (length(tD1)~=length(tD2)) | norm(tD1-tD2) ~= 0 % if tD1 ~= tD2, do interpolation
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
end      

%% Subfunction: Prepare mask
function maskRep                   = prepMask(mask,VD)
    if ~isempty(VD)
%         maskSh                     = ifftshift(ifftshift(mask,1),2);
        maskSh                     = mask;      % bstrchanged: My Fourier operator does not produce fftshifted images, so the masks should not be fftshifted!
        maskRep                    = repmat(maskSh(:),1,size(VD,1));
    else
        maskRep                    = 1;
    end
end              

%% Subfunction: Prepare V, Biglambda, and x0 of designPara
function [VDTr,BigLambda,x0]       = prepDesignParaVbx(BigLambda,x0,VD,UD,SD,NtTr,NyD1,NxD1,NyD2,NxD2)
    if ~isempty(VD)
        VDTr                       = VD(:,1:NtTr);
        BigLambda                  = [BigLambda;diag(SD)];
        % zero-padding
        for ind = 1:size(VD,1)
            temp                   = reshape(UD(:,ind),[NyD1,NxD1]);
            tempZp                 = ifftshift(ifftshift(mzprec2d(temp,[NyD2,NxD2]),1),2);
            x0                     = [x0;tempZp(:)];
        end
    else
        VDTr                       = [];
    end
end

%% Subfunction to prepare output
function [dx,dxr, xt]              = prepOutput(xSol,opSign,designPara,params,NyD2,NxD2,NtD2,d)
    [dx2,xt]                       = subUAOp(xSol, opSign, designPara);
    temp                           = zeros([params.Operators.InDataSize(1:4),NtD2]);
    temp(designPara.sampleIndex)   = dx2(:);
%     temp                           = fftshift(fftshift(temp,1),2);        % bstr: Dont need fftshifts
    dx                             = reshape(temp(params.sampleIndex),size(d));
    dxr                            = d-dx;      % bstr: WHY IS HERE NOT DIRECTLY dx2 TAKEN?!?!?!
end
