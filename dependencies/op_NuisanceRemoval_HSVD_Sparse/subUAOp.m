function [y,yk] = subUAOp(x, opSign, params)
	NHr                       = params.NHr;
	NyHr                      = NHr(1);
	NxHr                      = NHr(2);
	Nt                        = max([size(params.VtWat,2),size(params.VtLip,2),size(params.VtMet,2)]);
	RWat                      = size(params.VtWat,1);
	RLip                      = size(params.VtLip,1);
	RMet                      = size(params.VtMet,1);
    Ns                        = length(params.sampleIndex(:));

	if strcmp(opSign,'notransp')
		y                     = AOp(x);
	elseif strcmp(opSign,'transp')
		y                     = AcOp(x);
    elseif strcmp(opSign,'notransp_df')
        Nt_1                  = max([size(params.VtWat_1,2),size(params.VtLip_1,2),size(params.VtMet_1,2)]);
        y                     = AOp_df(x);
    elseif strcmp(opSign,'transp_df')
        Nt_1                  = max([size(params.VtWat_1,2),size(params.VtLip_1,2),size(params.VtMet_1,2)]);
        y                     = AcOp_df(x);
	elseif strcmp(opSign,'recon')
		[y,yk]                = AOp(x);
	elseif strcmp(opSign,'reconWL') || strcmp(opSign,'reconW') ...
		   || strcmp(opSign,'reconL') || strcmp(opSign,'reconM') 
		[y,yk]                = AOpOut(x);
	else
		warning('opSign is either notransp or transp!')
	end
    
	% generate output
	function [yWLM,reconWLM]  = AOpOut(x)
		% obtain UxLip and UxMet from x
		x                     = reshape(x,NyHr*NxHr,RWat+RLip+RMet);
        % water signal
        if RWat               ~= 0
            UxWat             = params.waterMask.*x(:,1:RWat);
            temp              = reshape(params.PsiWat.*(UxWat*params.VtWat),[params.Operators.OutDataSize(1:3),Nt]);
            temp = temp .* conj(params.Operators.SensMap);
            reconW            = temp;
            
            % x to k
%             temp = reshape(temp,[prod(params.Operators.OutDataSize(1:2)) numel(temp)/prod(params.Operators.OutDataSize(1:2))]);
%             tktHr = params.Operators.sft2_Oper * temp;
%             tktHr = fftshift(fftshift(fft(fft(ifftshift(ifftshift(temp,1),2),[],1),[],2),1),2)/sqrt(NyHr)/sqrt(NxHr);
            tktHr = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(conj(temp),1),2),[],1),[],2),1),2);

            % truncation
            yW                = tktHr(params.sampleIndex);
        else
            reconW            = 0;
            yW                = 0;
        end
        % lipid signal
        if RLip               ~= 0
            UxLip             = params.lipMask.*x(:,(RWat+1):(RWat+RLip));
            temp              = reshape(params.PsiLip.*(UxLip*params.VtLip),[params.Operators.OutDataSize(1:3),Nt]);
            temp = temp .* conj(params.Operators.SensMap);
            reconL            = temp;
            
            % x to k
%             temp = reshape(temp,[prod(params.Operators.OutDataSize(1:2)) numel(temp)/prod(params.Operators.OutDataSize(1:2))]);
%             tktHr = params.Operators.sft2_Oper * temp;
%             tktHr = fftshift(fftshift(fft(fft(ifftshift(ifftshift(temp,1),2),[],1),[],2),1),2)/sqrt(NyHr)/sqrt(NxHr);
            tktHr = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(conj(temp),1),2),[],1),[],2),1),2);

            % truncation
            yL                = tktHr(params.sampleIndex);
        else
            reconL            = 0;
            yL                = 0;
        end
        % metabolite signal
        if RMet               ~= 0
            UxMet             = params.metaMask.*x(:,(RWat+RLip+1):end);
            temp              = reshape(params.PsiMet.*(UxMet*params.VtMet),[params.Operators.OutDataSize(1:3),Nt]);
            temp = temp .* conj(params.Operators.SensMap);
            reconM            = temp;
            
            % x to k
%             temp = reshape(temp,[prod(params.Operators.OutDataSize(1:2)) numel(temp)/prod(params.Operators.OutDataSize(1:2))]);
%             tktHr = params.Operators.sft2_Oper * temp;
%             tktHr = fftshift(fftshift(fft(fft(ifftshift(ifftshift(temp,1),2),[],1),[],2),1),2)/sqrt(NyHr)/sqrt(NxHr);
            tktHr = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(conj(temp),1),2),[],1),[],2),1),2);

            % truncation
            yM                = tktHr(params.sampleIndex);
        else
            reconM            = 0;
            yM                = 0;
        end

        % generate output
        switch opSign
            case 'reconW'
                yWLM          = yW;
                reconWLM      = reconW;
            case 'reconL'
                yWLM          = yL;
                reconWLM      = reconL;
            case 'reconM'
                yWLM          = yM;
                reconWLM      = reconM; % WAS: yM
            case 'reconWL'
                yWLM          = yW + yL;
                reconWLM      = reconW + reconL;
            otherwise
                yWLM          = yW + yL + yM;
                reconWLM      = reconW + reconL + reconM;
        end
	end % end of AOpOut

	% forward operator for updating U with given V
	function y                = AOp(x)
		% obtain UxLip and UxMet from x
		x                     = reshape(x,NyHr*NxHr,RWat+RLip+RMet);
        ImtHr                 = 0;
        if RWat               ~= 0
			UxWat             = params.waterMask.*x(:,1:RWat);
            ImtHr = ImtHr + params.PsiWat.*reshape(UxWat*params.VtWat,[params.Operators.OutDataSize(1:3) Nt]);
        end
        if RLip               ~= 0
            UxLip             = params.lipMask.*x(:,(RWat+1):(RWat+RLip));
            ImtHr             = ImtHr + params.PsiLip.*reshape(UxLip*params.VtLip,[params.Operators.OutDataSize(1:3) Nt]);
        end
        if RMet               ~= 0     
			UxMet             = params.metaMask.*x(:,(RWat+RLip+1):end);
            ImtHr             = ImtHr + params.PsiMet.*reshape(UxMet*params.VtMet,[params.Operators.OutDataSize(1:3) Nt]); 
        end

        % Apply sensitivity map
        ImtHr = ImtHr .* conj(params.Operators.SensMap);
        
% 		ImtHr                 = reshape(ImtHr,[NyHr*NxHr,numel(ImtHr)/NyHr/NxHr]);
		% downsampling
		% x to k
%         ktHr = params.Operators.sft2_Oper * ImtHr;
%         ktHr = fftshift(fftshift(fft(fft(ifftshift(ifftshift(ImtHr,1),2),[],1),[],2),1),2)/sqrt(NyHr)/sqrt(NxHr);
        ktHr = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(conj(ImtHr),1),2),[],1),[],2),1),2);

% 		ktHr                  = fft(fft(ImtHr,[],1),[],2)/sqrt(NyHr)/sqrt(NxHr);
		% truncation
		y                     = ktHr(params.sampleIndex);
	end % end of AOp

	% Adjoint of the forward operator for updating U with given V
	function y                = AcOp(x)
		% zeropadding
		ktHr                  = zeros([params.Operators.InDataSize(1:4) Nt]);
		ktHr(params.sampleIndex)    = x;
        ktHr = reshape(ktHr,[params.Operators.OutDataSize(1:3) Nt]);

%         ImtHr = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(ktHr,1),2),[],1),[],2),1),2)*sqrt(NyHr)*sqrt(NxHr);
        ImtHr = fftshift(fftshift(conj(fft(fft(ifftshift(ifftshift(ktHr,1),2),[],1),[],2)),1),2);

        % Coil Combination
        ImtHr = sum(ImtHr .* params.Operators.SensMap,5);
                
		% estimate UxLip and UxMet
		y                     = [];
        if RWat               ~= 0		
            % Undo B0 Effect
            UxWat = conj(params.PsiWat).*ImtHr;
            % Calculate U
            UxWat = reshape(UxWat,[numel(UxWat)/size(params.VtWat,2) size(params.VtWat,2)]);
            UxWat = UxWat * params.VtWat';
            % Mask U
            UxWat = params.waterMask .* UxWat;
            y                 = [y;UxWat(:)];
        end
        if RLip               ~= 0		
            % Undo B0 Effect
            UxLip = conj(params.PsiLip).*ImtHr;
            % Calculate U
            UxLip = reshape(UxLip,[numel(UxLip)/size(params.VtLip,2) size(params.VtLip,2)]);
            UxLip = UxLip * params.VtLip';
            % Mask U
            UxLip = params.lipMask .* UxLip;
            y                 = [y;UxLip(:)];
        end
		if RMet               ~= 0		
            % Undo B0 Effect
            UxMet = conj(params.PsiMet).*ImtHr;
            % Calculate U
            UxMet = reshape(UxMet,[numel(UxMet)/size(params.VtMet,2) size(params.VtMet,2)]);
            UxMet = UxMet * params.VtMet';
            % Mask U
            UxMet = params.metaMask .* UxMet;
            y                 = [y;UxMet(:)];
        end
	end % end of AcOp

    % forward operator for updating U with given V + difference along f
    function y                = AOp_df(x)
        % data consistence part
        y                     = AOp(x);

        % calculate diff_f
        x                     = reshape(x,NyHr*NxHr,RWat+RLip+RMet);
        ImtHr_diff            = [];
        if RWat ~= 0 & params.lambdaDfWat ~= 0
            UxWat             = params.waterMask.*x(:,1:RWat);
            temp              = params.PsiWat_1.*(UxWat*params.VtWat_1);
            temp              = params.lambdaDfWat*params.specMask.*diff(mfft1d(temp,[],2),1,2);
            ImtHr_diff        = [ImtHr_diff;temp(:)];
        end
        if RLip ~= 0 & params.lambdaDfLip ~= 0
            UxLip             = params.lipMask.*x(:,(RWat+1):(RWat+RLip));
            temp              = params.PsiLip_1.*(UxLip*params.VtLip_1);
            temp              = params.lambdaDfLip*params.specMask.*diff(mfft1d(temp,[],2),1,2);
            ImtHr_diff        = [ImtHr_diff;temp(:)];
        end
        if RMet ~= 0 & params.lambdaDfMet ~= 0
            UxMet             = params.metaMask.*x(:,(RWat+RLip+1):end);
            temp              = params.PsiMet_1.*(UxMet*params.VtMet_1); 
            temp              = params.lambdaDfMet*params.specMask.*diff(mfft1d(temp,[],2),1,2);
            ImtHr_diff        = [ImtHr_diff;temp(:)];
        end
        if ~isempty(ImtHr_diff)
            y                 = [y;ImtHr_diff];
        end
    end % end of AOp_df

    % Adjoint of the forward operator for updating U with given V + diffrence along f
    function y                = AcOp_df(x)
        % data consistence part
        tx                    = x(1:Ns);
        y                     = AcOp(tx);

        % y_df;
        y_df                  = [];
        tInd                  = Ns;
        if (RWat ~= 0) & (params.lambdaDfWat ~= 0)
            temp_df           = x((tInd+1):(tInd+NyHr*NxHr*(Nt_1-1)));
            tInd              = tInd + NyHr*NxHr*(Nt_1-1);
            temp_df           = params.lambdaDfWat*params.specMask.*reshape(temp_df,NyHr*NxHr,Nt_1-1);
            temp_ddf          = mifft1d(conjDiff(temp_df,NyHr*NxHr,Nt_1),[],2);
            UxWat             = params.waterMask.*((conj(params.PsiWat_1).*temp_ddf)*params.VtWat_1');
            y_df              = [y_df;UxWat(:)];
        elseif (RWat ~= 0) & (params.lambdaDfWat == 0)
            y_df              = [y_df;zeros(NyHr*NxHr*RWat,1)];
        end
        if (RLip ~= 0) & (params.lambdaDfLip ~= 0)
            temp_df           = x((tInd+1):(tInd+NyHr*NxHr*(Nt_1-1)));
            tInd              = tInd + NyHr*NxHr*(Nt_1-1);
            temp_df           = params.lambdaDfLip*params.specMask.*reshape(temp_df,NyHr*NxHr,Nt_1-1);
            temp_ddf          = mifft1d(conjDiff(temp_df,NyHr*NxHr,Nt_1),[],2);
            UxLip             = params.lipMask.*((conj(params.PsiLip_1).*temp_ddf)*params.VtLip_1');
            y_df              = [y_df;UxLip(:)];
        elseif (RLip ~= 0) & (params.lambdaDfLip == 0)
            y_df              = [y_df;zeros(NyHr*NxHr*RLip,1)];
        end
        if (RMet ~= 0) & (params.lambdaDfMet ~= 0)
            temp_df           = x((tInd+1):(tInd+NyHr*NxHr*(Nt_1-1)));
            tInd              = tInd + NyHr*NxHr*(Nt_1-1);
            temp_df           = params.lambdaDfMet*params.specMask.*reshape(temp_df,NyHr*NxHr,Nt_1-1);
            temp_ddf          = mifft1d(conjDiff(temp_df,NyHr*NxHr,Nt_1),[],2);
            UxMet             = params.metaMask.*((conj(params.PsiMet_1).*temp_ddf)*params.VtMet_1');
            y_df              = [y_df;UxMet(:)];
        elseif (RMet ~= 0) & (params.lambdaDfMet == 0)
            y_df              = [y_df;zeros(NyHr*NxHr*RMet,1)];
        end

        if ~isempty(y_df)
            y                 = y + y_df;
        end

    end % end of AcOp

end % end of subUAOp

function y = conjDiff(x,N1,N2)
    y                         = zeros(N1,N2);
    y(:,1)                    = -x(:,1);
    y(:,2:(end-1))            = -diff(x,1,2);
    y(:,end)                  = x(:,end);
end
