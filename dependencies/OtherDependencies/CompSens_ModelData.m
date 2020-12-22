function Output = CompSens_ModelData(Transj,Input,Operators)
%
% Spice_ModelSpiralData: Function to model spiral data given U*Sigma, or the opposite
%
% This function was written by Bernhard Strasser, December 2018.
%
% In case of Transj = 'NoTransj', models spiral data given V(=Temporal Basis) and UTS = U*Sigma, where [U,Sigma,V] = svd(Rho),  
% with Rho ... Cartesian image domain MRSI data. So the function calculates SpiralData = A(UTS,V).
% In case of Transj = 'Transj', models UTS from SpiralData and V. 
% So the function calculates UTS = A^H(SpiralData,V), where A^H is the adjoint operator of A.
% However, because of convergence reasons, A^H includes a densitiy compensation for spiral data, while A does not (so strictly speaking
% A^H is not the adjoint of A...).
% This function can be used for an iterative reconstruction in the SPICE framework.
%
% 
%
%
% Output = Spice_SynthesizeMeasData(U,V,Operators.B0CorrMat_Spec,Operators.SamplingOperator)
%
% Input: 
% -         U                    ...     
% Output:
% -         Output                         ...     



% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!

% Further remarks: 




%% 0. Preparations

if(~exist('Transj','var') || isempty(Transj))
   Transj = 'NoTransj'; 
end



%%


if(strcmpi(Transj,'NoTransj') || strcmpi(Transj,'notransp')) % Going from UTS to NonCart k-space data

    % Input = U * Sigma
    % Output = SamplingOperator .* OffResoEffect ( sFT_iTokSpace (((U * Sigma) * V) ° B0) ) = Spiral-k-Space-Data = ModelFunction (U * Sigma)
    
%     Output = conj(Input);
    Output = reshape(Input,Operators.OutDataSize);

    % Mask Data
    Output = Output .* Operators.Mask;
    

    
    % Apply B0-Effect
    Output = Output .* conj(Operators.B0CorrMat_Spec);

    % Apply Coil Sensitivity
    Output = Output .* conj(Operators.SensMap);
    
    
%     % Inverse NUFFT (i-Space --> NonCart k-Space)
%     Output = reshape(Output,[size(Operators.sft2_Oper,2) numel(Output)/size(Operators.sft2_Oper,2)]);
%     Output = Operators.sft2_Oper * Output;
%     Output = reshape(Output,[Operators.InDataSize]);
    Output = fftshift(fftshift(fft(fft(ifftshift(ifftshift(Output,1),2),[],1),[],2),1),2)/sqrt(Operators.InDataSize(1)*Operators.InDataSize(2));
%     Output = fftshift(fftshift(fft(fft(ifftshift(ifftshift(Output,1),2),[],1),[],2),1),2);
%     Output = fftshift(fftshift(fft(fft(ifftshift(ifftshift(Output,1),2),[],1),[],2),1),2)*sqrt(Operators.InDataSize(1)*Operators.InDataSize(2));
    Output = reshape(Output,Operators.InDataSize);

    
%     dummy = Output(:,1);
%     Output = fft(ifftshift(Operators.TiltTrajMat.*fftshift(ifft(Output,[],5),5),5),[],5);
%     Output(:,1) = dummy;

    % FoVShift
    Output = Output .* conj(Operators.FoVShift);
    
    Output = Operators.SamplingOperator .* Output;
    
%     Output = reshape(Output,[numel(Output) 1]);    
    

    
    

else    % Going from NonCart k-space data to UTS

    % Input = Spiral-k-Space-Data
    % Output = (sFT_kToiSpace(DensComp(OffResoCorr(Spiral-k-Space-Data .* SamplingOperator))) ° B0) * V' = More or less Cartesian image
    
    Output = reshape(Input,Operators.InDataSize);
    
    % FoVShift
    Output = Output .* (Operators.FoVShift);
    
%     Output = reshape(Input,[size_MultiDims(Operators.DCFPreG,[1 2]) size_MultiDims(Operators.SamplingOperator,[3 4]) size(Operators.B0CorrMat_Spec,3) size(Operators.SensMap,4)]);
    Output = Operators.SamplingOperator .* Output;
    
%     dummy = Output(:,1);
%     Output = fft(ifftshift(conj(Operators.TiltTrajMat).*fftshift(ifft(Output,[],5),5),5),[],5);
%     Output(:,1) = dummy;    
    
    % NUFFT (NonCart k-Space --> i-Space)
%     Output = Output .* Operators.DCFPreG; % Apply PreGridding Density Comp
%     Output = reshape(Output,[prod(Operators.InDataSize(1:2)) prod(Operators.InDataSize(3:end))] );
% 	Output = Operators.sft2_Oper' * Output * size(Operators.sft2_Oper,2); %
%     Output = reshape(Output,[Operators.OutDataSize Operators.InDataSize(end)] );
    Output = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(Output,1),2),[],1),[],2),1),2)*sqrt(Operators.InDataSize(1)*Operators.InDataSize(2));
%     Output = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(Output,1),2),[],1),[],2),1),2);
%     Output = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(Output,1),2),[],1),[],2),1),2)/sqrt(Operators.InDataSize(1)*Operators.InDataSize(2));
    Output = reshape(Output,Operators.OutDataSize);
    
    
    
    
    % Coil Combination
    Output = sum(Output .* Operators.SensMap,5);
    
    % Undo B0-Effect
%     Output = reshape(Output,size(Operators.B0CorrMat_Spec));
    Output = Output .* Operators.B0CorrMat_Spec;

    

    

    % Mask Data
    Output = Output .* Operators.Mask;
    
%     Output = reshape(Output,[numel(Output) 1]);
    
%     Output = conj(Output);
    
    
    

    
    
    
    
    
end


