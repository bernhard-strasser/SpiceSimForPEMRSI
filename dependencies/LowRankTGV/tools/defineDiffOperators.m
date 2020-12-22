function [fhdxm,fhdym,fhdzm,fhdxp,fhdyp,fhdzp] = defineDiffOperators()
% defineDiffOperators
% Defines the necessary forward and adjoint finite difference operators
% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUTS:
%       step_size: A (1x3) vector containing the step sizes to be used for 
%                  the calculation of finite differences along each 
%                  Cartesian dimension
%--------------------------------------------------------------------------
%   OUTPUTS:
%      fhdxm: Function handle for the forward difference operator along the 
%             x-dimension
%      fhdym: Function handle for the forward difference operator along the 
%             y-dimension
%      fhdzm: Function handle for the forward difference operator along the 
%             z-dimension
%      fhdxp: Function handle for the adjoint difference operator along the 
%             x-dimension
%      fhdyp: Function handle for the adjoint difference operator along the 
%             y-dimension
%      fhdzp: Function handle for the adjoint difference operator along the 
%             z-dimension
%--------------------------------------------------------------------------
% NOTES:
% This function implements the forward finite difference operator using
% Neumann boundary conditions. In this case, for interval [a,b], the
% condition is that f'(b) = 0;
% The adjoint operator (i.e., divergence) uses Dirichlet boundary 
% conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fhdxm = @(x) dxm(x);
fhdym = @(x) dym(x);
fhdzm = @(x) dzm(x);
fhdxp = @(x) dxp(x);
fhdyp = @(x) dyp(x);
fhdzp = @(x) dzp(x);

%--------------------------------------------------------------------------
    function dx = dxm(im)
        % Forward finite differences in the x-dimension
	%[M N P] = size(u);
	%dx = cat(2,u(:,1:end-1,:),zeros(M,1,P)) - cat(2,zeros(M,1,P),u(:,1:end-1,:));        

        [Ny, ~, Nz, Next1,Next2] = size(im);
	ExtraIndcs  = repmat({':'}, [1, (ndims(im)-3)]);

        dx =  (cat(2, im(:, 1:end-1,:, ExtraIndcs{:}), zeros([Ny, 1, Nz, Next1,Next2], class(im))) - cat(2,  zeros([Ny, 1, Nz, Next1,Next2], class(im)), im(:, 1:end-1, :, ExtraIndcs{:})));
    end
%--------------------------------------------------------------------------
    function dy = dym(im)
        % Forward finite differences in the y-dimension
        
	%[M N P] = size(u);
	%dy = cat(1,u(1:end-1,:,:),zeros(1,N,P)) - cat(1,zeros(1,N,P),u(1:end-1,:));
        [~, Nx, Nz, Next1,Next2] = size(im);
        ExtraIndcs  = repmat({':'}, [1, (ndims(im)-3)]);

        dy =  (cat(1, im(1:end-1, :,:, ExtraIndcs{:}),  zeros([1, Nx, Nz, Next1,Next2], class(im))) - cat(1, zeros([1, Nx, Nz, Next1,Next2], class(im)), im(1:end-1, :,:, ExtraIndcs{:})));
    end
%--------------------------------------------------------------------------
    function dz = dzm(im)
        % Forward finite differences in the z-dimension
        
        [Ny, Nx, ~, Next1,Next2] = size(im);
        ExtraIndcs  = repmat({':'}, [1, (ndims(im)-3)]);
        dz =  (cat(3, im(:, :, 1:end-1,  ExtraIndcs{:}), zeros([Ny, Nx, 1, Next1,Next2], class(im))) - cat(3,  zeros([Ny, Nx, 1, Next1,Next2], class(im)), im(:, :, 1:end-1,  ExtraIndcs{:})));
    end
%--------------------------------------------------------------------------
    function dx = dxp(im)
        % Backward finite differences in the x-dimension
	%[M N P] = size(u); % numSamplesOnSpoke, numSpokes, nCh
	%dx = cat(2,u(:,2:end,:),u(:,end,:)) - u;
        ExtraIndcs  = repmat({':'}, [1, (ndims(im)-3)]);
        dx =  (cat(2, im(:, 2:end, :, ExtraIndcs{:}),  im(:, end, :, ExtraIndcs{:})) - im);
    end
%--------------------------------------------------------------------------
    function dy = dyp(im)
        % Backward finite differences in the y-dimension
	%[M N P] = size(u);
	%dy = cat(1,u(2:end,:,:),u(end,:,:)) - u;

        ExtraIndcs  = repmat({':'}, [1, (ndims(im)-3)]);
        dy =  (cat(1, im(2:end, :,:, ExtraIndcs{:}),   im(end, :,:, ExtraIndcs{:})) - im);
    end
%--------------------------------------------------------------------------
    function dz = dzp(im)
        % Backward finite differences in the y-dimension
        ExtraIndcs  = repmat({':'}, [1, (ndims(im)-3)]);
        dz =  (cat(3, im(:,:, 2:end,  ExtraIndcs{:}), im(:,:, end,  ExtraIndcs{:})) - im);
    end
%--------------------------------------------------------------------------

end
