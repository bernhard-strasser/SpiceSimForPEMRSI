function VisualizeTGV( Vol_Before, Vol_After ,name_fig)
%VISUALIZETGV Summary of this function goes here
%   Detailed explanation goes here
warning('off','MATLAB:DELETE:FileNotFound');
nDims=ndims(Vol_After);
SizeData=size(Vol_After);
s=sprintf('%s_ComparisonTGV.ps',name_fig);
delete(s); 
ImSiC=2048;ImSiR=2048;
    for k=1:SizeData(end)
        figs(k)=figure('visible', 'off','Position', [100, 200, ImSiC,ImSiR]); 
       if nDims ==4
           image2plot=(sum(Vol_Before(:,:,:,k),3));
       else
           image2plot=(Vol_Before(:,:,k));
       end
        lim=quantile(abs(image2plot(:)),0.99);%max(abs(image2plot(:)));
        subplot(2,2,1),imagesc(abs(image2plot),[0 lim]),colorbar, title('Initial Spatial Component ABS');  
        axis('off');
        colormap default 
        subplot(2,2,2),imagesc(angle(image2plot)),colorbar, title('Initial Spatial Component Phase');  
        axis('off');
        colormap default 
        if nDims ==4
           image2plot=(sum(Vol_After(:,:,:,k),3));
       elseif nDims ==3
           image2plot=(Vol_After(:,:,k));
        end
        lim=quantile(abs(image2plot(:)),0.99);%max(abs(image2plot(:)));%
        subplot(2,2,3),imagesc(abs(image2plot),[0 lim]),colorbar,title('Recon Spatial Component ABS');  
        axis('off');
        colormap default
        subplot(2,2,4),imagesc(angle(image2plot)),colorbar,title('Recon Spatial Component Phase');  
        axis('off');
        colormap default    
        print(figs(k),'-bestfit', '-append', '-dpsc2', s); 
        close(figs(k))
    end




end

