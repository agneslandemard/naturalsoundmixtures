function PlotAllSlicesBasic(X,param,titles,Color)
% will plot weights for all slices 
% X should be voxels * feats to plot
% colorbar will be centered on 0 and adjusted to the values of X for each
% feature separately
if nargin > 3
    cm = Color.cm;
    axi = Color.ColorAxis;
    fixedaxi = axi;
else
    cm = cmap_from_name('cbrewer-blue-red'); 
end
if nargin < 3
    titles = [];
end

for ft = 1:size(X,2)
    W = Pixs2Mat(X(:,ft),param.msk);
    m = max(abs(W(:)));
    if ~exist('fixedaxi','var')
        axi = 0.8*[-m m];
    end
    figure('Position', [440 607 922 191]);
    for ii = 1:size(W,3)
        subplot(ceil(sqrt(size(W,3))),ceil(sqrt(size(W,3))),ii)
        imagesc(W(:,:,ii),'AlphaData',~isnan(W(:,:,ii)));
        axis equal tight; colormap(cm);caxis(axi);
        set(gca,'XTick',[],'YTick',[])
    end
    if ~isempty(titles)
        sgtitle(titles{ft})
    end
    
end
end
