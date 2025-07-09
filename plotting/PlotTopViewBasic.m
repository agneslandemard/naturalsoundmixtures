function PlotTopViewBasic(toPlot,param,Anat,Color)

toPlot = Pixs2Mat(toPlot,param.msk);
func = @nanmedian;
if nargin > 3 && isfield(Color,'alpha') && Color.alpha
    nb_vx = (squeeze(sum(~isnan(toPlot),1)) - 3) ./ (prctile(mat2vec(sum(~isnan(toPlot),1)),95) - 3);
    alpha =  min(nb_vx,1);
    alpha = max(alpha, 0);
else
    nb_vx = squeeze(sum(~isnan(toPlot),1));
    alpha = nb_vx > 3;
end

imagesc(squeeze(func(toPlot,1)),'AlphaData',alpha)
if nargin < 4
    bd = prctile(abs(toPlot(:)),99);
    clim([-bd bd])
    colormap('redblue')
else
    clim(Color.ColorAxis)
    colormap(Color.cm)
end
hold on
vals = [prctile(mat2vec(snm(Anat,1)),90) prctile(mat2vec(snm(Anat,1)),95)];
contour(snm(Anat,1),vals)
set(gca,'DataAspectRatio',[1 4 1])
set(gca,'XTick',[],'YTick',[])
cbh = colorbar;
if isfield(Color,'colorbar_ticks')
    cbh.TickLabels = Color.colorbar_ticks;
end

end
