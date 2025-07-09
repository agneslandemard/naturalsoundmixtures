function PlotMapWithAnat(toPlot,Anat,Color)

color_range=[-Inf linspace(Color.ColorAxis(1),Color.ColorAxis(end),size(Color.cm,1)-1) Inf];

% overlay color info on anat
ma = prctile(Anat(:),98);
mi = prctile(Anat(:),2);
norm_Anat = (Anat - mi)./(ma-mi);
X = repmat(norm_Anat,[1 1 3]);
ns = toPlot(:);
ns_colors = nan(length(ns),3);
for x = 1:length(ns)
    for n = 1:length(color_range)-1
        if ns(x) >= color_range(n) && ns(x) < color_range(n+1)
            ns_colors(x,:) = Color.cm(n,:);
            break;
        end
    end
end

NanPx = isnan(sum(ns_colors,2));
PixtoPlot = ~isnan(toPlot);
PixtoPlot_NN = find(PixtoPlot);
PixtoPlotNew = false(size(PixtoPlot));
PixtoPlotNew(PixtoPlot_NN) = 1;

X(repmat(PixtoPlotNew,[1 1 3])) = ns_colors(~NanPx,:);

imagesc(X)
axis equal tight
set(gca,'XTick',[],'YTick',[])
colormap(Color.cm)

end