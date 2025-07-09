function PlotAllMaps(toPlot,X,Color,save_info)
if_save = nargin > 3;

if ~if_save
    save_info.name = '';
end

if isempty(Color)
    m = prctile(abs(toPlot),95);
    Color.ColorAxis = m * [-1 1];
    Color.cm = cbrewer('div','RdBu',64);
    Color.cm = Color.cm(end:-1:1,:);
    Color.alpha = 0;
end

PlotAllSlicesBasic(toPlot,X.param,{regexprep(save_info.name,'_',' ')},Color)
if  if_save && ~isempty(save_info.ext)
    saveas(gcf,mkpdir([save_info.path '/' save_info.name '.' save_info.ext]));
end

figure;
PlotTopViewBasic(toPlot,X.param,X.Anat,Color);
title(regexprep(save_info.name,'_',' '))
if  if_save && ~isempty(save_info.ext)
    saveas(gcf,mkpdir([save_info.path '/Topview_' save_info.name '.' save_info.ext]));
end

end