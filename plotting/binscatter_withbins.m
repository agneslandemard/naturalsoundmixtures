function binscatter_withbins(x,y,bin_width,on_fig)

if nargin < 4
    on_fig = 0;
end
if nargin < 3
    bin_width = 0.1;
elseif isempty(bin_width)
    bin_width = 0.1;
end

if ~iscell(y)
    y1{1} = y;
else 
    y1 = y;
end
if ~on_fig
    figure;
end
binscatter(x,y1{1},40)
hold all
line(ylim,ylim)
vline(0,'k'); hline(0,'k')
thresh_range = round(prctile(x,1),1): bin_width : round(prctile(x,99),1);
save_tmp = nan(length(thresh_range)-1,length(y1)+1);
for t = 1:length(thresh_range)-1
    vxs = logical((x >= thresh_range(t)).*(x < thresh_range(t+1)));
    save_tmp(t,1) = nanmedian(x(vxs));
    for k = 1 : length(y1)
        scatter(nanmedian(x(vxs)),nanmedian(y1{k}(vxs)),60,'r*')
        save_tmp(t,k+1) = nanmedian(y1{k}(vxs));
    end
   
end

for k = 1 : length(y1)
    p1(k) = plot(save_tmp(:,1),save_tmp(:,k+1),'--');
end
%legend(p1)
end