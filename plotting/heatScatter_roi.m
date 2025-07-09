function dens = heatScatter_roi(x,y,c,P)
N = 100;
if nargin < 4
    P.n_bins = 40;
elseif ~isfield(P,'n_bins')
    P.n_bins = 40;
end
if ~isfield(P,'min')
    P.min = 1;
end
if ~isfield(P,'alpha')
    P.alpha = 0;
else
    if ~isfield(P,'alpha_threshold')
        P.alpha_threshold = 10;
    end

end
n_rois = size(P.roi_colors,1);

if isempty(c)
    d = 1;
    n_vxs = length(x);
    M = 10* n_vxs / P.n_bins.^2;
else
    d = 0;
    if isfield(P,'clim')
        M = P.clim;
    else
        
        M = nanmax(c(:));
        if M > 0.2 
             M = round(M,1);
        end
    end
end

bins1 = min(x) :(max(x)-min(x))/(P.n_bins-1) :max(x);
bins2 = min(y) :(max(y)-min(y))/(P.n_bins-1) :max(y);

% Color.cm = cbrewer('div','RdBu',N+1);
% Color.cm = Color.cm(end:-1:1,:);
% n = 0.8; K = floor((N+1)/4);
% 
% b = Color.cm(K,:);
% r = Color.cm(end-K,:);
%cscale = cat(1,Color.cm(1:K,:),linspaceNDim(b,n*[1 1 1],(N/2-K))',linspaceNDim(n*[1 1 1],r,(N/2-K))',Color.cm(end-K:end,:));

hold all
cbins = linspace(-M,M,N);
for b1 = 1:P.n_bins-1
    grp1 = (x > bins1(b1)).*(x < bins1(b1+1));
    for b2 = 1:P.n_bins-1
        grp2 = (y > bins2(b2)).*(y < bins2(b2+1));
        vxs_bin = logical(grp1.*grp2);
        
        if d
            clr = length(find(vxs_bin));
            dens.density(b1,b2) = length(find((vxs_bin)));
        
        else
            clr = nanmean(c(vxs_bin));
            dens.density(b1,b2) = length(find(~isnan(c(vxs_bin))));
        
        end
    end
end
b0 = min(bins1(1),bins2(1));
be = max(bins1(end),bins2(end));

if P.alpha
   
  %  patch([b0 be be b0],[b0 b0 be be],0.9*[1 1 1],'EdgeAlpha',0);

end              
for b1 = 1:P.n_bins-1
    grp1 = (x > bins1(b1)).*(x < bins1(b1+1));
    for b2 = 1:P.n_bins-1
        grp2 = (y > bins2(b2)).*(y < bins2(b2+1));
        vxs_bin = logical(grp1.*grp2);
        
        if d
            clr = length(find(vxs_bin));

        else
            tmp_c = c(vxs_bin);
            tmp_c(isnan(tmp_c))=[];
            n_vx_by_roi =  histcounts(tmp_c,1:(n_rois + 1),'Normalization','probability');
            
            clr = n_vx_by_roi*P.roi_colors;
           

        end
        if P.alpha
            a = min(1,dens.density(b1,b2)./P.alpha_threshold);
        else
            a = 1;
        end
        
        if any(clr)
         
            if length(find(vxs_bin)) >= P.min
                
                patch([bins1(b1) bins1(b1+1) bins1(b1+1) bins1(b1)],[bins2(b2) bins2(b2) bins2(b2+1) bins2(b2+1)],clr,'EdgeAlpha',0,'FaceAlpha',a);
                
            end
        end
    end
end
l = line([b0 be],[b0 be]);
l.LineStyle = ':';
axis equal tight

hline(0,'k:'); vline(0,'k:')
ylim([b0 be]); xlim([b0 be])
colormap(P.roi_colors)
cbb = colorbar;


hold on;


dens.bins1 = bins1;
dens.bins2 = bins2;
end