%% Load Betas and format

species = 'ferret'; % 'ferret' or 'human'

[Betas, roi_idx, rois] = load_model_weights(species);

%% Look at weights across rois

for r = 1:length(rois)
    
    to_plot = nanmedian(Betas(roi_idx == r,:),1);
    m = prctile(mat2vec(abs(to_plot)),99);
    [fig1, fig2] = plot_in_CM_space(to_plot,['Weights ' rois{r}]);
    figure(fig1)
    clim([-m, m])
    figure(fig2)
    clim([-m, m])

    if  ~isempty(ext)
        saveas(fig1,mkpdir([figures_path 'CM_weights_'  rois{r} '_' species '.' ext]));
        saveas(fig2,mkpdir([figures_path 'CM_weights_' rois{r}  '_full_' species '.' ext]));
    end
end

%% Difference in weights across regions

for r = 2:length(rois)

    to_plot = nanmedian(Betas(roi_idx==r,:),1) - nanmedian(Betas(roi_idx==1,:),1);
    [fig1, fig2] = plot_in_CM_space(to_plot, ['Diff weights ' rois{r} ' - MEG']);
    figure(fig1)
    colormap('redblue')
    m = prctile(to_plot,99);
    clim([-m, m])
    figure(fig2)
    colormap('redblue')
    clim([-m, m])
      
    if  ~isempty(ext)
        saveas(fig1,mkpdir([figures_path 'CM_weights_diff_' rois{1} '_'  rois{r} '_full' species '.' ext]));
        saveas(fig2,mkpdir([figures_path 'CM_weights_diff_' rois{1} '_' rois{r} '_' species '.' ext]));
    end    
end

%% Differences in spectrotemporal model energy foregrounds vs backgrounds
switch species
    case 'ferret'
        load([metrics_path 'cortical_model_ferrets.mat'], 'patterns_CM');
        [vals, feats] = LoadCMparams;
        
        groups = {1,[4 5],[2 3]};
        diff_fg_bg = snm(patterns_CM(:,:,:,:,groups{1},:),2:6) -...
            snm(patterns_CM(:,:,:,:,groups{2},:),2:6);

    case 'human'
        load([metrics_path 'cortical_model_humans.mat'], 'all_CM');
        CM = all_CM.subj04;
        diff_fg_bg = snm(CM(:,:,1),2) - snm(CM(:,:,2),2);
end

[fig1, fig2] = plot_in_CM_space(diff_fg_bg);
if  ~isempty(ext)
    figure(fig1)
    colormap('redblue')
    m = prctile(diff_fg_bg,99);
    clim([-m, m])

    figure(fig2)
    colormap('redblue')
    m = prctile(diff_fg_bg,99);
    clim([-m, m])
    saveas(fig1,mkpdir([figures_path 'CM_' species '_diff_fgbg.' ext]));
    saveas(fig2,mkpdir([figures_path ['CM_' species '_diff_fgbgfull' ...
        '.'] ext]));
end



