function [Betas, roi_idx, rois] = load_model_weights(species)
global metrics_path human_path
switch species
    case 'ferret'
        load([metrics_path 'full_prediction.mat'],...
            'B','sigma','vxs_to_keep');
        rois = {'MEG','dPEG','VP'};
        ri = load([metrics_path 'roi_mapping.mat'],'roi_idx','list_rois');
        roi_idx = nan(size(ri.roi_idx));
        for r = 1 : length(rois)
            roi_idx(ri.roi_idx == find(strcmp(ri.list_rois,rois{r}))) = r;
        end
        roi_idx = roi_idx(logical(vxs_to_keep));

        strt = size(B,1)-630+1;
        n_freqs = (size(B,1)-strt+1)/(7*9);
        B = snm(B(strt:end,:,:),3);
        B = B.*sigma';

        Betas = reshape(B',[size(B,2) 7 9 n_freqs]);

    case 'human'
        load([metrics_path 'prediction_humans.mat'], 'Betas', 'roi_idx', 'rois')
end
end