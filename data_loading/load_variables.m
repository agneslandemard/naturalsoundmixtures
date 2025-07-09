function metrics = load_variables(data_type, params)
global metrics_path 
if nargin < 2
    params = struct();
end
if ~ isfield(params, 'filter_type')
    params.filter_type = 'only_test_retest_soft';
end
if ~ isfield(params, 'threshold')
    params.threshold = 0.3;
end
if ~ isfield(params, 'remove_nans')
    params.remove_nans = true;
end
if ~ isfield(params, 'suffix')
    params.suffix = '';
end
if ~ isfield(params, 'animals')
    params.animals = {'B', 'R', 'L'};
end
switch data_type
    case 'predicted_by_cm'
       
        file_name = ['voxelwise_metrics_from_prediction' params.suffix '.mat'];
        % load prediction accuracy
        load([metrics_path 'full_prediction' params.suffix  '.mat'],'all_data_pred','all_data_true','vxs_to_keep');
        metrics.pred_accuracy = nan(length(vxs_to_keep),1);
        metrics.pred_accuracy(logical(vxs_to_keep)) = nanfastcorr(all_data_pred,all_data_true)';

    case 'real'
        file_name = ['voxelwise_metrics' params.suffix '.mat'];

        % load permutation test results
        test = load([metrics_path 'voxelwise_metrics_permtest.mat'],...
            'corr_mix_solo_nc_all','corr_mix_solo_all','test_retest_all');
        metrics.perm_invariance_nc = test.corr_mix_solo_nc_all;
end


% load invariance
load([metrics_path file_name],...
    'corr_mix_solo_nc_all','corr_mix_solo_all','test_retest_all')
metrics.invariance = corr_mix_solo_all;
metrics.invariance_nc = corr_mix_solo_nc_all;
metrics.test_retest = test_retest_all;

% load animal idx
load(mkpdir([metrics_path 'voxelwise_metrics.mat']), 'an_idx', 'av_resp_all');
metrics.av_resp = av_resp_all;
metrics.an_idx = an_idx;

% get ROI info 
rois = {'MEG','dPEG','VP'};
ri = load([metrics_path 'roi_mapping.mat'],'roi_idx','list_rois');
metrics.roi_idx = nan(size(ri.roi_idx));
for r = 1 : length(rois)
    metrics.roi_idx(ri.roi_idx == find(strcmp(ri.list_rois,rois{r}))) = r;
end

% get slice index
metrics.slice_idx = [];
ct = 0;
for a = 1: length(params.animals)
    anat_info = load([metrics_path  '/anat_params_' params.animals{a} '.mat']);
    for hemi = 1:length(anat_info.hemis)
        slice_ref = anat_info.(anat_info.hemis{hemi}).SliceRef;
        metrics.slice_idx = [metrics.slice_idx; slice_ref + ct];
        ct = ct + max(slice_ref);
    end
end

% get tuning
if exist([metrics_path 'voxel_tuning' params.suffix '.mat'], 'file')
    load([metrics_path 'voxel_tuning' params.suffix '.mat'],'voxels_tuning','rel_tuning');
else
    load([metrics_path 'voxel_tuning.mat'],'voxels_tuning','rel_tuning');
end
metrics.tuning = voxels_tuning;

% only keep reliable tuning (consistent across cv folds)
for ft = 1 : size(metrics.tuning,2)
   metrics.tuning(rel_tuning(:,ft) < 0.5,ft) = nan;
end

%Filter voxels to keep for analysis
filter = LoadVoxelFilter(params.filter_type, params.threshold);
assert(length(filter) == size(corr_mix_solo_nc_all,1))

metrics = correct_struct(metrics, filter, params.remove_nans);

