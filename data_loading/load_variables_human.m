function metrics = load_variables_human(data_type)
global metrics_path

metrics = struct();
switch data_type
    case 'real'
        suff = '';
    case 'predicted_by_cm'
        suff = '_from_prediction';

        load([metrics_path 'prediction_accuracy_humans.mat'],'pred_accuracy');
        metrics.pred_accuracy = pred_accuracy; 
end

load([metrics_path 'voxelwise_metrics' suff '_humans.mat'], ...
    'corr_mix_solo_nc_all','corr_mix_solo_all','test_retest_all', ...
    'an_idx', 'av_resp_all');
metrics.invariance_nc = corr_mix_solo_nc_all;
metrics.invariance = corr_mix_solo_all;
metrics.test_retest = test_retest_all;
metrics.an_idx = an_idx;
metrics.av_resp = av_resp_all; 

load([metrics_path 'roi_mapping_humans.mat'],'roi_idx');
metrics.roi_idx =  roi_idx;

load([metrics_path 'voxel_tuning_humans.mat'],'voxels_tuning');
metrics.tuning = voxels_tuning; 


end

