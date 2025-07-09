function  filter = LoadVoxelFilter(filter_type,threshold)
global metrics_path
% returns boolean array indicating voxels to keep for analyses
% Currently, two criteria are applied:
% - av resp to mixtures should be > 5%CBV
% - nc corr mix solo should exist for at least fg or bg (could be more
% strict and demand that both exist) - itself determined by test-retest
% being positive for both mix and solo

if nargin < 1
    filter_type = 'only_test_retest_soft';
end
if nargin < 2
    threshold = 0.3;
end
disp(threshold)
ref = load(mkpdir([metrics_path 'voxelwise_metrics.mat']),...
    'av_resp_all','corr_mix_solo_nc_all','an_idx','test_retest_all');
resp_threshold = 0.03;

switch filter_type
    case 'test_retest_mix'
        filter = logical((snm(ref.av_resp_all(:,[3 4],:),[2 3]) > resp_threshold ).*(snm_z(ref.test_retest_all(:,[3:4]),[2]) > threshold));
    case 'test_retest_bycd'
        filter = logical((snm(ref.av_resp_all(:,[3 4],:),[2 3]) > resp_threshold ).*(ref.test_retest_all(:,[3:4]) > threshold));       
    case 'onlyresp'
        filter = logical((snm(ref.av_resp_all(:,[3 4],:),[2 3]) > resp_threshold ));
    case 'only_test_retest'
        filter = logical((snm_z(ref.test_retest_all(:,[3 4]),[2]) > threshold));
    case 'only_test_retest_soft'
        filter = logical(any(ref.test_retest_all > threshold,2));
    case 'only_test_retest_hard'
        filter = logical(all(ref.test_retest_all > threshold,2));
    case 'forboth'
        filter = logical((snm(ref.av_resp_all(:,[3 4],:),[2 3]) > resp_threshold ).*all(~isnan(ref.corr_mix_solo_nc_all),2));
    case 'none'
        filter = true(size(ref.av_resp_all,1), 1);
    otherwise
        error('Filter not recognized')
end


end