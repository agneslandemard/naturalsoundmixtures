%% Code to re-compute metrics. Can use only subset of sounds.
sounds_subset = 'vocalisations';

an_idx = [];
test_retest_all = [];
corr_mix_solo_all = [];
av_resp_all = [];
corr_mix_solo_nc_all = [];
corr_sum_solo_nc_all = [];
for a = 1:length(animals)
    animal = animals{a};

    keep_data_mat = 0;
    LoadData
    
    metrics = compute_metrics(patterns, fg_id_n, sounds_subset);
    
    %  build final files for all animals together
    test_retest_all = cat(1,  test_retest_all, metrics.test_retest);
    corr_mix_solo_all = cat(1,  corr_mix_solo_all, metrics.corr_mix_solo);
    corr_sum_solo_nc_all = cat(1,  corr_sum_solo_nc_all, metrics.corr_sum_solo_nc);
    corr_mix_solo_nc_all = cat(1,  corr_mix_solo_nc_all, metrics.corr_mix_solo_nc);
    av_resp_all = cat(1, av_resp_all, metrics.av_resp);
    an_idx = cat(1, an_idx, a*ones(size(metrics.av_resp,1),1));
    
end

save(mkpdir([metrics_path 'voxelwise_metrics_' sounds_subset '.mat']),...
    'test_retest_all','av_resp_all','an_idx','corr_mix_solo_all','corr_mix_solo_nc_all',...
    'corr_sum_solo_nc_all')