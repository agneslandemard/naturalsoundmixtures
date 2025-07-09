%% Compute voxelwise metrics

function metrics = compute_metrics(patterns, fg_id_n, sounds_subset)

n_vxs = size(patterns,1);

if ~ exist('blocks_to_use','var')
    blocks_to_use = 2:11; % blocs to take into account
end

n_snips = length(blocks_to_use) * 6 * 2;
respwin = [16:24; 4:12];
LoadStimsInfo

if exist('sounds_subset', 'var')
    switch sounds_subset
        case 'speech'
            snips = stims.fgs_category(mat2vec(fg_id_n(:,blocks_to_use,:))) == 2;
        case 'vocalisations'
            snips = stims.fgs_category(mat2vec(fg_id_n(:,blocks_to_use,:))) == 1;
        case 'nospeech'
            snips = stims.fgs_category(mat2vec(fg_id_n(:,blocks_to_use,:))) ~= 2;
        case 'music'
            snips = stims.fgs_category(mat2vec(fg_id_n(:,blocks_to_use,:))) == 3;
        case 'others'
            snips = stims.fgs_category(mat2vec(fg_id_n(:,blocks_to_use,:))) > 3;
    end
else
    snips  = 1:n_snips;
end

test_retest = nan(n_vxs,4);
corr_mix_solo = nan(n_vxs,2,2);
av_resp = nan(n_vxs,4,2);
corr_sum_solo_nc = nan(n_vxs,2);
corr_mix_solo_nc = nan(n_vxs,2);


fgs = reshape(snm(patterns(:,ismember(resp_wins,respwin(1,:)),:,blocks_to_use,[1 1],:),2),n_vxs,n_snips,2);
bgs = reshape(snm(patterns(:,ismember(resp_wins,respwin(2,:)),:,blocks_to_use,[4 5],:),2),n_vxs,n_snips,2);
mixs_fgs = reshape(snm(patterns(:,ismember(resp_wins,respwin(1,:)),:,blocks_to_use,[2 3],:),2),n_vxs,n_snips,2);
mixs_bgs = reshape(snm(patterns(:,ismember(resp_wins,respwin(2,:)),:,blocks_to_use,[2 3],:),2),n_vxs,n_snips,2);


% Average response
av_resp(:,1,:) =  snm(fgs(:,snips,:),2);
av_resp(:,2,:) =  snm(bgs(:,snips,:),2);
av_resp(:,3,:) =  snm(mixs_fgs(:,snips,:),2);
av_resp(:,4,:) =  snm(mixs_bgs(:,snips,:),2);

% Test retest
test_retest(:,1) = nancorr_pairwise(fgs(:,snips,1)',fgs(:,snips,2)');
test_retest(:,2) = nancorr_pairwise(bgs(:,snips,1)',bgs(:,snips,2)');
test_retest(:,3) = nancorr_pairwise(mixs_fgs(:,snips,1)',mixs_fgs(:,snips,2)');
test_retest(:,4) = nancorr_pairwise(mixs_bgs(:,snips,1)',mixs_bgs(:,snips,2)');

% Invariance (non noise-corrected)
for rep = 1:2
    corr_mix_solo(:,1,rep) = nancorr_pairwise(mixs_fgs(:,snips,rep)',fgs(:,snips,setdiff(1:2,rep))');
    corr_mix_solo(:,2,rep) = nancorr_pairwise(mixs_bgs(:,snips,rep)',bgs(:,snips,setdiff(1:2,rep))');
end


% Invariance (noise-corrected)
if all(round(test_retest(:,:),4) == 1)
    % if test-retest if perfect, no point in doing
    % noise-correction. Happens with fake data (predicted..)
    corr_mix_solo_nc(:,:)  = corr_mix_solo(:,:,1,pm);
    corr_sum_solo_nc(:,:) = nan;
else

    % noise corrected corr mix solo
    vxs_for_nc = snm(test_retest(:,1),2) > 0 & snm(test_retest(:,3),2) > 0;
    % this condition should ensure that obtained values will be real

    corr_mix_solo_nc(vxs_for_nc,1) = noise_corrected_correlation_pairwise(mixs_fgs(vxs_for_nc,snips,1)',mixs_fgs(vxs_for_nc,snips,2)',fgs(vxs_for_nc,snips,1)',fgs(vxs_for_nc,snips,2)');
    % what would happen if just summing both contributions
    corr_sum_solo_nc(vxs_for_nc,1) = noise_corrected_correlation_pairwise(fgs(vxs_for_nc,snips,1)'+bgs(vxs_for_nc,snips,1)',fgs(vxs_for_nc,snips,2)'+bgs(vxs_for_nc,snips,2)',fgs(vxs_for_nc,snips,1)',fgs(vxs_for_nc,snips,2)');

    vxs_for_nc = snm(test_retest(:,2),2) > 0 & snm(test_retest(:,4),2) > 0;
    corr_mix_solo_nc(vxs_for_nc,2) = noise_corrected_correlation_pairwise(mixs_bgs(vxs_for_nc,snips,1)',mixs_bgs(vxs_for_nc,snips,2)',bgs(vxs_for_nc,snips,1)',bgs(vxs_for_nc,snips,2)');
    corr_sum_solo_nc(vxs_for_nc,2) = noise_corrected_correlation_pairwise(fgs(vxs_for_nc,snips,1)'+bgs(vxs_for_nc,snips,1)',fgs(vxs_for_nc,snips,2)'+bgs(vxs_for_nc,snips,2)',bgs(vxs_for_nc,snips,1)',bgs(vxs_for_nc,snips,2)');
end

metrics = struct('test_retest', test_retest, 'av_resp',av_resp,...
    'corr_mix_solo', corr_mix_solo_nc, 'corr_mix_solo_nc', corr_mix_solo_nc,...
    'corr_sum_solo_nc', corr_sum_solo_nc);

end

