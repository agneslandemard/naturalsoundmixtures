%% Load and format fUS data for all animals

all_data_true = [];
for a = 1:length(animals)
    animal = animals{a};
    keep_data_mat = 1;
    LoadData
    clear patterns
    blocs = stims.blocs_ids(1:12:end,1:2);
    non_baseline = sum(blocs,2) ~= 0;
    [n_vxs_or, n_tps_or, n_blcs_or, n_reps_or] = size(norm_data_mat);

    all_data_true = cat(1, all_data_true, norm_data_mat(:,:,:,:));
    clear norm_data_mat
end

% keep only responsive voxels, could use other filters!
params.filter_type = 'onlyresp'; 
params.threshold = 0.25;
filter = LoadVoxelFilter(params.filter_type,params.threshold); 

all_data_true = all_data_true(filter,:,:,:);
resp_wins = -5:40;
patterns =  FormatIntoPatterns(all_data_true,stims,param,resp_wins);
patterns_r = patterns - nanmean(patterns,2:6); % recenter across sounds
clear patterns all_data_true

load(mkpdir([metrics_path 'voxelwise_metrics.mat']), 'an_idx');
an_idx = an_idx(filter);

%% Compute cross-correlation matrices

corr_across = 'sounds'; % can be 'sounds' or 'voxels'

n_vxs = size(patterns_r,1);
bl = 2:11; % blocs to take into account. removing first and last by default (edge effects)
n_snips = length(bl) * 6 * 2;
n_tps = size(patterns_r,2);
win_size = 0; % smoothing if >0

fgs = nan(n_vxs,n_tps,n_snips,2);
bgs = nan(n_vxs,n_tps,n_snips,2);
mixs_fgs = nan(n_vxs,n_tps,n_snips,2);
mixs_bgs = nan(n_vxs,n_tps,n_snips,2);

for t = 1:n_tps
    intw = intersect(t-win_size:t+win_size,1:n_tps);
    fgs(:,t,:,:) = reshape(snm(patterns_r(:,intw,:,bl,[1 1],:),2),n_vxs,1,n_snips,2);
    bgs(:,t,:,:)  = reshape(snm(patterns_r(:,intw,:,bl,[4 5],:),2),n_vxs,1,n_snips,2);
    mixs_fgs(:,t,:,:)  = reshape(snm(patterns_r(:,intw,:,bl,[2 3],:),2),n_vxs,1,n_snips,2);
    mixs_bgs(:,t,:,:)  = reshape(snm(patterns_r(:,intw,:,bl,[2 3],:),2),n_vxs,1,n_snips,2);
end

% run correlations
switch corr_across
    case 'voxels'
        n = n_snips;
    case 'sounds'
        n = n_vxs;
end
crossc = nan(n_tps,n_tps,2,n);
sf_crossc = nan(n_tps,n_tps,2,n);
nc_crossc = nan(n_tps,n_tps,2,n);
crosstrt = nan(n_tps,n_tps,2,n);
sf_crosstrt = nan(n_tps,n_tps,2,n);
for k = 1:n
    switch corr_across
        case 'voxels'

            crossc(:,:,1,k) = corr(snm(fgs(:,:,k,:),4),snm(mixs_fgs(:,:,k,:),4),'rows','pairwise');
            crossc(:,:,2,k) = corr(snm(bgs(:,:,k,:),4),snm(mixs_bgs(:,:,k,:),4),'rows','pairwise');

            crosstrt(:,:,1,k) = corr(snm(fgs(:,:,k,1),4),snm(fgs(:,:,k,2),4),'rows','pairwise');
            crosstrt(:,:,2,k) = corr(snm(bgs(:,:,k,1),4),snm(bgs(:,:,k,2),4),'rows','pairwise');
            crosstrt(:,:,3,k) = corr(snm(mixs_bgs(:,:,k,1),4),snm(mixs_bgs(:,:,k,2),4),'rows','pairwise');

            shuf = randperm(n_snips);
            sf_crossc(:,:,1,k) = corr(snm(fgs(:,:,shuf(k),:),4),snm(mixs_fgs(:,:,k,:),4),'rows','pairwise');
            sf_crossc(:,:,2,k) = corr(snm(bgs(:,:,shuf(k),:),4),snm(mixs_bgs(:,:,k,:),4),'rows','pairwise');

            sf_crosstrt(:,:,1,k) = corr(snm(fgs(:,:,shuf(k),1),4),snm(fgs(:,:,k,2),4),'rows','pairwise');
            sf_crosstrt(:,:,2,k) = corr(snm(bgs(:,:,shuf(k),1),4),snm(bgs(:,:,k,2),4),'rows','pairwise');
            sf_crosstrt(:,:,3,k) = corr(snm(mixs_bgs(:,:,shuf(k),1),4),snm(mixs_bgs(:,:,k,2),4),'rows','pairwise');


        case 'sounds'
            crossc(:,:,1,k) = corr(snm(fgs(k,:,:,:),4)',snm(mixs_fgs(k,:,:,:),4)','rows','pairwise');
            crossc(:,:,2,k) = corr(snm(bgs(k,:,:,:),4)',snm(mixs_bgs(k,:,:,:),4)','rows','pairwise');

            crosstrt(:,:,1,k) = corr(snm(fgs(k,:,:,1),4)',snm(fgs(k,:,:,2),4)','rows','pairwise');
            crosstrt(:,:,2,k) = corr(snm(bgs(k,:,:,1),4)',snm(bgs(k,:,:,2),4)','rows','pairwise');
            crosstrt(:,:,3,k) = corr(snm(mixs_bgs(k,:,:,1),4)',snm(mixs_bgs(k,:,:,2),4)','rows','pairwise');

            shuf = randperm(n_snips);
            sf_crosstrt(:,:,1,k) = corr(snm(fgs(k,:,shuf,1),4)',snm(fgs(k,:,:,2),4)','rows','pairwise');
            sf_crosstrt(:,:,2,k) = corr(snm(bgs(k,:,shuf,1),4)',snm(bgs(k,:,:,2),4)','rows','pairwise');
            sf_crosstrt(:,:,3,k) = corr(snm(mixs_bgs(k,:,shuf,1),4)',snm(mixs_bgs(k,:,:,2),4)','rows','pairwise');

            sf_crossc(:,:,1,k) = corr(snm(fgs(k,:,shuf,:),4)',snm(mixs_fgs(k,:,:,:),4)','rows','pairwise');
            sf_crossc(:,:,2,k) = corr(snm(bgs(k,:,shuf,:),4)',snm(mixs_bgs(k,:,:,:),4)','rows','pairwise');

    end
end

%% Test-retest reliability (Figure 1)

cm = cbrewer('seq','OrRd',64);
f1 = figure('Position',[1000 632 1161 706]);
f2 = figure;
widths = nan(3,1);
wins = {3:24, 12+ (3:24), 3:36};
for cd = 1:3

    to_plot_true = snm_z(crosstrt(:,:,cd,:),4);
    to_plot_shuf = snm_z(sf_crosstrt(:,:,cd,:),4);

    figure(f1)
    subplot(2,3,cd)
    imagesc(resp_wins./param.exp.SR,resp_wins./param.exp.SR,to_plot_true)
    hline(stims.events,'k-')
    vline(stims.events,'k-')
    m = prctile(mat2vec(snm_z(crosstrt,4)),99);
    mi = min(mat2vec(snm_z(crosstrt,4)));
    xlabel('Time in solo (s)')
    ylabel('Time in solo (s)')
    title(['Test retest ' groups_names{cd}])
    clim([mi m])
    colormap(cm)
    axis equal tight
    colorbar

    subplot(2,3,cd+3)
    imagesc(resp_wins./param.exp.SR,resp_wins./param.exp.SR,to_plot_shuf)
    hline(stims.events,'k-')
    vline(stims.events,'k-')
    xlabel('Time in solo (s)')
    ylabel('Time in solo (s)')
    title(['Test retest ' groups_names{cd} ' with shuffled sounds'])
    clim([mi m])
    colormap(cm)
    axis equal tight
    colorbar


    figure(f2)
    hold all
    lags = -30:30;
    t_range = find(ismember(resp_wins,wins{cd}));
    diagw = nan(length(t_range),length(lags));
    for t = 1:length(t_range)
        [int,i] = intersect(t_range(t)+lags,1:length(resp_wins));
        diagw(t,i) = to_plot_true(t_range(t),int);
    end
    decay = snm(diagw,1)./snm(diagw(:,lags==0),1);
    decays_neg = decay(lags<=0);
    decays_pos = decay(lags>=0);
    lags_pos = lags(lags>=0)./param.exp.SR;
    sym_decay = (decays_neg(end:-1:1)+decays_pos)/2;
    plot(lags_pos,sym_decay,'LineWidth',2)
    legend(groups_names)
    xlabel('Lag (s)')
    ylabel('Corr test-retest')
    hline(0.8)


end
saveas(f1,[figures_path 'CrossTestRetest_v' version '_corr_across' corr_across '.' ext])
saveas(f2,[figures_path 'TestRetestWidth_v' version '_corr_across' corr_across '.' ext])


%% Invariance (Figure S1)
cm = cbrewer('seq','OrRd',64);

f1 = figure;
f2 = figure;

m = max(mat2vec(snm_z(crossc,4)));
mi = -0.01;
for cd = 1:2

    to_plot_true = snm_z(crossc(:,:,cd,:),4);
    to_plot_shuf = snm_z(sf_crossc(:,:,cd,:),4);

    figure(f1)
    subplot(2,2,cd)
    imagesc(resp_wins./param.exp.SR, resp_wins./param.exp.SR, to_plot_true)
    hline(stims.events,'k-')
    vline(stims.events,'k-')
    xlabel('Time in solo (s)')
    ylabel('Time in mix (s)')
    title(['Corr mix-' groups_names{cd}])
    clim([mi m])
    colormap(cm)
    axis equal tight
    colorbar

    subplot(2,2,cd+2)
    imagesc(resp_wins./param.exp.SR, resp_wins./param.exp.SR, to_plot_shuf)
    hline(stims.events,'k-')
    vline(stims.events,'k-')
    xlabel('Time in solo (s)')
    ylabel('Time in mix (s)')
    title(['Corr mix-' groups_names{cd} ' with shuffled sounds'])
    clim([mi m])
    colormap(cm)
    axis equal tight
    colorbar


    figure(f2)
    hold all
    plot(resp_wins./param.exp.SR, diag(to_plot_true),'color', squeeze(cmaps.cm(1,cd,:))','LineWidth',2)
    plot(resp_wins./param.exp.SR, diag(to_plot_shuf),'color', squeeze(cmaps.cm(2,cd,:))','LineWidth',2)
    vline(stims.events,'k-') 
    xlabel('Time in solo (s)')
    ylabel('Corr mix solo')
    ylim([mi m])

end

saveas(f1,[figures_path 'CrossCorr_v' version '_corr_across' corr_across '.' ext])
saveas(f2,[figures_path 'CrossCorrDiags_v' version '_corr_across' corr_across '.' ext])


