% Contains analyses to compare background or foreground invariance, as well
% as other metrics, across ROIS, for both ferrets and humans.

%% Load metrics for chosen species, including true and predicted invariance

% Choose species
species = 'human';

all_metrics = struct();
save_suffix = '';

switch species
    case 'ferret'
        pred_in_isolation = false;
        

        rois = {'MEG', 'dPEG', 'VP'};
        
        all_metrics.true = load_variables('real');
       % all_metrics.true = load_variables('real', struct('suffix', '_vocalisations'));

        if pred_in_isolation
            all_metrics.pred = load_variables('predicted_by_cm', struct('suffix', '_solos'));
        else
            all_metrics.pred = load_variables('predicted_by_cm');
        end
        
        
        xpos_rois = [1 3 4];
        xpos_rois2 = [1 3 4];
        roi_grps = {1, 2, 3};
        subj_ids = animals;

    case 'human'
        
        rois = {'primary', 'nonprimary'};
        
        all_metrics.true = load_variables_human('real');
        all_metrics.pred = load_variables_human('predicted_by_cm');
        
        xpos_rois = [1 2];
        xpos_rois2 = [1 2];
        roi_grps = {1, 2};
        subj_ids =  arrayfun(@(x)[ 'subj ' num2str(x)], unique(all_metrics.true.an_idx)', 'UniformOutput', false);
end

roi_idx = all_metrics.true.roi_idx;
an_idx = all_metrics.true.an_idx;

[vals, feats] = LoadCMparams;
n_feats = length(feats);
n_rois = length(rois);

a_range = unique(an_idx(~isnan(roi_idx)))';

mean_func = @nanmedian; % function to compute statistic across voxels (usually median)
vxs_keep = all(~isnan(all_metrics.true.invariance_nc),2);
all_metrics.true.av_resp_all = snm(all_metrics.true.av_resp, [2, 3]);
all_metrics.true.test_retest_all = snm(all_metrics.true.test_retest, 2);

%% Invariance across subjects and ROIs - Figure 2F, I - Figure 4D, H - Figure 5C, E, G, I

types = {'true','pred'};
plotWhat = 'invariance_nc';

figure('Position', [94 38 328 958]);
clear ax
n = 3; 

all_vals = nan(length(roi_grps),n,max(an_idx),length(types));
n_vx = zeros(length(roi_grps),max(an_idx),length(types));  
for t = 1:length(types)
   if ~isfield(all_metrics.(types{t}), plotWhat) || isempty(all_metrics.(types{t}).(plotWhat))
       continue
   end
    n = min(3,size(all_metrics.(types{t}).(plotWhat),2));
    % show data for each animal
    for a = a_range
        voxel_subset = (an_idx == a);
       
        for r = 1: length(roi_grps)
            vxs_roi = find(logical(ismember(roi_idx,roi_grps{r}).*voxel_subset.*vxs_keep));
            n_vx(r,a,t) = length(vxs_roi);
            
            for g = 1:n
                all_vals(r,g,a,t) = squeeze(mean_func(snm(all_metrics.(types{t}).(plotWhat)(vxs_roi,g,:),3)));
            end
            
        end

        if n==2
            all_vals(:,3,:,t) = all_vals(:,1,:,t) - all_vals(:,2,:,t);
        end

        for g = 1:n
            subplot(length(types),n,g+(t-1)*n) 
            hold all
            sz = (n_vx(:,a,t)' - 100)/(1700 - 100);
            scatter(xpos_rois2,all_vals(:,g,a,t),sz*50,0.8*[1 1 1],'filled')
            if strcmp(species,'ferret') && strcmp(subj_ids{a}, 'L')
                plot(xpos_rois2,all_vals(:,g,a,t),'color',0.7*[1 1 1],'LineWidth',1)
            else
                plot(xpos_rois2,all_vals(:,g,a,t),'color',0.7*[1 1 1])
            end
        end

    end

    for r = 1: n_rois
        sz = (sum(n_vx(r,:,t),2)' - 100)/(1700 - 100);
        for g = 1:n

            ax(g,t) = subplot(length(types),n,g+(t-1)*n);
            hold all
            vxs_roi = logical((roi_idx == r).*vxs_keep);
            to_plot = squeeze(mean_func(all_metrics.(types{t}).(plotWhat)(vxs_roi,g,:),[1, 3]))';
            scatter(xpos_rois(r),to_plot, sz*50 ...
                ,cmaps.roi_colors(r,:),'o','LineWidth',2);
            title([types{t} ' ' groups_names{g}])
            xticks(xpos_rois)
            xticklabels(rois)
            xlim([0 max(xpos_rois)+1])
             
            if strcmp(plotWhat, 'av_resp') || strcmp(plotWhat, 'av_resp_all')
                 ylim([0 0.2])
          
            else
                ylim([-0.2 1])
            end
            
            ylabel([regexprep(plotWhat, '_', ' ') groups_names{g}])
        end
    end

end

if ~ isempty(ext)
    saveas(gcf,[figures_path plotWhat  '_' species save_suffix '_by_ROI.' ext])
end



%% True vs predicted for bins regardless of ROI - Figure 4C, G - Figure 5D, H

plotWhat = 'invariance_nc';

bins = -0.2:0.1:1;
figure('Position',[1000 729 1.1363e+03 528],'Renderer','Painters');
for g = 1 : 2
   subplot(1,2,g)
   hold all
   if g==1
       P.bins1 = 0:0.02:1;
       P.bins2 = 0:0.02:1;
   else
       P.bins1 = -0.23:0.02:0.8;
       P.bins2 = -0.23:0.02:0.8;
   end

   P.n_bins=length(P.bins1);

   P.r = [0.5,0.5,0.5];
   P.b = [1, 1, 1];
   heatScatter(all_metrics.true.(plotWhat)(vxs_keep,g),...
       all_metrics.pred.(plotWhat)(vxs_keep,g),[], P);


   for a = 1:length(a_range)
       voxel_subset = logical((an_idx == a_range(a)).*vxs_keep);

       [~,~,loc] = histcounts(all_metrics.true.(plotWhat)(voxel_subset,g),bins);
       val_bins = nan(length(bins),1);
       tmp = all_metrics.pred.(plotWhat)(voxel_subset,g);
       for b = 1:length(bins)
           if sum(loc==b)>50
               val_bins(b) = nanmedian(tmp(loc==b));
           end
       end
       plot(bins +0.05, val_bins,'Color',squeeze(cmaps.cm(1,g,:)),'LineWidth',1)

       xlabel('Measured')
       ylabel('Predicted')
       
   end
end

sgtitle(regexprep(plotWhat, '_', ' '))
if ~ isempty(ext)
    saveas(gcf,[figures_path plotWhat '_' species save_suffix '_compare_true_predicted_bins.' ext])
end


%% Invariance vs prediction accuracy for each ROI - Figure S5D, E, H, I

plotWhat = 'invariance_nc';

bins = -0.2:0.1:1;
f1 = figure('Position',[1000 729 1.1363e+03 528],'Renderer','Painters');
f2 = figure('Position', [1000 914.3333 827 423.6667]);
for g = 1 : 2
    for  r = 1:n_rois
        
        if g==1
            P.bins1 = 0:0.02:1;
            P.bins2 = 0:0.02:1;
        else
            P.bins1 = -0.23:0.02:0.8;
            P.bins2 = -0.23:0.02:0.8;
        end
        P.n_bins=length(P.bins1);

        P.r = [0.5,0.5,0.5];
        P.b = [1, 1, 1];
        vxs_roi = vxs_keep & (roi_idx == r);
        figure(f1)
        subplot(2,n_rois,(g-1)*n_rois +r)
        hold all
        heatScatter(all_metrics.pred.pred_accuracy(vxs_roi),all_metrics.pred.(plotWhat)(vxs_roi,g),[], P);

        for a = 1:length(a_range)
            voxel_subset = logical((an_idx == a_range(a)).*vxs_roi);

            [~,~,loc] = histcounts(all_metrics.pred.pred_accuracy(voxel_subset),bins);
            val_bins = nan(length(bins),1);
            tmp = all_metrics.pred.(plotWhat)(voxel_subset,g);
            for b = 1:length(bins)
                if sum(loc==b)>10
                    val_bins(b) = nanmedian(tmp(loc==b));
                end
            end
            figure(f1)
            subplot(2,n_rois,(g-1)*n_rois +r)
            hold all
            plot(bins +0.05, val_bins,'Color',squeeze(cmaps.cm(1,g,:)),'LineWidth',1)
            title([rois{r} ' ' groups_names{g}])

            xlabel('Prediction accuracy')
            ylabel(['Predicted ' regexprep(plotWhat, '_', ' ') ' ' groups_names{g}])

            figure(f2)
            subplot(1,2,g)
            hold all
            plot(bins +0.05, val_bins,'Color',cmaps.roi_colors(r,:),'LineWidth',0.5)
            xlabel('Prediction accuracy')
            ylabel(['Predicted ' regexprep(plotWhat, '_', ' ') ' ' groups_names{g}])
            linkaxes
        end

        voxel_subset = logical(vxs_roi);

        [~, ~, loc] = histcounts(all_metrics.pred.pred_accuracy(voxel_subset),bins);
        val_bins = nan(length(bins),1);
        tmp = all_metrics.pred.(plotWhat)(voxel_subset,g);
        for b = 1:length(bins)
            if sum(loc==b)>50
                val_bins(b) = nanmedian(tmp(loc==b));
            end
        end
        figure(f2)
        subplot(1,2,g)
        hold all
        plot(bins +0.05, val_bins,'Color',cmaps.roi_colors(r,:),'LineWidth',2)
        xlabel('Prediction accuracy')

    end
end

sgtitle(regexprep(plotWhat, '_', ' '))
if ~ isempty(ext)
    saveas(f1,[figures_path plotWhat '_' species save_suffix '_compare_vs_pred_accuracy_by_roi_bins.' ext])
    saveas(f2,[figures_path plotWhat '_' species save_suffix '_compare_vs_pred_accuracy_by_roi.' ext])
end

%% True vs predicted for bins for each ROI

plotWhat = 'invariance_nc';

bins = -0.2:0.1:1;
f1 = figure('Position', [1000 729 1.1363e+03 528], 'Renderer', 'Painters');
f2 = figure();
for g = 1 : 2
    for  r = 1:n_rois
        
        if g==1
            P.bins1 = 0:0.02:1;
            P.bins2 = 0:0.02:1;
        else
            P.bins1 = -0.23:0.02:0.8;
            P.bins2 = -0.23:0.02:0.8;
        end

        P.n_bins=length(P.bins1);

        P.r = [0.5, 0.5, 0.5];
        P.b = [1, 1, 1];
        vxs_roi = vxs_keep & (roi_idx == r);
        figure(f1)
        subplot(2,n_rois,(g-1)*n_rois +r)
        hold all
        heatScatter(all_metrics.true.(plotWhat)(vxs_roi,g),...
            all_metrics.pred.(plotWhat)(vxs_roi,g),[], P);


        for a = 1:length(a_range)
            voxel_subset = logical((an_idx == a_range(a)).*vxs_roi);

            [~,~,loc] = histcounts(all_metrics.true.(plotWhat)(voxel_subset,g),bins);
            val_bins = nan(length(bins),1);
            tmp = all_metrics.pred.(plotWhat)(voxel_subset,g);
            for b = 1:length(bins)
                if sum(loc==b)>10
                    val_bins(b) = nanmedian(tmp(loc==b));
                end
            end
            figure(f1)
            subplot(2,n_rois,(g-1)*n_rois +r)
            hold all
            plot(bins +0.05, val_bins,'Color',squeeze(cmaps.cm(1,g,:)),'LineWidth',1)
            title([rois{r} ' ' groups_names{g}])

            xlabel('Measured')
            ylabel('Predicted')

            figure(f2)
            subplot(1,2,g)
            hold all
            plot(bins +0.05, val_bins,'Color',cmaps.roi_colors(r,:),'LineWidth',1)
            xlabel('Measured')
            ylabel('Predicted')

        end
    end
end

sgtitle(regexprep(plotWhat, '_', ' '))
if ~ isempty(ext)
    saveas(f1,[figures_path plotWhat '_' species save_suffix '_compare_true_predicted_by_roi_bins.' ext])
    saveas(f2,[figures_path plotWhat '_' species save_suffix '_compare_true_predicted_by_roi.' ext])
end

%% Invariance vs test-retest for bins for each ROI

plotWhat = 'invariance_nc';
as_function_of = 'test_retest';

bins = -0.2:0.1:1;
f1 = figure('Position', [1000 729 1.1363e+03 528], 'Renderer', 'Painters');
f2 = figure();
for g = 1 : 2
    for  r = 1:n_rois
        
  
        P.bins1 = -0.2:0.02:1;
        P.bins2 = -0.2:0.02:1;
        P.n_bins=length(P.bins1);
        P.r = [0.5, 0.5, 0.5];
        P.b = [1, 1, 1];
        vxs_roi = vxs_keep & (roi_idx == r);
        figure(f1)
        subplot(2,n_rois,(g-1)*n_rois +r)
        hold all
        heatScatter(snm(all_metrics.true.(as_function_of)(vxs_roi,:),2),...
            all_metrics.true.(plotWhat)(vxs_roi,g),[], P);


        for a = 1:length(a_range)
            voxel_subset = logical((an_idx == a_range(a)).*vxs_roi);

            [~,~,loc] = histcounts(snm(all_metrics.true.(as_function_of)(voxel_subset,:),2),bins);
            val_bins = nan(length(bins),1);
            tmp = all_metrics.true.(plotWhat)(voxel_subset,g);
            for b = 1:length(bins)
                if sum(loc==b)>10
                    val_bins(b) = nanmedian(tmp(loc==b));
                end
            end
            figure(f1)
             subplot(2,n_rois,(g-1)*n_rois +r);
            hold all
            plot(bins +0.05, val_bins,'Color',squeeze(cmaps.cm(1,g,:)),'LineWidth',1)
            title([rois{r} ' ' groups_names{g}])

            xlabel(as_function_of)
            ylabel(plotWhat)

            figure(f2)
            ax(g) =subplot(1,2,g);
            hold all
            plot(bins +0.05, val_bins,'Color',cmaps.roi_colors(r,:),'LineWidth',0.5)
            xlabel(as_function_of)
            ylabel(plotWhat)
            ylim([-0.2, 1])

        end


        figure(f2)
        voxel_subset = logical(vxs_roi); 
        [~,~,loc] = histcounts(snm(all_metrics.true.(as_function_of)(voxel_subset,:),2),bins);
        val_bins = nan(length(bins),1);
        tmp = all_metrics.true.(plotWhat)(voxel_subset,g);
        for b = 1:length(bins)
            if sum(loc==b)>10
                val_bins(b) = nanmedian(tmp(loc==b));
            end
        end
      
        plot(bins +0.05, val_bins,'Color',cmaps.roi_colors(r,:),'LineWidth',2)
        linkaxes(ax)
         
        axis equal
        ylim([-0.2, 1.2])
    end
end

sgtitle(regexprep(plotWhat, '_', ' '))
if ~ isempty(ext)
    saveas(f1,[figures_path plotWhat '_' species save_suffix '_compare_' plotWhat '_' as_function_of '_by_roi_bins.' ext])
    saveas(f2,[figures_path plotWhat '_' species save_suffix '_compare_' plotWhat '_' as_function_of '_by_roi.' ext])
end

%% STATS: background invariance vs foreground invariance

if isfield(all_metrics.true, 'perm_invariance_nc')
    pvals = nan(n_rois,1);
    diff_fg_bg_true = all_metrics.true.invariance_nc(:, 1, :) - all_metrics.true.invariance_nc(:,2, :);
    diff_fg_bg_shuffled = all_metrics.true.perm_invariance_nc(:, 1, :) - all_metrics.true.perm_invariance_nc(:,2, :);

    for r = 1:n_rois
        vxs_roi = find(logical(ismember(roi_idx,r).*vxs_keep));
        n_total = size(diff_fg_bg_shuffled,3) + 1;
        pvals(r) = (1+ sum(abs(nanmean(diff_fg_bg_shuffled(vxs_roi,:,:),1)) >=  abs(nanmean(diff_fg_bg_true(vxs_roi,:,:),1))))/n_total;
        disp([rois{r} ': p=' num2str(round(pvals(r),4))])
    end
end

%% STATS: invariance across ROIs (pairwise) by shuffling ROI labels within each subject

r1 = 1; % region 1
r2 = [2]; % region 2 (can be multiple)
n_perms = 1000;

types = {'true','pred'};
metric_tested = 'invariance_nc';


for t = 1:length(types)

    dtype = types{t};
    metrics = all_metrics.(dtype);   
    figure('Position', [1.1103e+03 191.6667 1.1927e+03 1086]);

    shuffled_invariance_by_animal = nan(n_perms, length(a_range), 2);
    invariance_by_animal = nan(length(a_range), 2,2);
    diff_invariance_by_animal = nan(length(a_range),2);

    vxs_all = logical(ismember(roi_idx , [r1, r2]).*vxs_keep);%
    animal_invariance_all = metrics.(metric_tested)(vxs_all,:);
    an_idx_all = an_idx(vxs_all);
    true_roi_all = roi_idx(vxs_all);
    fake_roi_all = nan(length(true_roi_all), n_perms);

    invariance_total  = mean_func(animal_invariance_all(ismember(true_roi_all, r2),:),1) - mean_func(animal_invariance_all(true_roi_all == r1,:),1);
    
    for a = 1:length(a_range)

        vxs_an = logical(ismember(roi_idx , [r1, r2]).*vxs_keep.*(an_idx==a_range(a)));%
        animal_invariance = metrics.(metric_tested)(vxs_an,:);
        true_roi = roi_idx(vxs_an);

        val_r2 = mean_func(animal_invariance(ismember(true_roi, r2),:),1) ;
        val_r1 = mean_func(animal_invariance(true_roi == r1,:),1);
        invariance_by_animal(a,: ,1) = val_r2; 
        invariance_by_animal(a,: ,2) = val_r1; 
        diff_invariance_by_animal(a,:) =  invariance_by_animal(a, :,1)- invariance_by_animal(a, :,2);
        for sf = 1:n_perms
            fake_roi = true_roi(randperm(length(true_roi)));
            shuffled_invariance_by_animal(sf, a, :) = mean_func(animal_invariance(ismember(fake_roi, r2),:),1) - mean_func(animal_invariance(fake_roi == r1,:),1);
            fake_roi_all(an_idx_all == a_range(a), sf) = fake_roi;
        end

        shuffled_invariance_by_animal(end, a, :)  = diff_invariance_by_animal(a,:);

        for cd = 1:2
            subplot(length(a_range)+1,2,cd+(a-1)*2)
            histogram(shuffled_invariance_by_animal(:,a,cd))

            pval = sum(abs(shuffled_invariance_by_animal(:, a, cd)) >= abs(diff_invariance_by_animal(a, cd)))/(n_perms+1);
            disp([dtype ' ' groups_names{cd} ' ' subj_ids{a} ])
            disp([rois{r2} ': ' num2str(invariance_by_animal(a, cd,1),2) ', ' rois{r1} ': ' num2str(invariance_by_animal(a, cd,2),2) ', p=' num2str(pval)]);
          
            vline(diff_invariance_by_animal(a, cd))
            text(diff_invariance_by_animal(a, cd),20,['Md= ' num2str(diff_invariance_by_animal(a, cd),3)])
            title([subj_ids{a} ', ' groups_names{cd} ', p = ' num2str(pval,3)])
        end

    end

    shuffled_invariance_total = nanmean(shuffled_invariance_by_animal,2);
    invariance_total = nanmean(diff_invariance_by_animal,1);

    for cd = 1:2
        subplot(length(a_range)+1,2,cd+(a)*2)
        histogram(shuffled_invariance_total(:,cd))

        pval = sum(abs(shuffled_invariance_total(:, cd)) >= abs(invariance_total(cd)))/(n_perms+1);
        disp([dtype ' ' groups_names{cd} ' all' ])
        disp([rois{r2} ': ' num2str(nanmean(invariance_by_animal(:, cd,1)),2) ', ' rois{r1} ': ' num2str(nanmean(invariance_by_animal(:, cd,2)),2) ', p=' num2str(pval)]);
          

        vline(invariance_total(cd))
        text(invariance_total(cd),20,['Md= ' num2str(invariance_total(cd),3)])
        title(['All, ' groups_names{cd} ', p = ' num2str(pval,3)])
    end

    
    sgtitle(['ROIs: ' rois{r2} ' - ' rois{r1} ' for ' dtype])
    saveas(gcf,[figures_path  'permutation_test_rois_by_animal_shuffle_animal' species save_suffix '_' metric_tested '_' num2str(r1) '_' num2str(r2) '_' dtype '.png'])

end



%% Invariance by tuning to rates, for each region - Figure 3H

type = 'true';
ft = 3; % to look at rates
lim = 5; % limit between low and high rates

tuning = all_metrics.true.tuning;
an_idx = all_metrics.true.an_idx;

figure('Position', [94 38 899 554.3333]);
clear ax
xpos = @(r) 1+ (r-1)*n_rois;
all_vals = nan(n_rois, 2, length(subj_ids));
for g = 1:2
    subplot(1,2,g)
    hold all
    for r = 1: n_rois
        
        % show data for each animal
        for a = a_range
            vxs_roi = (roi_idx ==r).*vxs_keep.*(an_idx == a).*(tuning(:,ft) < lim );
            if ~ any(vxs_roi)
                continue
            end
            to_plot1 = squeeze(mean_func(all_metrics.(type).invariance_nc(logical(vxs_roi),g),1))';
            scatter(xpos(r),to_plot1,length(find(vxs_roi))*0.1,0.8*[1 1 1],'filled')

            vxs_roi = (roi_idx ==r).*vxs_keep.*(an_idx == a).*(tuning(:,ft) > lim);
            if any(vxs_roi)
                to_plot2 = squeeze(mean_func(all_metrics.(type).invariance_nc(logical(vxs_roi),g),1))';
                scatter(xpos(r)+1,to_plot2,length(find(vxs_roi))*0.1,0.8*[1 1 1],'filled')
                plot(xpos(r) + (0:1),[to_plot1, to_plot2],'color',0.7*[1 1 1])
    
                all_vals(r, 1, a) = to_plot1;
                all_vals(r, 2, a) = to_plot2;
            end
        end
        % all subj_ids together
     
        vxs_roi = logical((roi_idx == r).*vxs_keep.*(tuning(:,ft) < lim));
        to_plot = squeeze(mean_func(all_metrics.(type).invariance_nc(vxs_roi,g),1))';
        scatter(xpos(r),to_plot,length(find(vxs_roi))*0.1 ...
            ,cmaps.roi_colors(r,:),'o','LineWidth',2);
         
        vxs_roi = logical((roi_idx == r).*vxs_keep.*(tuning(:,ft) > lim));
        to_plot = squeeze(mean_func(all_metrics.(type).invariance_nc(vxs_roi,g),1))';
        scatter(xpos(r)+1,to_plot,length(find(vxs_roi))*0.1 ...
            ,cmaps.roi_colors(r,:),'o','LineWidth',2);

        title([groups_names{g} ' low vs high' feats{ft}])
        xticks(xpos(1:n_rois))
        xticklabels(rois)
        xlim([0 xpos(n_rois)+2])
       
        ylabel(['Corr mix ' groups_names{g}])
        
    end
end


if ~ isempty(ext)
    saveas(gcf,[figures_path 'InvariancebyTuning_by_ROI_filter_' species  save_suffix '_' type '_' feats{ft} '.' ext])
end


%% STATS: Test for differences in invariance between voxels tuned to low or high rates, with shuffling by animal

dtype = 'true';
metrics = all_metrics.(dtype);


ft = 3;
tuning = all_metrics.true.tuning(:, ft);
for r = 1:n_rois
    figure('Position', [1.1103e+03 191.6667 1.1927e+03 1086]);
    shuffled_invariance_by_animal = nan(n_perms, length(a_range), 2);
    invariance_by_animal = nan(length(a_range), 2);
    for a = 1:length(a_range)

        vxs_an = logical(ismember(roi_idx , [r]) & vxs_keep &(an_idx==a_range(a)));%.*vxs_keep
        animal_invariance = metrics.invariance_nc(vxs_an,:);
        true_tuning = tuning(vxs_an);

        invariance_by_animal(a,: ) = mean_func(animal_invariance(true_tuning>lim,:),1) - mean_func(animal_invariance(true_tuning<lim,:),1);

        for sf = 1:n_perms
            fake_tuning = true_tuning(randperm(length(true_tuning)));
            shuffled_invariance_by_animal(sf, a, :) = mean_func(animal_invariance(fake_tuning>lim,:),1) - mean_func(animal_invariance(fake_tuning<lim,:),1);

        end
        shuffled_invariance_by_animal(end, a, :)  = invariance_by_animal(a, :);

        for cd = 1:2
            subplot(length(a_range)+1,2,cd+(a-1)*2)
            histogram(shuffled_invariance_by_animal(:,a,cd))
           
            pval = sum(abs(shuffled_invariance_by_animal(:, a, cd)) >= abs(invariance_by_animal(a, cd)))/(n_perms+1);

            if invariance_by_animal(cd) > 0
                disp([dtype ' ' groups_names{cd} ': ' rois{r2} ' > ' rois{r1} ', ' a  ': p=' num2str(pval)]);
            else
                disp([dtype ' ' groups_names{cd} ': ' rois{r2} ' < ' rois{r1} ', ' a  ': p=' num2str(pval)]);
            end
            vline(invariance_by_animal(a, cd))
            text(invariance_by_animal(a, cd),20,['Md= ' num2str(invariance_by_animal(a, cd),3)])
            title([a ', ' groups_names{cd} ', p = ' num2str(pval,3)])
        end
    end

    shuffled_invariance_total = nanmean(shuffled_invariance_by_animal,2);
    invariance_total = nanmean(invariance_by_animal,1);

    for cd = 1:2
        subplot(length(a_range)+1,2,cd+(a)*2)
        histogram(shuffled_invariance_total(:,cd))
        pval = sum(abs(shuffled_invariance_total(:, cd)) >= abs(invariance_total(cd)))/(n_perms+1);
        if invariance_total(cd) > 0
            disp([dtype ' ' groups_names{cd} ': high > low rates for' rois{r} ' ' a ': p=' num2str(pval)]);
        else
            disp([dtype ' ' groups_names{cd} ': high > low rates for ' rois{r} ' ' a  ': p=' num2str(pval)]);
        end
        vline(invariance_total(cd))
        text(invariance_total(cd),20,['Md= ' num2str(invariance_total(cd),3)])
        title(['All, ' groups_names{cd} ', p = ' num2str(pval,3)])
    end
    sgtitle([rois{r} ': low vs high rates'])
end

%% STATS: Test for differences in invariance across ROIS in voxels tuned to low or high rates

r1 = 1; % region 1
r2 = [2]; % region 2
n_perms = 1000;
ft = 3; 
lim = 5;
groups_to_test ={ all_metrics.true.tuning(:,ft) < lim, all_metrics.true.tuning(:,ft) > lim};
groups_to_test_names = {'low rates', 'high rates'};
tables = [];
for gr = 1:length(groups_to_test)
    disp(groups_to_test_names{gr})
    dtype = 'true';
    metrics = all_metrics.(dtype);   
    figure('Position', [1.1103e+03 191.6667 1.1927e+03 1086]);
    shuffled_invariance_by_animal = nan(n_perms, length(a_range), 2);
    invariance_by_animal = nan(length(a_range), 2,2);
    diff_invariance_by_animal = nan(length(a_range),2);
    
    for a = 1:length(a_range)

        vxs_an = logical(ismember(roi_idx , [r1, r2]) & vxs_keep & groups_to_test{gr} &(an_idx==a_range(a)));%.*vxs_keep
        animal_invariance = metrics.invariance_nc(vxs_an,:);
        true_roi = roi_idx(vxs_an);

        val_r2 = mean_func(animal_invariance(ismember(true_roi, r2),:),1) ;
        val_r1 = mean_func(animal_invariance(true_roi == r1,:),1);
        invariance_by_animal(a,: ,1) = val_r2; 
        invariance_by_animal(a,: ,2) = val_r1; 
        diff_invariance_by_animal(a,:) =  invariance_by_animal(a, :,1)- invariance_by_animal(a, :,2);
        for sf = 1:n_perms
            fake_roi = true_roi(randperm(length(true_roi)));
            shuffled_invariance_by_animal(sf, a, :) = mean_func(animal_invariance(ismember(fake_roi, r2),:),1) - mean_func(animal_invariance(fake_roi == r1,:),1);
         
        end

        shuffled_invariance_by_animal(end, a, :)  = diff_invariance_by_animal(a,:);

        for cd = 1:2
            subplot(length(a_range)+1,2,cd+(a-1)*2)
            histogram(shuffled_invariance_by_animal(:,a,cd))

            pval = sum(abs(shuffled_invariance_by_animal(:, a, cd)) >= abs(diff_invariance_by_animal(a, cd)))/(n_perms+1);
          
            disp([dtype ' ' groups_names{cd} ' ' subj_ids{a} ])
            disp([rois{r2} ': ' num2str(invariance_by_animal(a, cd,1),2) ', ' rois{r1} ': ' num2str(invariance_by_animal(a, cd,2),2) ', p=' num2str(pval)]);
          
            vline(diff_invariance_by_animal(a, cd))
            text(diff_invariance_by_animal(a, cd),20,['Md= ' num2str(diff_invariance_by_animal(a, cd),3)])
            title([subj_ids{a} ', ' groups_names{cd} ', p = ' num2str(pval,3)])
        end

    end

    shuffled_invariance_total = nanmean(shuffled_invariance_by_animal,2);
    invariance_total = nanmean(diff_invariance_by_animal,1);

    for cd = 1:2
        subplot(length(a_range)+1,2,cd+(a)*2)
        histogram(shuffled_invariance_total(:,cd))
       
        pval = sum(abs(shuffled_invariance_total(:, cd)) >= abs(invariance_total(cd)))/(n_perms+1);
          
        disp([dtype ' ' groups_names{cd} ' all' ])
        disp([rois{r2} ': ' num2str(nanmean(invariance_by_animal(:, cd,1)),2) ', ' rois{r1} ': ' num2str(nanmean(invariance_by_animal(:, cd,2)),2) ', p=' num2str(pval)]);
       
        vline(invariance_total(cd))
        text(invariance_total(cd),20,['Md= ' num2str(invariance_total(cd),3)])
        title(['ALL , ' groups_names{cd} ', p = ' num2str(pval,3)])
    end

    permute(round(invariance_by_animal,2),[3,2,1])
    sgtitle(['ROIs: ' rois{r2} ' - ' rois{r1} ' for voxels tuned to ' groups_to_test_names{gr}])
    saveas(gcf,[figures_path  'permutation_test_rois_by_animal_shuffle_animal_tuning' species '_' num2str(r1) '_' num2str(r2) '_' dtype '.png'])

end

