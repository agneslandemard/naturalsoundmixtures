%% Compare predictions accuracy and its relationship to other metrics, across species

all_metrics = struct();

all_metrics.true_f = load_variables('real');
all_metrics.pred_f = load_variables('predicted_by_cm');
        
all_metrics.true_h = load_variables_human('real');
all_metrics.pred_h = load_variables_human('predicted_by_cm');

[vals, feats] = LoadCMparams;
n_feats = length(feats);


%% Compare overall accuracy 

metric_to_compare = 'pred_accuracy'; % here and elsewhere could compare other metrics 
preds = {'pred_f', 'pred_h'}; % choose what to compare among true/predicted/ferret/human
figure; hold all
for pred_name = preds
    v = snm(all_metrics.(pred_name{1}).(metric_to_compare), 2);
    histogram(v(~isnan(v)),'BinWidth', 0.1,'Normalization','probability');
end
xlim([-0.5,1])
legend(preds)
saveas(gcf,[figures_path 'global_comparison_' metric_to_compare '.' ext])

%% Compare overall accuracy by ROI

metric_to_compare = 'pred_accuracy';

figure('Position',[1000 661.6667 394.3333 676.3333])
subplot(211); hold all

v = snm(all_metrics.pred_f.(metric_to_compare),2);
for r = 1:3
    v_roi = v(all_metrics.pred_f.roi_idx==r);
    histogram(v_roi,'BinWidth', 0.1,'Normalization','probability','FaceColor',cmaps.roi_colors(r,:));
end
xlim([-0.7,1])
subplot(212); hold all
v = snm(all_metrics.pred_h.(metric_to_compare),2);
for r = 1:2
    v_roi = v(all_metrics.true_h.roi_idx==r);
    histogram(v_roi,'BinWidth', 0.1,'Normalization','probability','FaceColor',cmaps.roi_colors(r,:));
end
xlim([-0.7,1])
xlabel('Prediction accuracy')
ylabel('Voxel density')
saveas(gcf,[figures_path  'comparison_by_roi_' metric_to_compare '.' ext])


%% Compute overall accuracy by subject and compare across species

metric_to_compare = 'pred_accuracy';

v = snm(all_metrics.pred_f.(metric_to_compare),2);
acc_f = [];
subj_idx_f = unique(all_metrics.true_f.an_idx);
for a = 1:length(subj_idx_f)
    acc_f = [acc_f, nanmedian(v(all_metrics.true_f.an_idx==subj_idx_f(a)))];
end


v = snm(all_metrics.pred_h.(metric_to_compare),2);
acc_h = [];
subj_idx_h = unique(all_metrics.true_h.an_idx);
for a = 1:length(subj_idx_h)
    acc_h = [acc_h, nanmedian(v(all_metrics.true_h.an_idx==subj_idx_h(a)))];
end

pval = ranksum(acc_h, acc_f);

%% Compare prediction accuracy to test retest

figure ('Position', [1000 925 1.1277e+03 413],'Renderer','Painters')

subplot(121)
title('Ferrets')
binscatter(snm(all_metrics.pred_f.(metric_to_compare),2),snm(all_metrics.true_f.test_retest,2),50,'XLimits',[-0.6,0.9])
xlabel('Prediction accuracy')
ylabel('Average test-retest')
xlim([-0.7,1])
subplot(122)
title('Humans')
binscatter(snm(all_metrics.pred_h.(metric_to_compare),2),snm(all_metrics.true_h.test_retest,2),50,'XLimits',[-0.6,0.9])
xlabel('Prediction accuracy')
ylabel('Average test-retest')
linkaxes
xlim([-0.7,1])
saveas(gcf,[figures_path  'comparison_' metric_to_compare '_vs_test_retest.' ext])


