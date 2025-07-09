%% Load data

metrics_true = load_variables('real',struct('remove_nans', false));
metrics_true_no_filt = load_variables('real', struct('filter_type','none'));
metrics_pred = load_variables('predicted_by_cm',struct('remove_nans', false));

%% Plotting
% Plot maps to put show voxelwise metrics back in the original voxel space
% ('3D') or in a flattened/surface view ('surface')

% plotting parameters, to choose
plot_version = '3D'; % 'surface' for top view or '3D' for all full slices

toplotlist = {{'test_retest','av_resp','rois'},...
    {'invariance_nc_fg','invariance_nc_bg'},...
    {'invariance_nc_fg_pred','invariance_nc_bg_pred','pred_accuracy'},...
    {'freqs', 'rates', 'scales'}};


example_hemi = 'BLL';
example_slices = [4, 5, 10];

% Colormaps for model features
[vals, feats] = LoadCMparams;
n_feats = length(feats);
colormaps_feats.freqs =  cbrewer('div','Spectral',length(vals.freqs));
colormaps_feats.rates =  cbrewer('div','Spectral',length(vals.rates));
colormaps_feats.scales =  cbrewer('div','Spectral',length(vals.scales));

f1 = figure;
if strcmp(plot_version, '3D')
    f2 = figure;
    set(f2,'Position',[271 572 402 276])
end
for tp =  1:length(toplotlist)
    plotWhat = toplotlist{tp};
    for tpp = 1:length(plotWhat)
        what = plotWhat{tpp};

        Color = [];
        switch what
            case 'test_retest'
                metric = snm(metrics_true_no_filt.test_retest(:,3:4),2);
                m = 0.5;
                Color.ColorAxis = m * [-1 1];
                Color.alpha = 0;
                Color.cm = cbrewer('div','RdBu',64);
                Color.cm = Color.cm(end:-1:1,:);
            case 'av_resp'
                metric = snm(metrics_true_no_filt.av_resp,[2 3 4]);
                m =0.1;
                Color.ColorAxis = m * [-1 1];
                Color.alpha = 0;
                Color.cm = cbrewer('div','RdBu',64);
                Color.cm = Color.cm(end:-1:1,:);
            case feats
                n_ft = find(strcmp(feats,what));
                metric = metrics_pred.tuning(:,n_ft);
                Color.ColorAxis = [1 length(vals.(what))];
                Color.cm = colormaps_feats.(what);
                Color.cm = Color.cm(end:-1:1,:);
                metric(isnan(snm(metrics_true.test_retest,[2 3]))) = nan;
                Color.alpha = 1;
            case 'rois'
                metric = metrics_true.roi_idx;
                Color.ColorAxis = [1 max(metrics_true.roi_idx)];
                switch plot_version
                    case 'surface'
                        Color.cm = cmaps.roi_colors;
                    case '3D'
                        Color.cm = cat(1,[0 0 0],cmaps.roi_colors);
                end
                metric(isnan(snm(metrics_true.test_retest,[2 3]))) = nan;

            case 'invariance_nc_fg'
                metric = metrics_true.invariance_nc(:,1);
                vxs_good_pred = all(~isnan(metrics_true.invariance_nc),2);
                metric(~vxs_good_pred) = nan;
                Color.cm = cmap_from_name('lightblue-to-yellow1');
                Color.ColorAxis = [0, 1];
                Color.alpha = 1;
            case 'invariance_nc_bg'
                metric = metrics_true.invariance_nc(:,2);
                vxs_good_pred = all(~isnan(metrics_true.invariance_nc),2);
                metric(~vxs_good_pred) = nan;
                Color.alpha = 1;
                Color.ColorAxis = [0, 1];
                Color.cm = cmap_from_name('lightblue-to-yellow1');
            case 'invariance_nc_fg_pred'
                metric = metrics_pred.invariance_nc(:,1);
                vxs_good_pred = all(~isnan(metrics_true.invariance_nc),2);
                metric(~vxs_good_pred) = nan;
                Color.alpha = 0;
                Color.cm = cmap_from_name('lightblue-to-yellow1');
                Color.ColorAxis = [0, 1];
                Color.alpha = 1;

            case 'invariance_nc_bg_pred'
                metric = metrics_pred.invariance_nc(:,2);
                vxs_good_pred = all(~isnan(metrics_true.invariance_nc),2);
                metric(~vxs_good_pred) = nan;
                Color.cm = cmap_from_name('lightblue-to-yellow1');
                Color.ColorAxis = [0, 1];
                Color.alpha = 1;

            case 'pred_accuracy'

                metric = metrics_pred.pred_accuracy;
                vxs_good_pred = all(~isnan(metrics_true.invariance_nc),2);
                metric(~vxs_good_pred) = nan;
                m =0.6;
                Color.ColorAxis = m * [-1 1];

            otherwise
                metric = metrics_true.(what);


        end
        m = 1;

        % default color parameters
        if ~isfield(Color,'ColorAxis')
            Color.ColorAxis = m * [-1 1];
        end
        if ~isfield(Color,'cm')
            Color.cm = cbrewer('div','RdBu',64);
            Color.cm = Color.cm(end:-1:1,:);

            n = 0.8; N = 17;

            b = Color.cm(N,:);
            r = Color.cm(end-N,:);
            Color.cm = cat(1,Color.cm(1:N,:),linspaceNDim(b,n*[1 1 1],(32-N))',linspaceNDim(n*[1 1 1],r,(32-N))',Color.cm(end-N:end,:));
        end
        if ~isfield(Color,'alpha')
            Color.alpha = 1;
        end

        switch plot_version
            case 'surface'

                f1 = figure;
                for a = 1:n_animals
                    XA = load([metrics_path  '/anat_params_' animals{a} '.mat']);
                    toPlotan = metric(metrics_true_no_filt.an_idx == a);
                    for h = 1 : length(XA.hemis)
                  
                        r = double(strcmp(XA.hemis{h}(3),'L')) + 1; % 1 if left, 2 if right
                        figure(f1)
                        subplot(2,n_animals,(r-1)*n_animals + a)
                        toPlot = toPlotan(XA.si == h);
                        PlotTopViewBasic(toPlot,XA.(XA.hemis{h}).param,XA.(XA.hemis{h}).Anat,Color);
                        title(XA.hemis{h})

                    end
                end
                sgtitle(f1,regexprep(what,'_',' '))
                saveas(f1,[figures_path 'Maps_' what '_' plot_version '.' ext])
                

            case '3D'
                offset = 0;
                Anat = nan(100,100,20,n_animals*2);
                W = nan(100,100,20,n_animals*2); 
                ct = 0; 
                list_hemis = [];

                for a = 1:n_animals
                    XA = load([metrics_path  '/anat_params_' animals{a} '.mat']);

                    toPlotan = metric(metrics_true_no_filt.an_idx == a);
                    tmph = cell(1,2);
                    for h = 1:length(XA.hemis)
                        r = double(strcmp(XA.hemis{h}(3),'R')) + 1; % 1 if left, 2 if right

                        toPlot = Pixs2Mat(toPlotan(XA.si == h),XA.(XA.hemis{h}).param.msk);
                        sz = size(toPlot);

                        W(1:sz(1),1:sz(2),offset + (1:sz(3)),(a-1)*2+r) = toPlot;
                        Anat(1:sz(1),1:sz(2),offset + (1:sz(3)),(a-1)*2+r) = XA.(XA.hemis{h}).Anat;

                        tmph{r} = regexprep(XA.hemis{h},'_1','');
                    end
                    list_hemis = [list_hemis, tmph];
                end
                Anat(all(isnan(W),[2 3 4]),:,:,:)=[];
                W(all(isnan(W),[2 3 4]),:,:,:)=[];

                Anat(:,all(isnan(W),[1 3 4 ]),:,:)=[];
                W(:,all(isnan(W),[1 3 4 ]),:,:)=[];

                Anat(:,:,all(isnan(W),[1 2 4]),:)=[];
                W(:,:,all(isnan(W),[1 2 4]),:)=[];
                sz = size(W);

                Wc = W(:,:,:,strcmp(list_hemis,example_hemi));
                Wc = Wc(:,:,example_slices);
                szc = size(Wc);
                Wc = reshape(Wc,szc(1),szc(2)*szc(3));

                Anatc = Anat(:,:,:,strcmp(list_hemis,example_hemi));
                Anatc = reshape(Anatc(:,:,example_slices),szc(1),szc(2)*szc(3));
                Wc(:,all(isnan(Anatc),1))= [];
                Anatc(:,all(isnan(Anatc),1))= [];


                W = reshape(permute(W,[1 3 2 4]),sz(1)*sz(3),sz(2)*sz(4));

                Anat = reshape(permute(Anat,[1 3 2 4]),sz(1)*sz(3),sz(2)*sz(4));
                Anat(:,all(isnan(W),1))= [];
                W(:,all(isnan(W),1))= [];

                figure(f1)
                subplot(1,length(plotWhat),tpp)
                PlotMapWithAnat(W,sqrt(Anat),Color);
                list_hemis(cellfun(@(x) isempty(x),list_hemis)) = [];
                xlabel(regexprep(cell2str(list_hemis),'"',""))
                ylabel('A -> P')

                figure(f2)
                subplot(length(plotWhat),1,tpp)
                PlotMapWithAnat(Wc,sqrt(Anatc),Color);
                sgtitle(what)
        end

    end

    if strcmp(plot_version,'3D')
        name = regexprep(cell2str(plotWhat),"[',{,},;,]","");
        name = regexprep(name,' ','_');
        sgtitle(f1,name)
        saveas(f1,[figures_path 'Maps_' name '_' plot_version '.' ext])
        saveas(f2,[figures_path 'Maps_' name '_' plot_version '_ex' example_hemi '.' ext])
    end
end
