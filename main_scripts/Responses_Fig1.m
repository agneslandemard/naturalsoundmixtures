%% Load data for each animal and keep average resp

animals = {'B','R','L'};

an_idx = []; 
bl = 2 : 11; 
n_snips = 6*length(bl)*2;    
keep_data_mat = 0;
av_resp = nan(35,n_snips,length(groups),length(animals));
for a = 1:length(animals)

    animal = animals{a};
    LoadData
    groups = {[1 1],[4 5],[2 3]};
    for g = 1:length(groups)
        av_resp(:,:,g,a) =  reshape(snm(patterns(:,:,:,bl,groups{g},:),[1 6]),length(resp_wins),n_snips);  
    end
    clear patterns
end

%% Plot average response timecourse for each 

groups_names = {'fg','bg','mix'};
clrs = cat(1,squeeze(cmaps.cm(1,1:2,:)),snm(cmaps.cm(1,1:2,:),2)');
figure; 
hold all
x = resp_wins./param.exp.SR;
y = snm(av_resp,[2 4]);
e = std(snm(av_resp,4),[],2)./sqrt(n_snips); % error across sounds
for g = 1:length(groups_names)
    boundedline(x,y(:,g)*100,e(:,g)*100,'cmap',clrs(g,:),'alpha')
end
title('All ferrets')
xlabel('Time from background change (s)')
xlim(([resp_wins(1) resp_wins(end)])/param.exp.SR)
ylabel('%CBV')
vline(stims.events,'k:')
saveas(gcf,mkpdir([figures_path 'AvResp.' ext]));

figure('Position', [325 796 188 481])
for a = 1:length(animals)
    x = resp_wins./param.exp.SR;
    y = snm(av_resp(:,:,:,a),[2 4]);
    e = nanstd(snm(av_resp(:,:,:,a),4),[],2)./sqrt(n_snips); % error across sounds

    for g = 1:length(groups_names)
        subplot(length(animals),1,a)
        hold all
        boundedline(x,y(:,g)*100,e(:,g)*100,'cmap',clrs(g,:),'alpha')
    end

    xlabel('Time from background change(s)')
    xlim(([resp_wins(1) resp_wins(end)])/param.exp.SR)
    ylabel('%CBV')
    vline(stims.events,'k:')
    title(animals{a})
end
saveas(gcf,mkpdir([figures_path 'AvResp_byAnimal.' ext]));


