
function [patterns, or_idx] =  FormatIntoPatterns(norm_data_mat, stims, param, resp_wins, correct)

plot_figs = 0;

n_steps = length(resp_wins);

if nargin < 5
    correct = true;
end

% size should be variables (eg vx) * time in bloc * blocs * reps
[n_vxs, n_tps, ~, n_reps] = size(norm_data_mat);

if correct 
    correction_lag = - 1; 
else
    disp('not correcting')
    correction_lag = 0;
end

who = 'fg'; id = 1;
who_other = 'bg';

times = stims.([who 's_times']);
othertimes = stims.([ who_other 's_times']);
fulllist = stims.([who 's']);
other = stims.([ who_other 's']);
bid = stims.blocs_ids(:,id);
oid = stims.blocs_ids(:,setdiff(1:2,id));
list = unique(fulllist(1:12:end));

sd_length = 9.6;
lag = (sd_length/2)*param.exp.SR ;

n_fgs = length(list) - 1;
patterns = nan(n_vxs, n_steps, n_fgs, 12, 5, n_reps);
or_idx = nan(n_fgs, 12, 5);

for sd = 2:length(list)

    % find all occurences of this sound presented in mixture, 
    % align activity to change in background
    sd_mix_idx = find(strcmp(fulllist,list{sd}).*(~strcmp(other,'')));
    if strcmp(param.exp.Exp,'BGFG_v1')
        assert(length(sd_mix_idx) == 4)
        sd_mix_idx = [sd_mix_idx(1) + (0:11), sd_mix_idx(3) + (0:11)];
    end
    assert(all(bid(sd_mix_idx) ~= 0))

    for k = 1:2

        subset = find(oid(sd_mix_idx)==k);
        assert(length(subset)==12)

        kb = unique(oid(sd_mix_idx(subset)));
        krun = unique(stims.runs(sd_mix_idx(subset)));
        sd_solo_idx_bg = find((stims.runs == krun).*(oid == kb).*(bid == 0));
        assert(length(sd_solo_idx_bg)==12)
        
        if plot_figs
            figure;
            plot(0:318,snm(norm_data_mat(:,:,unique(stims.blocs_ids(sd_mix_idx(subset),3)),:),[1,4]))
            vline(times(sd_mix_idx(subset))*param.exp.SR+ correction_lag,'red')
            vline(times(sd_mix_idx(subset))*param.exp.SR+ 12+correction_lag,'green')
            title('Mix')

            figure;
            plot(0:318,snm(norm_data_mat(:,:,stims.blocs_ids(sd_solo_idx_bg(1),3),:),[1, 4]))
            vline((2-k)*12+times(sd_solo_idx_bg)*param.exp.SR+ correction_lag)
            title(['Solo ' num2str(k)])
        end
        
        for s = 1:12

            int = times(sd_mix_idx(subset(s)))*param.exp.SR+resp_wins(1) + lag + correction_lag: times(sd_mix_idx(subset(s)))*param.exp.SR+resp_wins(end)+ lag+ correction_lag;
            [int_t,idx_t] = intersect(round(int),1:n_tps);
            assert(~isempty(int_t))
            patterns(:,idx_t,sd-1,s,k+1,:) = norm_data_mat(:,int_t,stims.blocs_ids(sd_mix_idx(subset(s)),3),:);

            or_idx(sd-1,s,k+1) = sd_mix_idx(subset(s));

            %  same sound in isolation
            int = othertimes(sd_solo_idx_bg(s))*param.exp.SR+resp_wins(1) + correction_lag  : othertimes(sd_solo_idx_bg(s))*param.exp.SR+resp_wins(end) + correction_lag;
            [int_t,idx_t] = intersect(round(int),1:n_tps);

            patterns(:,idx_t,sd-1,s,k+3,:) = norm_data_mat(:,int_t,stims.blocs_ids(sd_solo_idx_bg(s),3),:);

            or_idx(sd-1,s,k+3) = sd_solo_idx_bg(s);
        end
    end


    %find all occurences of this fg sound presented in isolation, time
    %lock to change in bg
    sd_solo_idx = find(strcmp(fulllist,list{sd}).*strcmp(other,''));
    sd_solo_idx = sd_solo_idx(1) + (0:11);
    
    for s = 1:length(sd_solo_idx)

        int = times(sd_solo_idx(s))*param.exp.SR+resp_wins(1) + lag+ correction_lag  : times(sd_solo_idx(s))*param.exp.SR+resp_wins(end) + lag+ correction_lag;
        [int_t,idx_t] = intersect(round(int),1:n_tps);

        patterns(:,idx_t,sd-1,s,1,:) = norm_data_mat(:,int_t,stims.blocs_ids(sd_solo_idx(s),3),:);

        or_idx(sd-1,s,1) = sd_solo_idx(s);
    end
end

end
