%% LOAD DATA AND FORMAT IT TO BE READY TO BE ANALYZED
% the output will be a matrix:
% patterns:  (N_voxels x N_timepoints_in_snippet x N_foreground_streams x N_snippets_in_stream x N_stream_combinations x N_repetitions)
% for N_stream_combinations  dimension: 1 = isolated foreground stream, 
% 2 = foreground + background stream 1, 3 = foreground + background stream 2, 
% 4 = background stream 1, 5 = background stream 2.
% See manuscript for description of the stream structure.

%% load data
if ~exist('keep_data_mat','var')
    disp('Will not keep unformatted data')
    keep_data_mat = 0;
end

disp (['Loading data for ' animal])
file_name = [metrics_path 'full_data/Data_' animal '.mat'];
if exist(file_name,'file')
    load(file_name,'norm_data_mat','stims','hemis','param','si','AllSliceRef');
else
    error('Data missing...')
end

LoadStimsInfo

sd_length = 9.6;
resp_wins = -5 : sd_length * param.exp.SR + 5;
[patterns, or_idx] =  FormatIntoPatterns(norm_data_mat, stims, param, resp_wins);
fg_id = stims.fgs(or_idx(:, :, [2 3]));
fg_id_n = arrayfun( @(x) find(strcmp(stims.fgs_list, x)), fg_id) - 1;
bg_id = stims.bgs(or_idx(:, :, [2 3]));
bg_id_n = arrayfun( @(x) find(strcmp(stims.bgs_list, x)),bg_id) -1;

if ~keep_data_mat
    clear norm_data_mat % do relieve memory load
end
