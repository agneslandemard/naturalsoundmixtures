%% define a few useful parameters

% time window for analysis
sd_length = 9.6;
resp_wins = -5 : sd_length * 2.5 + 5;

% sound categories
n_sds = 36;
stims.fgs_categories = {'ferret', 'speech', 'music', 'other animals', 'random'};
stims.fgs_category = [4 1 1 1 1 4 4 4 1 1 1 4 5 4 3 5 4 5 5 2 3 3 3 5 3 3 5 2 2 2 2 2 2 3 5 5];
stims.events = [0 sd_length/2 sd_length];

% grouping into broad categories
groups = {1, [4 5], [2 3]};
groups_names = {'fg','bg','mix'};

% animals
animals = {'B','R','L'};
n_animals = length(animals);

% colors
clrs =  [0.949    0.604	0.722  ;...
    0.1882    0.7490    0.6196 ;...
    0.0549    0.3020    0.5843 ;...
    0.7407    0.3447    0.2843;...
    0.5    0.45    0.3];
cmaps.bg = cbrewer('qual','Set2',n_sds);
cmaps.fg = clrs(stims.fgs_category,:);
cmaps.cm = permute(cat(3,[69 170 153; 169 204 161],[159 73 141;222 121 166]),[1 3 2])./255;
cmaps.roi_colors = [230 149 50; 159 190 67; 41 99 56]./255;


