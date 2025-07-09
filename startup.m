global metrics_path figures_path glob_path

glob_path = 'naturalsoundmixtures/'; %% UPDATE THAT TO YOUR PATH

metrics_path = [glob_path 'data/'];

figures_path = mkpdir([glob_path 'figures/']);

addpath(genpath(glob_path))

LoadStimsInfo

% extension to save figures
ext = 'png';


