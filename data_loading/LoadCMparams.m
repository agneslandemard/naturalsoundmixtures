
function [vals, feats] = LoadCMparams()

vals.rates = [0.5000 1 2 4 8 16 32 64 128]; % in Hz
vals.scales = [0 0.2500 0.5000 1 2 4 8]; % cyc/oct
vals.freqs = [27 50  91 167 307 562 1032 1892 3470 6363]; % Hz
feats = {'freqs','scales','rates'};

end