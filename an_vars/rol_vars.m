% Parameters for single trial response onset latencies using linear method
rol_buff     = [0 0.7]; % buffer around stim and RT (stim onset to max(RT)+0.7)
quant_thresh = 0.75;    % cut off in distribution of power for initial start point
%   (Foster2015 uses 100ms, Bartoli HBM uses 40ms)
win_len      = 0.05;    % window for regression method in sec
rol_lim_s = [-0.2 0.1];   % window around peak to find onset latency (Foster 2015 = [-0.2 0.1])
