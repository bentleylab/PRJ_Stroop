% Parameters for single trial response onset latencies using linear method
rol_buff     = [0 0.7]; % buffer around stim and RT (stim onset to max(RT)+0.7)
quant_thresh = 0.75;    % cut off in distribution of power for initial start point
rol_lim_s = [-0.2 0.1];   % window around peak to find onset latency (Foster 2015 = [-0.2 0.1])

% Filter settings
sgfilt_ord = 9;
sgfilt_win = 501;

%   (Foster2015 uses 100ms, Bartoli HBM uses 40ms)
reg_win_len_s = 0.05;    % window for regression method in sec
reg_win_stp_s = 0.01;
n_big_slopes  = 5;      % number of biggest slopes to examine

% Plotting Parameters:
trl_plt_perc = 0.05; %
%deriv_scale  = 20;      % scaling factor for plotting the derivative
clim_perc    = [5 95];    % percentile of power for color limits
deriv_marker = '+';
lin_marker   = 'o';
rt_marker    = '*';

