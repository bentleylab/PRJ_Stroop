rol_lab  = {'der','lin'};
evnt_lab = {'S','R'};

% Parameters for single trial response onset latencies using linear method
rol_trl_lim_s = [0 0.4];        % buffer around stim and RT (stim onset to max(RT))
remove_postRT = 1;              % NaN out data past the rol_trl_lim_s(2) on single trial basis
min_actv_s    = 0.1;            % Minimum activation  length (in sec)
rol_lim_s     = [-0.2 0.1];     % window around peak to find onset latency (Foster 2015 = [-0.2 0.1])
quant_thresh  = 0.75;           % cut off in distribution of power for initial start point

% Filter settings
sgfilt_ord = 5;
sgfilt_win = 201;

% Linear Regression Method Settings
%   (Foster2015 uses 100ms, Bartoli HBM uses 40ms)
reg_win_len_s = 0.05;    % window for regression method in sec
reg_win_stp_s = 0.01;
n_big_slopes  = 5;      % number of biggest slopes to examine

% Stimulus/Response modeling
rol_outlier_std = 2;
%rol_outlier_s   = 0.5;
%plot_model_fits = 0;

% PLOTTING:
% QA settings
trl_plt_perc = 0.05; %
%deriv_scale  = 20;      % scaling factor for plotting the derivative
n_hist_bins  = 30;      % for rol window sizes, activation lengths, and activation amplitudes

% Stack settings
clim_perc    = [5 95];    % percentile of power for color limits
evnt_colors  = {'b','r'};
rol_mrkrs    = {'*','o'};
rol_style    = {'--',':'};
rol_colors   = [250 159 181; 221 52 151]./255;    % der = pink; lin = magenta
rol_mrkr_sz  = 10;

rt_marker    = '+';
rt_color     = 'k';
rt_mrkr_sz   = 5;

