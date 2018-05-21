event_type  = 'stim';           % event around which to cut trials
% trial_lim_s will be expanded in SBJ09a by t_ftimwin/2 on front and back to avoid NaNs within real trial_lim_s
%   cfg_tfr.method = 'mtmconvol': expand by t_ftimewin/2
%   cfg_tfr.method = 'wavelet: expand by cfg_tfr.width/foi_center(1)*2
trial_lim_s = [-0.3 1];       % window in SEC for cutting trials
%plt_lim     = [-0.2 1];         % window for plotting these results

% ERP vars
erp_vars.demean_yn   = 'yes';
erp_vars.bsln_evnt   = 'stim';
erp_vars.bsln_type   = 'demean';
erp_vars.bsln_lim    = [-0.25 -0.05];    % window RELATIVE TO BSLN_EVNT in SEC for baseline correction
erp_vars.lp_yn       = 'yes';
erp_vars.lp_freq     = 10;
erp_vars.hp_yn       = 'yes';
erp_vars.hp_freq     = 0.5;


stat_lim    = [0 1];            % window in SEC for stats
n_boots     = 1000;             % Repetitions for non-parametric stats
actv_win    = 0.05;               % minimum window in sec for consecutive difference from null distribution

% cfg_stat = [];
% cfg_stat.latency          = stat_lim;
% cfg_stat.channel          = 'all';
% cfg_stat.parameter        = 'trial';
% cfg_stat.method           = 'montecarlo';
% cfg_stat.statistic        = 'ft_statfun_indepsamplesT';
% cfg_stat.correctm         = 'cluster';
% cfg_stat.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
% cfg_stat.clusterstatistic = 'maxsum';
% cfg_stat.clustertail      = 0;
% cfg_stat.tail             = 0; %two sided
% cfg_stat.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
% cfg_stat.alpha            = 0.05;
% cfg_stat.numrandomization = n_boots;
% cfg_stat.neighbours       = [];%neighbors;
% % cfg_stat.minnbchan        = 0;
% cfg_stat.ivar             = 1;  %row of design matrix containing independent variable
% % cfg_stat.uvar             = 2;  %row containing dependent variable, not needed for indepsamp

