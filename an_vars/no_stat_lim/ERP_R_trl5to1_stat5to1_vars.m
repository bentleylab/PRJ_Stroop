error('needs update to move this to stat_vars');
an.evnt_lab    = 'R';              % event around which to cut trials
an.trial_lim_s = [-0.5 1];         % window in SEC for cutting trials
an.demean_yn   = 'yes';
an.bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction
an.stat_lim    = [0 1];            % window in SEC for stats
an.n_boots     = 1000;             % Repetitions for non-parametric stats
an.lp_yn       = 'yes';
an.lp_freq     = 20;
an.hp_yn       = 'yes';
an.hp_freq     = 0.5;

cfg_stat = [];
cfg_stat.latency          = stat_lim;
cfg_stat.channel          = 'all';
cfg_stat.parameter        = 'trial';
cfg_stat.method           = 'montecarlo';
cfg_stat.statistic        = 'ft_statfun_indepsamplesT';
cfg_stat.correctm         = 'cluster';
cfg_stat.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.clustertail      = 0;
cfg_stat.tail             = 0; %two sided
cfg_stat.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_stat.alpha            = 0.05;
cfg_stat.numrandomization = n_boots;
cfg_stat.neighbours       = [];%neighbors;
% cfg_stat.minnbchan        = 0;
cfg_stat.ivar             = 1;  %row of design matrix containing independent variable
% cfg_stat.uvar             = 2;  %row containing dependent variable, not needed for indepsamp

