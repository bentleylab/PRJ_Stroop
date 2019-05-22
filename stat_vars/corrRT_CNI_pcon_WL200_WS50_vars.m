% Model Factors and Levels
model_lab   = 'corrRT_CNI_pcon';

% ANOVA Parameters
regress_rt  = 1;    % regrees reaction time off before running ANOVA
groups      = {'CNI', 'pcon'};
%levels      = {{'con','neu','inc'},{'mcon','same','minc'}};

% Stats parameters
n_boots = 1000;

% Sliding Window Parameters
win_len  = 0.2;%200;
win_step = 0.05;%50;

% RT Correlation Parameters
rt_correlation = 1;
cfg_rt = [];
cfg_rt.parameter        = 'powspctrm';
cfg_rt.statistic        = 'ft_statfun_correlationT';
cfg_rt.method           = 'montecarlo';
cfg_rt.numrandomization = n_boots;
cfg_rt.correctm         = 'cluster';
cfg_rt.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_rt.clusterstatistic = 'maxsum';
cfg_rt.clustertail      = 0;
cfg_rt.tail             = 0; %two sided
cfg_rt.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_rt.computestat      = 'yes';
cfg_rt.computeprob      = 'yes';
cfg_rt.alpha            = 0.05;
cfg_rt.neighbours       = [];

