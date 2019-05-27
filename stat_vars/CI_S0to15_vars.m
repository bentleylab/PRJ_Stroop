st.conditions = 'CI';

% Stats parameters
st.stat_lim    = [0 1.5];            % window in SEC for stats
st.n_boots     = 500;             % Repetitions for non-parametric stats
st.alpha    = 0.05;

cfg_stat = [];
cfg_stat.latency          = st.stat_lim;
cfg_stat.channel          = 'all';
cfg_stat.parameter        = 'powspctrm';
cfg_stat.method           = 'montecarlo';
cfg_stat.statistic        = 'ft_statfun_indepsamplesT';
cfg_stat.correctm         = 'cluster';
cfg_stat.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.clustertail      = 0;
cfg_stat.tail             = 0; %two sided
cfg_stat.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_stat.alpha            = st.alpha;
cfg_stat.numrandomization = st.n_boots;
cfg_stat.neighbours       = [];%neighbors;
% cfg_stat.minnbchan        = 0;
cfg_stat.ivar             = 1;  %row of design matrix containing independent variable
% cfg_stat.uvar             = 2;  %row containing dependent variable, not needed for indepsamp

