% Time Parameters
st.evnt_lab = 'R';
st.stat_lim = [-0.5 1.0];
st.alpha    = 0.05;

% Sliding Window Parameters
st.win_len  = 0.2;%200;
st.win_step = 0.05;%50;

% ANOVA Parameters
st.model_lab   = 'crRT_CNI_PC';
st.regress_rt  = 1;    % regrees reaction time off before running ANOVA
st.groups      = {'CNI', 'PC'};
st.n_boots     = 1000;

% RT Correlation Parameters
st.rt_corr       = 1;
cfg_rt = [];
cfg_rt.parameter        = 'powspctrm';
cfg_rt.statistic        = 'ft_statfun_correlationT';
cfg_rt.method           = 'montecarlo';
cfg_rt.numrandomization = st.n_boots;
cfg_rt.correctm         = 'cluster';
cfg_rt.clusteralpha     = 0.05;   %threshold for a single comparison (time point) to be included in the clust
cfg_rt.clusterstatistic = 'maxsum';
cfg_rt.clustertail      = 0;
cfg_rt.tail             = 0; %two sided
cfg_rt.correcttail      = 'alpha'; %correct the .alpha for two-tailed test (/2)
cfg_rt.computestat      = 'yes';
cfg_rt.computeprob      = 'yes';
cfg_rt.alpha            = st.alpha;
cfg_rt.neighbours       = [];

