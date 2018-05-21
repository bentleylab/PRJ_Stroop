stat_vars.only_actv_ch = 1;             % Limit analysis to only activ channels
stat_vars.actv_epochs_overlap = 0;      % Limit correlation pairs to channels that have overlapping active periods? 0/1
stat_vars.an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';    % R-locked analysis to limit to active channels
stat_vars.an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';  % S-locked analysis to limit to active channels
stat_vars.actv_win = 100;
stat_vars.event_lab = 'stim';
stat_vars.stat_lim = [0 2.5];
stat_vars.conditions = 'CI';    % just for plotting RTs in S-locked

stat_vars.lp_flag = 1;
stat_vars.lp_freq = 6;
stat_vars.hp_flag = 0;
stat_vars.hp_freq = [];

stat_vars.n_boots = 1000;
%now in plt_vars: stat_vars.sig_cut = 0.05;
