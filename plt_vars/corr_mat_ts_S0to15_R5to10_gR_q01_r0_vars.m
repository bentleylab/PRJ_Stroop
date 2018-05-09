% Correlataion Matrix Parameters
plt_vars.sort_vec = [3, 2];  % gROI, ROI
plt_vars.sig_cut = 0.01;
plt_vars.corr_cut = 0;
plt_vars.exclude_FWM = 1;
plt_vars.exclude_OUT = 1;

plt_vars.double_sig_dot = 1;
plt_vars.sig_dot_size = 25;
plt_vars.sig_dot_color = 'w';
plt_vars.sig_dot_size2 = 15;
plt_vars.sig_dot_color2 = 'r';
plt_vars.thin_labels = [1,1];   % x and y labels
plt_vars.plt_maj_div = 1;
plt_vars.plt_min_div = 1;
plt_vars.div_line_width = 1;
plt_vars.font_sz    = 14;
plt_vars.plt_min_div = 1;
plt_vars.subplot_pos = [0 0 0.5 1; 0.5 0 0.5 1];  % [left bottom width height]

plt_vars.fig_dim = [0 0 0.8 1];

% Histogream parameters
plt_vars.y_jitter = 0.1;    % proportion of histogram height to jitter scatter
plt_vars.hist_bins = 25;

% Time Series Plot Parameters
plt_vars.ylim_fudge = 0.1;
plt_vars.y_tick_int = 0.2;
plt_vars.x_step_sz  = 0.2;
plt_vars.legend     = 1;
plt_vars.legend_loc_S = 'northwest';
plt_vars.legend_loc_R = 'northeast';
plt_vars.legend_loc = 'northeast';
plt_vars.errbar_alpha = 0.2;
plt_vars.errbar_color  = [0.5 0.5 0.5];
plt_vars.main_style = {'-','-'};
plt_vars.main_colors = {[0.9290, 0.6940, 0.1250], [0, 0.4470, 0.7410]};

plt_vars.evnt_width = 2;
plt_vars.evnt_color = 'k';
plt_vars.evnt_style = '--';
