%% Pipeline Processing Variables: "main_ft"

% Data Preprocessing
proc_vars.plot_psd      = '1by1';         % type of plot for channel PSDs
proc_vars.resample_freq = 1000;
proc_vars.demean_yn     = 'yes';
proc_vars.hp_freq       = 0.5;            % [] skips this step
proc_vars.lp_freq       = 300;            % [] skips this step
proc_vars.notch_type    = 'bandstop';     % method for nothc filtering out line noise

% Behavioral Processing
proc_vars.rt_bounds = [0.3 2.0];          % bounds on a reasonable RT to be detected with KLA algorithm

% Trial Cut Parameteres
proc_vars.event_type    = 'stim';         % 'stim'/'resp': lock trial to these event
proc_vars.trial_lim_sec = [-0.25 2];      % data segments (in seconds) to grab around events
proc_vars.RT_std_thresh = 3;              % rejection threshold for RTs

% Varaince-Based Trial Rejection Parameters
proc_vars.var_std_warning_thresh = 3;
