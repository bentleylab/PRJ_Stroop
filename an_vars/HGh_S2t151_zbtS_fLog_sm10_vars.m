an.evnt_lab    = 'S';              % event around which to cut trials
% trial_lim_s will NOT be full of data! the first and last t_ftimwin/2 epochs will be NaNs
an.trial_lim_s = [-0.25 1.51];      % window in SEC for cutting trials
an.demean_yn   = 'no';             % z-score for HFA instead
an.bsln_evnt   = 'stim';
an.bsln_type   = 'zboot';
an.bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction
an.bsln_boots  = 500;             % Repetitions for non-parametric stats

% HFA Calculations
an.HFA_type   = 'hilbert';
an.foi_lim = [70 150]; % min and max of desired frequencies
an.n_foi   = 8;
an.min_exp = log(an.foi_lim(1))/log(2); % that returns the exponents
an.max_exp = log(an.foi_lim(2))/log(2);
an.fois    = 2.^[linspace(an.min_exp,an.max_exp,an.n_foi)];
an.foi_bws = fn_semilog_bws(an.fois);     % semilog bandwidth spacing to match Erik Edwards & Chang lab
an.bp_lim  = zeros([numel(an.fois) 2]);
for f = 1:numel(an.fois)
    an.bp_lim(f,:) = fn_freq_lim_from_CFBW(an.fois(f), an.foi_bws(f));
end

cfg_hfa = [];
cfg_hfa.hilbert  = 'abs';
cfg_hfa.bpfilter = 'yes';
cfg_hfa.bpfreq   = [];      % to be filled by looping through foi_center
cfg_hfa.channel  = 'all';

% Cleaning up power time series for plotting
an.smooth_pow_ts = 1;
an.lp_yn       = 'yes';
an.lp_freq     = 10;
an.hp_yn       = 'no';
an.hp_freq     = 0.5;

% Resampling
an.resample_ts   = 0;
% an.resample_freq = 250;

