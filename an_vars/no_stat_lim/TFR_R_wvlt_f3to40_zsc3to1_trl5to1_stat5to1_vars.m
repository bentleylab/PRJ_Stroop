an.evnt_lab    = 'R';              % event around which to cut trials
% trial_lim_s will be expanded in SBJ09a by t_ftimwin/2 on front and back to avoid NaNs within real trial_lim_s
%   cfg_tfr.method = 'mtmconvol': expand by t_ftimewin/2
%   cfg_tfr.method = 'wavelet: expand by cfg_tfr.width/foi_center(1)*2
an.trial_lim_s = [-0.5 1];         % window in SEC for cutting trials
an.demean_yn   = 'no';
an.bsln_evnt   = 'stim';
an.bsln_type   = 'zscore';
an.bsln_lim    = [-0.3 -0.1];    % window in SEC for baseline correction

% TFR Calculations
an.foi_center = [3:0.5:6 7:20 25:3:40];

cfg_tfr = [];
cfg_tfr.output       = 'pow';
cfg_tfr.channel      = 'all';
cfg_tfr.method       = 'wavelet';
cfg_tfr.width        = 3;
cfg_tfr.pad          = 'maxperlen'; %add time on either side of window
cfg_tfr.padtype      = 'zero';
cfg_tfr.foi          = an.foi_center;  % analysis 2 to 30 Hz in steps of 2 Hz 
cfg_tfr.toi          = 'all';%-buff_lim(1):0.1:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg_tfr.keeptrials   = 'yes';       %must be 'yes' for stats

