an.evnt_lab    = 'S';              % event around which to cut trials
% trial_lim_s will be expanded in SBJ09a by t_ftimwin/2 on front and back to avoid NaNs within real trial_lim_s
%   cfg_tfr.method = 'mtmconvol': expand by t_ftimewin/2
%   cfg_tfr.method = 'wavelet: expand by cfg_tfr.width/foi_center(1)*2
% 3 Hz min, 1 cycles = 1s win, 0.5s before 1st point of interest
an.trial_lim_s = [-0.5 2.0];         % window in SEC for cutting trials
an.bsln_evnt   = 'stim';
an.bsln_type   = 'zboot';
an.bsln_boots  = 1000;
an.bsln_lim    = [-0.5 -0.2];    % window in SEC for baseline correction
an.demean_yn   = 'no';

% TFR Calculations
an.foi_center = 2.^((8:40)/8);

cfg_tfr = [];
cfg_tfr.output       = 'pow';
cfg_tfr.channel      = 'all';
cfg_tfr.method       = 'wavelet';
cfg_tfr.width        = 3;
cfg_tfr.pad          = 'nextpow2'; %add time on either side of window
cfg_tfr.padtype      = 'zero';
cfg_tfr.foi          = an.foi_center;  % analysis 2 to 32 Hz, log spaced
cfg_tfr.toi          = 'all';
cfg_tfr.keeptrials   = 'yes';       %must be 'yes' for stats