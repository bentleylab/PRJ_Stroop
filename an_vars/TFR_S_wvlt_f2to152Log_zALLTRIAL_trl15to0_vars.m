an.evnt_lab    = 'S';              % event around which to cut trials
an.trial_lim_s = [-1.5 1.5];         % window in SEC for cutting trials
an.bsln_evnt   = 'stim';
an.bsln_type   = 'zboot';
an.bsln_boots  = 1000;
an.bsln_lim    = [-1.5 1.5]; 
an.demean_yn   = 'no';

% TFR Calculations
an.foi_center = 2.^((8:58)/8);

cfg_tfr = [];
cfg_tfr.output       = 'pow';
cfg_tfr.channel      = 'all';
cfg_tfr.method       = 'wavelet';
cfg_tfr.width        = 6;
cfg_tfr.pad          = 'nextpow2'; %add time on either side of window
cfg_tfr.padtype      = 'zero';
cfg_tfr.foi          = an.foi_center;  
cfg_tfr.toi          = 'all';
cfg_tfr.keeptrials   = 'yes';       %must be 'yes' for stats