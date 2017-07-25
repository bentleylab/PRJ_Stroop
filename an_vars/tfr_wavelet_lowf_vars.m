event_type  = 'stim';           % event around which to cut trials
trial_lim_s = [-0.5 2.5];       % window in SEC for cutting trials
demean_yn   = 'yes';
bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction
stat_lim    = [0 2];            % window in SEC for stats
n_boots     = 1000;             % Repetitions for non-parametric stats

foi_center = [2:0.5:6 7:20 25:5:60];
%%2.^[2:1/4:8]; % Center frequencies
%octave = 3/4;              % Frequency resolution
%foi_min = 2^(-octave/2)*foi_center;
%foi_max = 2^(octave/2)*foi_center;
%foi = (foi_min+foi_max)/2;
%delta_freq = foi_max-foi_min;
%delta_time = 0.5;
%n_taper_all = max(1,round(delta_freq.*delta_time-1));
%foi_center  = round(foi_center*10)/10;
%delta_freq_true = (n_taper_all+1)./delta_time; % total badn width around

cfg_tfr = [];
cfg_tfr.output       = 'pow';
cfg_tfr.channel      = 'all';
cfg_tfr.method       = 'wavelet';
cfg_tfr.width        = 4;
cfg_tfr.pad          = 'maxperlen'; %add time on either side of window
cfg_tfr.padtype      = 'zero';
cfg_tfr.foi          = foi_center;%2:5:150;                         % analysis 2 to 30 Hz in steps of 2 Hz 
%cfg_tfr.t_ftimwin    = ones(length(cfg_tfr.foi),1).*0.5;   % length of time window = 0.5 sec
cfg_tfr.toi          = 'all';%-buff_lim(1):0.1:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cfg_tfr.keeptrials   = 'yes';
% cfg.t_ftimwin    = ones(1,length(cfg.tapsmofrq))*delta_time;
% cfg.width        = 7;
% cfg.gwidth       = 3;

