event_type  = 'stim';           % event around which to cut trials
trial_lim_s = [-0.5 2.5];       % window in SEC for cutting trials
demean_yn   = 'yes';
bsln_lim    = [-0.25 -0.05];    % window in SEC for baseline correction
stat_lim    = [0 2];            % window in SEC for stats
n_boots     = 1000;             % Repetitions for non-parametric stats

foi_center  = [2:0.5:6 7:16 18 20 25 30 40 50 60 80 100 125 150 175 200 250];%2.^[2:1/4:8]; % Center frequencies
octave      = 3/4;              % Frequency resolution
foi_min     = 2^(-octave/2)*foi_center;
foi_max     = 2^(octave/2)*foi_center;
foi         = (foi_min+foi_max)/2;
delta_freq  = foi_max-foi_min;
delta_time  = 0.5;
n_taper_all = max(1,round(delta_freq.*delta_time-1));   %number of tapers for each frequency
foi_center  = round(foi_center*10)/10;          %convert to float?
delta_freq_true = (n_taper_all+1)./delta_time; % total bandwidth around

cfg_tfr = [];
cfg_tfr.output       = 'pow';
cfg_tfr.channel      = 'all';
cfg_tfr.method       = 'mtmconvol';
cfg_tfr.taper        = 'dpss';
cfg_tfr.tapsmofrq    = delta_freq_true./2;                  %ft wants half bandwidth around the foi
cfg_tfr.keeptapers   = 'no';
cfg_tfr.pad          = 'maxperlen';                         %add time on either side of window
cfg_tfr.padtype      = 'zero';
cfg_tfr.foi          = foi_center;%2:5:150;                 % analysis 2 to 30 Hz in steps of 2 Hz 
cfg_tfr.t_ftimwin    = ones(length(cfg_tfr.foi),1).*0.5;    % length of time window; 0.5 sec, could be n_cycles./foi for n_cylces per win
cfg_tfr.toi          = 'all';%-buff_lim(1):0.1:1.5;         % time window centers
cfg_tfr.keeptrials   = 'yes';                               % must be 'yes' for stats
% cfg.t_ftimwin    = ones(1,length(cfg.tapsmofrq))*delta_time;

