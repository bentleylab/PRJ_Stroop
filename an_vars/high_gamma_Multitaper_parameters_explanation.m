foi_center  = [70:10:150];      % Center frequencies that you want, though the exact numbers will change
octave      = 3/4;              % Frequency resolution; e.g., for 40 Hz + 1 octave = 80 Hz, so 3/4 would be 40 +/- 15, or 25:55 Hz
foi_min     = 2^(-octave/2)*foi_center; % move these up by log spacing
foi_max     = 2^(octave/2)*foi_center;
foi         = (foi_min+foi_max)/2;      % find the actual center frequencies after doing that
delta_freq  = foi_max-foi_min;          % find the bandwidths around those actual frequencies
delta_time  = 0.5;                      % time window length (in sec) over which to do the analysis
n_taper_all = max(1,round(delta_freq.*delta_time-1));   % # tapers for each frequency; this formula is from Persaran/Mirta paper
foi_center  = round(foi_center*10)/10;          %convert to float?
delta_freq_true = (n_taper_all+1)./delta_time; % total bandwidth around; this is what fieldtrip wants, and then will redo the above computations under the hood

% Some intuitions to check on from Arjen:
% always need an odd number of tapers
% need at least 3 tapres to work well


cfg_hfa = [];
cfg_hfa.output       = 'pow';
cfg_hfa.channel      = 'all';
cfg_hfa.method       = 'mtmconvol';
cfg_hfa.taper        = 'dpss';
cfg_hfa.tapsmofrq    = delta_freq_true./2;                  %ft wants half bandwidth around the foi
cfg_hfa.keeptapers   = 'no';
cfg_hfa.pad          = 'maxperlen';                         %add time on either side of window
cfg_hfa.padtype      = 'zero';
cfg_hfa.foi          = foi_center;                          % analysis 2 to 30 Hz in steps of 2 Hz 
cfg_hfa.t_ftimwin    = ones(length(cfg_hfa.foi),1).*delta_time;    % length of time window; 0.5 sec, could be n_cycles./foi for n_cylces per win
cfg_hfa.toi          = 'all';%-buff_lim(1):0.1:1.5;         % time window centers
cfg_hfa.keeptrials   = 'yes';                               % must be 'yes' for stats

