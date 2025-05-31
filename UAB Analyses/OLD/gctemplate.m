cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, roi_trl);

%%
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'hanning';
cfg.output    = 'fourier';
% cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, roi_trl);
%%
cfg                     = [];
cfg.method              = 'granger';
cfg.granger.sfmethod    = 'bivariate';
cfg.granger.conditional = 'no';
cfg.granger.channelcmb  = fn_create_channel_combinations(MPFC,LPFC);
granger       = ft_connectivityanalysis(cfg, freq);

cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 0.1];
cfg.xlim      = [0 40];
ft_connectivityplot(cfg, granger);

%% Energy in tapers

nTrials = length(freq.cumtapcnt); % Number of trials (605 in your case)
nTapers = 3; % Number of tapers per trial
nChannels = size(freq.fourierspctrm, 2); % Number of channels (7 in your case)
nFrequencies = size(freq.fourierspctrm, 3); % Number of frequencies (252 in your case)

fourierspctrm_reshaped = reshape(freq.fourierspctrm, [nTrials, nTapers, nChannels, nFrequencies]);

spectral_energy = abs(fourierspctrm_reshaped).^2;

energy_per_taper_channel = squeeze(sum(spectral_energy, 4)); % Size: [nTrials, nTapers, nChannels]

total_energy_per_trial_channel = sum(energy_per_taper_channel, 2); % Sum across tapers; Size: [nTrials, nChannels]

percentage_energy_per_taper_channel = (energy_per_taper_channel ./ total_energy_per_trial_channel) * 100;

valid_tapers_channel = percentage_energy_per_taper_channel >= 90;

