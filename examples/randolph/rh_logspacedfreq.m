function freq = rh_logspacedfreq(data, trialnb, toi)

% might need some tweaking, is optimized for the ECoG context experiment

%% freq analysis in log spacing

% prep log spacing for freq analysis on the data
% toi = -4.5:0.1:1;
foi_center = 2.^[0:1/4:8]; % Center frequencies
octave = 3/4;              % Frequency resolution
foi_min = 2^(-octave/2)*foi_center;
foi_max = 2^(octave/2)*foi_center;
foi = (foi_min+foi_max)/2;
delta_freq = foi_max-foi_min;
delta_time = 0.5;
n_taper_all = max(1,round(delta_freq.*delta_time-1))
foi_center  = round(foi_center*10)/10;
delta_freq_true = (n_taper_all+1)./delta_time;

% compute freq analysis
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.toi     = toi;
cfg.foi     = foi_center;
cfg.tapsmofrq = delta_freq_true./2; 
cfg.keeptrials = 'yes';
cfg.keeptapers = 'no';
cfg.t_ftimwin = ones(1,length(cfg.tapsmofrq))*delta_time;
cfg.taper = 'dpss';
cfg.pad = 'maxperlen';

if strcmp(trialnb, 'all') ==1
    cfg.trials = 'all';
else
    cfg.trials = find(data.trialinfo(:,1) == trialnb);
end
freq = ft_freqanalysis(cfg,data);

end