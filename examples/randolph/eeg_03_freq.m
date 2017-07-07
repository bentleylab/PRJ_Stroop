%% freq analysis
clc
clear all
close all

pathname = '/home/knight/rhelfric/Data/EEG_10Hzcontext/';
cd(pathname)

subjects = {'S01', 'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S09', 'S10', 'S11', ...
            'S12', 'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S20'};

for s = 1:length(subjects)
    
    loadname = [pathname, subjects{s}, '_ica'];
	load(loadname)

    % calculate total power and csd for connectivity analyses
    cfg            = [];
    cfg.method     = 'mtmconvol';
    cfg.output     = 'pow';
    cfg.toi        = -0.5:0.05:4;
    cfg.foi        = 1:1:30;
    cfg.keeptrials = 'yes';
    cfg.t_ftimwin  = ones(1,length(cfg.foi))* 0.5;
    cfg.taper      = 'hanning';
    cfg.pad        = 'maxperlen';    
    cfg.trials     = 'all';
    
    freq = ft_freqanalysis(cfg,data);
    
    % z correct
    freq = rh_zbaseline(freq, -0.2, 0);
    freq.powspctrm = freq.zspctrm; 
    freq.zspctrm = [];
    
    savename = [pathname, subjects{s}, '_freq'];
    save(savename, 'freq', '-v7.3')     
    
end