% OS4 :: Convert EDF, Setup Format, Filter/Resamp, Define Trials, Reference ::  
% Strips (3 x 6chans, 1 X 4chans) = 22
% XX, XX, XX


%% Load Data

clear all
addpath('/home/knight/julia/Scripts/Extraction/');
addpath('/home/knight/julia/fieldtrip-20151020');
addpath('/home/knight/julia/ECoG_DirCog');

ft_defaults

cd(['/home/knight/julia/ECoG_DirCog/OS4/raw_data/']);
[data_info,rawdata]=edfread('OS4_task.edf');


%% Setup Data

cd(['/home/knight/julia/ECoG_DirCog/OS4/analysis/task/']);

% ~~ enter relevant info 
data.trial{1,1} = rawdata;
data.label = data_info.label;
data.fsample = 1024;
data.nSamples = size(data.trial{1,1},2);
data.time{1,1}  = (1/data.fsample):(1/data.fsample):(data.nSamples/data.fsample);

% ~~ read into fieldtrip format
cfg = [];
cfg.dataset = 'data';
cfg.continuous = 'yes';
[data] = ft_preprocessing(cfg,data);

data.nTrials = 1;
data.nChans = dataInfo.nChans;
data.rawhdr.orig_srate = 1024;
data.subj = 'OS4';
data.taskdir = ('/home/knight/julia/ECoG_DirCog/OS4/analysis/task/');
data.eeg_nChans = 22;    % total number of EEG channels
data.numTrials = 728; 

save('data','data','-v7.3');


%% 1. Downsample

cfg = [];
cfg.resamplefs = 1000;
[data_resamp] = ft_resampledata(cfg, data);

clear data; 


%% 2. Filter
% Lowpass filter @ 160Hz
% Bandstop filter @ 60Hz, 120Hz, 180Hz

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 160;
cfg.lpfiltwintype = 'hann';
[data_filt] = ft_preprocessing(cfg,data_resamp);

cfg.bsfilter = 'yes'; 
cfg.bsfreq = [59 61];
cfg.bsfiltwintype = 'hann';
[data_filt] = ft_preprocessing(cfg,data_filt);

cfg = [];
cfg.bsfilter = 'yes'; 
cfg.bsfreq = [119 121];
cfg.bsfiltwintype = 'hann';
[data_filt] = ft_preprocessing(cfg,data_filt);

cfg = [];
cfg.bsfilter = 'yes'; 
cfg.bsfreq = [179 181];
cfg.bsfiltwintype = 'hann';
[data_filt] = ft_preprocessing(cfg,data_filt);

save('data_filt','data_filt','-v7.3')
clear data_resamp; 


%% 3. Define trials

cd(['/home/knight/julia/ECoG_DirCog/OS4/raw_data/']);

% ~~ trigger channel = 23; threshold set at 2000;
cfg            = [];
cfg.dataset    = 'OS4_task.edf';
cfg.trialfun   = 'identifyEvents_edf_OSLO'; % set chanindx=45 for trigger 1
cfg            = ft_definetrial(cfg);


%%%%%%%%%%% ECOG DATA SAMP RATE == TRIALINFO SAMP RATE == 1000HZ %%%%%%%%%%%

cd(['/home/knight/julia/ECoG_DirCog/OS4/analysis/task/']);

% ~~ extract trial information
events_timeInfo = cfg.trl; 
events_timeInfo = events_timeInfo/(hdr.Fs/1000);
events_timeInfo=round(events_timeInfo);

[stimInfo,trialInfo] = identify_trial_conditions(events_timeInfo,data_filt.fsample);
trialInfo = round(trialInfo);
data_filt.trialInfo = trialInfo;
save('trial_info.mat','events_timeInfo','stimInfo','trialInfo','eventlog');


% ~~ epoch data
cfg = [];
cfg.trl = trialInfo(:,1:3);
cfg.trialdef.pretone = 1; % in secs
cfg.trialdef.posttone = 2; % in secs
cfg.trialdef.preswitch = 1; % in secs
cfg.trialdef.postswitch = 4; % in secs
data_epoched = ft_redefinetrial(cfg,data_filt);

cfg = [];
cfg.detrend = 'yes';
cfg.demean = 'yes';
[data_epoched] = ft_preprocessing(cfg, data_epoched);
data_epoched.trialInfo = trialInfo;

clear data_filt

%% 4. Remove bad channels & trials 

% ~~ first determine which trials contained bad epochs
load('artifact_epochs.mat'); % saved in secs
artifact_epochs=artifact_epochs*data_epoched.fsample; % change units to data's samp rate

cfg                         = [];
cfg.artfctdef.reject        = 'complete';
cfg.artfctdef.eog.artifact  = artifact_epochs;
data_goodtrials             = ft_rejectartifact(cfg, data_epoched);

% ~~ identify bad & non-EEG channels
for i=1:numel(chanInfo.excludeChannels)
    excludeChannels{1,i} = strcat('-',chanInfo.excludeChannels{i}); 
end

% TO DO: Check output  
% 1) all irrelevant channels not included (eg: DCs, EKGs, REF, depth ref chans, other bad channels not already excluded)
% 2) bad trials already removed
cfg                     = [];
cfg.method              = 'trial';
cfg.channel             = {'all', excludeChannels{:}};
cfg.keepchannel         = 'no'; 
[data_cleaned]          = ft_rejectvisual(cfg,data_goodtrials);
data_cleaned.trialInfo  = data_epoched.trialInfo(data_cleaned.cfg.trials); %save trialInfo for "good" trials only

% ~~ re-preprocess: delete above bad trials using ft_select data
badtrials = [95 334:335 459];
goodtrials = setxor(badtrials,[1:length(data_cleaned.trial)]);

cfg             = [];
cfg.trials      = goodtrials;
cfg.keeptrial   = 'no';
[data_cleaned]  = ft_selectdata(cfg,data_cleaned);

save('data_cleaned','data_cleaned','events_timeInfo','stimInfo','trialInfo','-v7.3');

clear data_epoched data_goodtrials

%% 5. Common Avg Reference
% Each strip has different level of noise. Thus, CAR within each strip. 

strip_chans = [];
strip_chans = ft_channelselection({'VST*'},data_cleaned.label); 
cfg = [];
cfg.channel = strip_chans;
cfg.reref = 'yes';
cfg.refchannel = 'all';
[vst] = ft_preprocessing(cfg,data_cleaned);

strip_chans = [];
strip_chans = ft_channelselection({'HST*'},data_cleaned.label); 
cfg = [];
cfg.channel = strip_chans;
cfg.reref = 'yes';
cfg.refchannel = 'all';
[hst] = ft_preprocessing(cfg,data_cleaned);

strip_chans = [];
strip_chans = ft_channelselection({'HTP*'},data_cleaned.label); 
cfg = [];
cfg.channel = strip_chans;
cfg.reref = 'yes';
cfg.refchannel = 'all';
[htp] = ft_preprocessing(cfg,data_cleaned);

strip_chans = [];
strip_chans = ft_channelselection({'HTO*'},data_cleaned.label); 
cfg = [];
cfg.channel = strip_chans;
cfg.reref = 'yes';
cfg.refchannel = 'all';
[hto] = ft_preprocessing(cfg,data_cleaned);

[data_reref] = ft_appenddata([], vst, hst, htp, hto);
clear vst hst htp hto


%% 6. Final Visual Inspection & Save Data (WM ref)

% ~~ visually inspect the referenced data
cfg = [];
cfg.channel    = {'all'}; 
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data_reref);

% ~~ SAVE DATA (WM ref)
data_preproc = data_reref;
[sample,allTrial,cleanTrialNum]=intersect(data_preproc.sampleinfo(:,1),data_cleaned.cfg.previous.previous.previous.previous.previous.trl(:,1));
data_preproc.trialInfo = trialInfo(cleanTrialNum,:);

% ~~ Sanity Check - Examine Power Spectrum
cfg = [];
cfg.output      = 'pow';
cfg.method      = 'mtmfft';
cfg.taper       = 'hanning';
cfg.foi         = [1:1:200];
cfg.tapsmofrq   = 0.5*cfg.foi;

[fft_data] = ft_freqanalysis(cfg,data_preproc);

cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'log10';
[fft_data] = ft_math(cfg,fft_data);

figure;
channum = min(length(fft_data.label),56);
if channum <= 56
    for chan = 1:channum
        cfg = [];
        cfg.channel = chan;
        subplot(7,8,chan);
        ft_singleplotER(cfg,fft_data);
        title(fft_data.label{chan})
    end
else
    for chan = 1:56
        cfg = [];
        cfg.channel = chan;
        subplot(7,8,chan);
        ft_singleplotER(cfg,fft_data);
        title(fft_data.label{chan})
    end
    for chan = 56:channum
        cfg = [];
        cfg.channel = chan;
        subplot(7,8,chan-56);
        ft_singleplotER(cfg,fft_data);
        title(fft_data.label{chan})
    end
end

save('data_preproc','data_preproc','cleanTrialNum','-v7.3');
clear data_reref data_preproc

