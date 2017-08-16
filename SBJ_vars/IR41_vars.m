%% IR41 Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR41';
SBJ_vars.raw_file = '2016051011_0001.besa';
SBJ_vars.block_prefix = '';

SBJ_vars.dirs.SBJ     = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
if ~exist(SBJ_vars.dirs.import,'dir')
    mkdir(SBJ_vars.dirs.import);
end
if ~exist(SBJ_vars.dirs.preproc,'dir')
    mkdir(SBJ_vars.dirs.preproc);
end
if ~exist(SBJ_vars.dirs.events,'dir')
    mkdir(SBJ_vars.dirs.events);
end
if ~exist(SBJ_vars.dirs.proc,'dir')
    mkdir(SBJ_vars.dirs.proc);
end

SBJ_vars.dirs.raw_filename = strcat(SBJ_vars.dirs.raw,SBJ_vars.raw_file);

%--------------------------------------
% Channel Selection
%--------------------------------------
hdr = ft_read_header(SBJ_vars.dirs.raw_filename);
SBJ_vars.orig_n_ch = length(hdr.label);
SBJ_vars.orig_n_samples = hdr.nSamples;
SBJ_vars.orig_srate = hdr.Fs;
clear hdr;

SBJ_vars.ch_lab.probes = {'LHH','LAC','LES','LOF','LSM','LPC','LIN',...
                          'RAC','RAM','ROF','RIN','RMT','RSM'};
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'LAC*','LES*','LOF*','LSM*','LPC*','LIN*',...
                          'RAC*','ROF*','RIN*','RMT*','RSM*'};

SBJ_vars.ch_lab.bad = {...
    'RAM1','RAM2','RAM3','LHH1','LHH2','LHH3',...%epileptic
    'RMT1','LPC7','LPC8','LPC9',...%noisy/bad, also watch out for spread from LPC7-9 in LPC6,11
    'LES1',...%removing because huge variance in LES1-2!
    'RAM9','RAM10','LHH10','LOF9','LOF10',...%out of brain
    'ROF11','ROF12','ROF13','ROF14','ROF15','LAC8','LAC9','LAC10',...%out of brain
    'RAC8','RAC9','RAC10','RAC11','RAC12','RMT9','RMT10',...%out of brain
    'LPC12','LPC13','LPC14','LPC15','LES8','LES9','LES10',...%out of brain
    'LIN13','LIN14','LIN15','LSM10','RSM9','RSM10','RIN12',...%out of brain
    'DC03','DC04','EYE1','EYE2','EYE3','EYE4',....% Not real data
    'E','LSH','LLE','RSH','V1','V2','V3','V4','V5','V6',...% Not real data
    'EKG*'...
    };
% edge of cortex: LOF1,ROF1,RMT8,LIN12,RSM8
%     'LHH4',...%noisy in Kata's data
% SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.eeg = {'F3','FZ','F4','T3','T4','CZ','O1','O2','OZ','C3','C4'};
SBJ_vars.ch_lab.photod = {'DC01'};
SBJ_vars.ch_lab.mic    = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% Looks like there is noise in 25-30 (???!), nothing at 60, then again at 180, 300 Hz...
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% events start ~155 or 160s to ~384; ~580 to end (~1360?)
SBJ_vars.analysis_time = {[140 400], [560 1360]};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 8;
SBJ_vars.artifact_params.hard_threshold_raw = 1000; % ~1000 and below is the main cloud of points

SBJ_vars.artifact_params.std_limit_diff = 8;
SBJ_vars.artifact_params.hard_threshold_diff = 100;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% These should be indices AFTER SBJ05 has run!
SBJ_vars.trial_reject_ix = [27 51 63 67 129 148 151 193 195];
