%% IR48 Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR48';
SBJ_vars.raw_file = '2016092915_0028.besa';
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

SBJ_vars.ch_lab.probes = {'FG','AF','MF'};
SBJ_vars.ref_types     = {'CAR','CAR','CAR'};
SBJ_vars.ch_lab.ROI    = {'MF*','FG*','AF*'};
SBJ_vars.ch_lab.eeg_ROI = {'FPZ'};

SBJ_vars.ch_lab.bad = {...
    'FG27','FG28','FG29','FG37','FG38','FG45','FG46','FG47','FG48',...%epileptic, Jack (via Jie) said FG47+48
    'FG43',...%slow ~1 Hz rhythm...
    'AF12',...%noisy
    'E','REF','V1','V2','V3','V4','V5','V6','C10','OFG*','PST*','LAM*','LHH*', ...% Not real data
    'EKG'...
    };
%     'RSH','LSH','LLE',...% EOG?
%     'FG12',...% slow rhythm like FG43 but lower amplitude
%     'FG25',...% lots of drift, can maybe be saved with preprocessing?
SBJ_vars.ref_exclude = {'FG44'}; % sometimes reflects artifact in FG45/46
SBJ_vars.ch_lab.eeg = {'FPZ','FP1','FP2','OZ','T5','T6'};
SBJ_vars.ch_lab.FPZ_lap_ref = {'FP1','FP2'};
% SBJ_vars.ch_lab.eeg_bad = {};
SBJ_vars.ch_lab.eog = {'RSH','LSH','LLE'}; % Janna says these channels are empty?
SBJ_vars.ch_lab.photod = {'DC02'};
SBJ_vars.ch_lab.mic    = {'DC04'};%burried in noise, so can't visually see responses in plot (can hear just barely)

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% task stats ~199.5s, ends ~1266s and goes to ~1280s
SBJ_vars.analysis_time = {[189 1276]};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
% SBJ_vars.artifact_params.std_limit_raw = 7;
% SBJ_vars.artifact_params.hard_threshold_raw = 1000;

% SBJ_vars.artifact_params.std_limit_diff = 7;
% SBJ_vars.artifact_params.hard_threshold_diff = 100;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% These should be indices AFTER SBJ05 has run!
% SBJ_vars.trial_reject_ix = [];
