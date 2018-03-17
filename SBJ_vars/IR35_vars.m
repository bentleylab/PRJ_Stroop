%% IR35 Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR35';
SBJ_vars.raw_file = '2016021711_0006.besa';
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

SBJ_vars.ch_lab.probes = {'RAM','RHH','RTH','RIN','ROF','LAM','LHH','LTH','LAC','LOF','LIN','LPC'};
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'RIN*','ROF*','LAC*','LOF*','LIN*','LPC*','-LPC6-7'};

SBJ_vars.ch_lab.bad = {...
    'LHH1','LHH2','LHH3','RHH1','RHH2','RHH3','LTH3','LTH4',...%epileptic
    'LAC4',...%crazy strong noise
    'RIN10','LAM8','LAM9','LAM10','LPC9','LPC10','LHH9','LHH10',...%out of brain
    'LTH1','LTH10','ROF8','ROF9','ROF10','RTH8','RTH9','RTH10',...%out of brain
    'RHH7','RHH8','RHH9','RHH10','RAM9','RAM10','LIN10',...%out of brain
    'LOF1','LOF8','LOF9','LOF10','LAC11','LAC12',...%out of brain
    'NULL','NULL-1','NILL','NULL-2','DC03','DC04'....% Not real data
    'E','LSH','LLE','RSH','V1','V2','V3','V4','V5','V6','xREF',...% Not real data
    'EKG*'...
    };
SBJ_vars.ref_exclude = {}; %exclude from the CAR
% SBJ_vars.ch_lab.eeg = {};
% SBJ_vars.ch_lab.eog = {};
SBJ_vars.ch_lab.photod = {'DC01'};
SBJ_vars.ch_lab.mic    = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% Big bumps at 60, 100, (small)120, 180, 200, 240, 300 Hz, some other small ones
% LAC4 has HUGE noise! evrything is way bigger, new peak at 150
% most have only regular harmonics with normal width, and 120 and 240 are weak
SBJ_vars.notch_freqs = [60 120 180 240 300]; %100, 200 should come out in re-referencing
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {[105 1140]};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 7;
SBJ_vars.artifact_params.hard_threshold_raw = 400;

SBJ_vars.artifact_params.std_limit_diff = 7;
SBJ_vars.artifact_params.hard_threshold_diff = 40;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% These should be indices AFTER SBJ05 has run!
SBJ_vars.trial_reject_n = [42];
% old index version: SBJ_vars.trial_reject_ix = [42];
