%% ${SBJ} Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR52';
SBJ_vars.raw_file = '2017010513_0022.besa';
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

SBJ_vars.ch_lab.probes = {'RAM','RHH','RTH','RAC','ROF','RSMA','RAIN','RPIN','RBT','LAM','LHH','LTH','LAC','LOF'};
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'RAC*','ROF*','RSMA*','RAIN*','RPIN*','RAM5','RAM6','RAM7','RAM8'};%RAM elecs in vaINS

SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'LAM2','LAM3','LAM4','LHH2','LHH3',...%epileptic accoridng to Jack
    'RHH1','RHH2','RHH6','RAM3','RTH4','RTH5','LHH4','LTH3',...%epileptic according to Bob
    'RAM1','RHH10','RTH1','RTH2','RAC10','LAM10','RAIN10'...%out of brain
    };
SBJ_vars.ch_lab.eeg = {'FZ','CZ','OZ','C3','C4','O1','O2','T3','T4','F3','F4'};
SBJ_vars.ch_lab.eog = {'LUE','LLE','RUE','RLE'};
SBJ_vars.ch_lab.photod = {'DC02'};
SBJ_vars.ch_lab.mic    = {'DC03'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% First event is at 6 seconds in! Very tight, be careful with filtering...
% Last event is over by 495s
SBJ_vars.analysis_time = {[0.001 505]};

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
