%% ${SBJ} Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR37';
SBJ_vars.raw_file = '2016031117_0051.besa';
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

SBJ_vars.ch_lab.probes = {'FG','OF','PT','AT','ID'};    %'CD' all bad
SBJ_vars.ref_types     = {'CAR','CAR','CAR','CAR','BP'};
SBJ_vars.ch_lab.ROI    = {'FG','OF','ID'};

SBJ_vars.ch_lab.mislabel = {{'FG17-1','FG27'}};
SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'AT7','AT1','AT2',...%Jack epileptic
    'PT7','PT8','PT9','PT10',...% Bob epileptic
    'OF1','FG4','FG5','FG6','FG7','FG8','FG15','FG16',...%bad channels
    'CD1','CD2','CD3','CD4','CD5','CD6',...%crazy line noise, but trying to save becuase aMCC!
    'C4','EKG','E','NULL','NULL-1','NULL-2','NULL-3','NULL-4','XREF','DC03','DC04'...%junk
    };
% SBJ_vars.ch_lab.eeg = {};
% SBJ_vars.ch_lab.eog = {};
SBJ_vars.ch_lab.photod = {'DC01'};
SBJ_vars.ch_lab.mic    = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% photodiode wasn't on for most of 1st block
% kicks on at 190.5, which is B1T26 (10 full trials left)
% last block over ~1182
SBJ_vars.analysis_time = {[180 1194]};

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
