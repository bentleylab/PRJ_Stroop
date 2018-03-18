%% IR21 Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR21';
SBJ_vars.raw_file = '2015063014_0002.edf';
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

SBJ_vars.ch_lab.probes = {'ROF','RAH','RC','LIN','LAH'};
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'ROF*','RC*','LIN*'};

SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'

SBJ_vars.ch_lab.bad = {...
    'RAH6','RAH7',...%epileptic
    'RC13','RC14','RC15','RC16',...%out of brain
    'LIN16',...% strongest HF nois bursts; causes excessive variance rejection in derivative data
    'DC02','DC03'....% Not real data
    'E','LSh ','LLE','RSh','V1','V2','V3','V4','V5','V6','REF',...% Not real data
    'EKG','C16','EDF Annotations'...
    };
SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.eeg = {};
SBJ_vars.ch_lab.eog = {};
SBJ_vars.ch_lab.photod = {'DC01'};
SBJ_vars.ch_lab.mic    = {'DC04'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% nothing at 60, most have samll 300 bump and maybe 180
% LIN has 100, 120, 180, 200, 240, 300 Hz, even 360
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
%   12 blocks, if you count the one split by the screen freeze as two blocks
%   events start 75, screen freezes ~168, kicks back in ~213, over !1365
%   not sure about excluding the freeze, which is ~45s
SBJ_vars.analysis_time = {[65 1375]};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 5;
SBJ_vars.artifact_params.hard_threshold_raw = 200;

SBJ_vars.artifact_params.std_limit_diff = 5;
SBJ_vars.artifact_params.hard_threshold_diff = 20;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% These should be indices AFTER SBJ05 has run!
SBJ_vars.trial_reject_ix = [];
