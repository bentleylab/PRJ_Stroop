%% IR31 Processing Variables
[root_dir, ft_dir] = fn_get_root_dir();
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR31';
SBJ_vars.raw_file = '2015111117_0005.besa';
SBJ_vars.block_prefix = '';

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Stroop/data/' SBJ_vars.SBJ '/'];
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
%hdr = ft_read_header(SBJ_vars.dirs.raw_filename);
%SBJ_vars.orig_n_ch = length(hdr.label);
%SBJ_vars.orig_n_samples = hdr.nSamples;
%SBJ_vars.orig_srate = hdr.Fs;
%clear hdr;

SBJ_vars.ch_lab.probes = {'ROF','RAM','RHH','RTH','RAC','LAC','LOF','LIN','LAM','LHH'};
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'RAC*','ROF*','LAC*','LOF*','LIN*'};

SBJ_vars.ch_lab.bad = {...
    'RHH7','RTH5','RTH6','RTH7',...%epileptic
    'RAC3',...% Noisy
    'RAC8','RAC9','RAC10',...%out of brain;  'RAC7' is half out, but keeping it because RAC6-7 PSD looks good
    'RAM10','RHH1','RHH9','RHH10','RTH8','RTH9','RTH10','ROF7','ROF8','ROF9','ROF10',...%out of brain
    'LIN8','LIN9','LIN10','LIN11','LIN12','LHH7','LHH8','LHH9','LHH10',...%out of brain
    'LAC6','LAC7','LAC8','LAC9','LAC10','LOF8','LOF9','LOF10','LAM7','LAM8','LAM9','LAM10',...%out of brain
    'DC03','DC04'....% Not real data
    'E','LSh ','LLE','RSh','V1','V2','V3','V4','V5','V6','REF',...% Not real data
    'EKG*',...
    };
SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.eeg = {};
SBJ_vars.ch_lab.eog = {};
SBJ_vars.ch_lab.photod = {'DC01'};
SBJ_vars.ch_lab.mic    = {'DC02'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% Big bumps at 60, 100, (small)120, 180, 200, 240, 300 Hz, some other small ones
% LAC4 has HUGE noise! evrything is way bigger, new peak at 150
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {[0 910]};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 7;
SBJ_vars.artifact_params.hard_threshold_raw = 500;

SBJ_vars.artifact_params.std_limit_diff = 7;
SBJ_vars.artifact_params.hard_threshold_diff = 25;

%--------------------------------------
% Trials to Reject
%--------------------------------------
SBJ_vars.trial_reject_n = [112, 247, 275, 302, 303];
% OLD ix SBJ_vars.trial_reject_ix = [71, 198, 225, 251, 252];
