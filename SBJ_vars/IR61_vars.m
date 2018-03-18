%% IR61 Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR61';
SBJ_vars.raw_file = '2009010101_0002.besa';
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

SBJ_vars.ch_lab.probes = {'LAM','LHH','LTH','LOF','LAC','RAM','RHH','RTH','ROF','RAC'};
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'LOF*','LAC*','ROF*','RAC*','RAM4','RAM5','RAM6',... % RAM5,6 inf. ant. Insula, RAM4 is WM nearby
                            'RHH5','RHH6','RHH7'}; % RHH5,6 in inf. post. Insula, RHH7 WM nearby
SBJ_vars.ch_lab.eeg_ROI = {'CZ','FZ','FPZ'};

SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'LAM1','LAM2','LHH1','LHH2','LHH3','LTH1','LTH2','RHH3','RHH4','RAM2','RTH2','RTH3',...%epileptic
    'LOF10','LAC10','RAM10','RHH1','RAC10','ROF10',...%out of brain
    'LAM3','LHH7',...%loose
    'Z','---(13)','---(14)','---(17)','---(18)',...% not real?
    'REF','EKG','DC01','DC04'...%non-neural
    };
    % BEWARE: LTH4 (ventricle), LHH1+RAM9 (border)
    % watch out for FZ, Janna said it was bad
SBJ_vars.ch_lab.eeg = {'FZ','FPZ','CZ','OZ','C3','C4'};
SBJ_vars.ch_lab.CZ_lap_ref = {'C3','C4'};
SBJ_vars.ch_lab.FZ_lap_ref = {'C3','C4'};
SBJ_vars.ch_lab.FPZ_lap_ref = {'C3','C4'};
SBJ_vars.ch_lab.eog = {'LLE','LUE','RLE','RUE'};
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
SBJ_vars.analysis_time = {[138 872], [1230 1720]};

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
