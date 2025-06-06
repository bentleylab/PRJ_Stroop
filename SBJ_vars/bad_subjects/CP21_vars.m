%% IR21 Processing Variables
error('Do not processing this patient. No microphone data so cannot proceed.');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'CP21';
SBJ_vars.raw_file = 'CP21_raw_ft.mat';
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
% load(SBJ_vars.dirs.raw_filename);
SBJ_vars.orig_n_ch = 74;
SBJ_vars.orig_n_samples = 1155000;
SBJ_vars.orig_srate = 1000;
% clear hdr;

SBJ_vars.ch_lab.probes = {'LIHA','LIHP','LLFP','LLT','LMT','LOF','RIHA','RIHP','RLFP','RLT','RMT','ROF'};
SBJ_vars.ref_types     = {'CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR','CAR'}; %consdire BP!
SBJ_vars.ch_lab.ROI    = {'LIHA*','LIHP*','LLFP*','LOF*','RIHA*','RIHP*','RLFP*','ROF*'};

SBJ_vars.ch_lab.bad = {...
    'FG27','FG28','FG29','FG37','FG38','FG45','FG46','FG47','FG48',...%epileptic, Jack (via Jie) said FG47+48
    'FG43',...%slow ~1 Hz rhythm...
    'AF12',...%noisy
    'E','REF','V1','V2','V3','V4','V5','V6','C10','OFG*','PST*','LAM*','LHH*', ...% Not real data
    'RSH','LSH','LLE',...% EOG?
    'EKG'...
    };
%     'FG12',...% slow rhythm like FG43 but lower amplitude
%     'FG25',...% lots of drift, can maybe be saved with preprocessing?
SBJ_vars.ch_lab.ref_exclude = {'FG44'}; % sometimes reflects artifact in FG45/46

SBJ_vars.ch_lab.eeg = {'FPZ','FP1','FP2','OZ','T5','T6'};
SBJ_vars.ch_lab.eog = {};
SBJ_vars.ch_lab.photod = {'DC02'};
SBJ_vars.ch_lab.mic    = {'DC04'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {'error'};

% %--------------------------------------
% % Save SBJ_vars
% %--------------------------------------
% out_filename = [SBJ_vars.dirs.SBJ SBJ '_vars.mat'];
% save(out_filename,'SBJ_vars');

