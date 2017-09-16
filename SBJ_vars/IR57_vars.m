%% IR57 Processing Variables
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR57';
SBJ_vars.raw_file = '2017032319_0023.besa';
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

SBJ_vars.ch_lab.probes = {'RSM','RAC','ROF','RIN','RTI','RAM','RHH','RTH',...
                         'LSMA','LAC','LOF','LIN','LTI','LAM','LTH'};%'LHH' doesn't count because all elecs are bad
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP',...
                          'BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'RSM*','RAC*','ROF*','RIN*','RTI*',...
                              'LAM4-5','LAC*','LOF*'}; %LAM4 is inferior anterior insula

SBJ_vars.ch_lab.bad = {...
    'RHH1','RHH2','RHH3','RHH4','ROF1','RAM1','RTH1','RTH2',...%epileptic
    'LHH1','LHH2','LHH3','LHH4','LHH5','LHH6','LHH7','LHH8','LHH9','LHH10',...%epileptic
    'LTH1','LTH2','LTH3','LAM1','LAM2','LAM3',...%epileptic
    'RSM2','RSM3','RIN9','RIN10','LTI2','LTI3','RHH6','LHH1','LHH2','LHH4','LHH6',...%bad line noise
    'RSM8','RSM9','RAC10','RAM10','RTH9','RTH10','LOF10','LAM10',...%out of brain
    'RHH7',...%marked as bad in visual for some reason?
    'DC01','DC03',...% not real data
    'LUC','LLC','RUC','RLC','XREF','E',...% not real data
    'EKG'...
    };
% SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.eeg = {'FPZ' 'CZ' 'OZ' 'C3' 'C4' 'Z' 'FP1' 'FP2' 'T3' 'T4' 'O1' 'O2'};
SBJ_vars.ch_lab.photod = {'DC02'};
SBJ_vars.ch_lab.mic    = {'DC04'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% most have only regular harmonics with normal width, and 120 and 240 are weak
% RBT, RPIN, RSMA,LAC have an extra peak at 200
% RHH6 has really bad at all harmonics, like LUE and other nonsense chan
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
SBJ_vars.analysis_time = {[70 1040]};

% %--------------------------------------
% % Save SBJ_vars
% %--------------------------------------
% out_filename = [SBJ_vars.dirs.SBJ SBJ '_vars.mat'];
% save(out_filename,'SBJ_vars');

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
% Add rejection of LOF1-2, maybe RTI2-3, RIN4-5
SBJ_vars.artifact_params.std_limit_raw = 7;
SBJ_vars.artifact_params.hard_threshold_raw = 700;

SBJ_vars.artifact_params.std_limit_diff = 7;
SBJ_vars.artifact_params.hard_threshold_diff = 35;

%--------------------------------------
% Trials to Reject
%--------------------------------------
% These should be indices AFTER SBJ05 has run!
SBJ_vars.trial_reject_n = [50 125 131 134 213 214 306 311];
% BADBADBAD SBJ_vars.trial_reject_ix = [22, 95, 174, 176, 180, 245, 254];
