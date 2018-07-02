%% IR65 Processing Variables
[root_dir, ft_dir] = fn_get_root_dir();
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR65';
SBJ_vars.raw_file = '2017120714_0015.besa';
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

SBJ_vars.ch_lab.probes = {'RAM','RHH','RTH','RAC','ROF','LAM','LHH','LTH','LAC','LOF'};
SBJ_vars.ref_types     = {'BP','BP','BP','BP','BP','BP','BP','BP','BP','BP'};
SBJ_vars.ch_lab.ROI    = {'RAC*','ROF*','LAC*','LOF*',...
                            '-RAC3-4','-ROF1-2'};% rejected because of epileptic spread in preproc
                            % BEWARE: RAC6-7, LAC1-2, and some other channels like ROF2-3 have small spread...
SBJ_vars.ch_lab.eeg_ROI = {'CZ'};

SBJ_vars.ch_lab.prefix = '';    % before every channel except 'EDF Annotations'
SBJ_vars.ch_lab.suffix = '';    % after every channel except 'EDF Annotations'
SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'RHH1','RHH2','RHH3','RTH1','RTH2','RTH3','RTH4','RAM1','RAM2',...% epileptic
    'LHH4','LHH5','LAM1','LAM2','LAM3','LAM4',...% noisy, big ERP-like jumps
    'AST*','PST*','AMG*','HH*','TH*',...% Not real data
    'C3-1','C4-1','CZ-1','OZ-1','FPZ','Z','LUE-1','LLE-1','RUE-1','RLE-1','EKG-1',...% Not real data
    'REF','EKG','GND','XREF'...% note real data
    };
SBJ_vars.ch_lab.eeg = {'FZ','CZ','C3','C4','OZ'};% copies of each of these plus FPZ and Z...
% SBJ_vars.ch_lab.CZ_lap_ref = {};
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
SBJ_vars.analysis_time = {[145 490], [610 1290]};

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
% SBJ_vars.trial_reject_n = [];
