%% IR32 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR32';
SBJ_vars.raw_file = {'2015121512_0002.edf'};
SBJ_vars.block_name = {''};
SBJ_vars.low_srate  = [0];

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Stroop/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
SBJ_vars.dirs.recon   = [SBJ_vars.dirs.SBJ '05_recon/'];
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
if ~exist(SBJ_vars.dirs.recon,'dir')
    mkdir(SBJ_vars.dirs.recon);
end

SBJ_vars.dirs.raw_filename = strcat(SBJ_vars.dirs.raw,SBJ_vars.raw_file);

SBJ_vars.recon.surf_l     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_lh.mat'];
SBJ_vars.recon.surf_r     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_cortex_rh.mat'];
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_fr.mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_frv.mat'];
SBJ_vars.recon.elec_mni_s = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_fsavg_frs.mat'];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'FPG','IHL','IHR','AG'};
SBJ_vars.ch_lab.probe_type = {'ecog','ecog','ecog','ecog'};
SBJ_vars.ch_lab.ref_type   = {'CAR','CAR','CAR','CAR'};
SBJ_vars.ch_lab.ROI        = {'IHL*','IHR*','AG*','FPG*'};
SBJ_vars.ch_lab.eeg_ROI    = {};

SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ch_lab.ref_exclude = {'FPG11','FPG12','FPG19','FPG30','FPG34','IHR20','IHR21','IHR28'}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'FPG27','FPG28','FPG29','FPG35','FPG36','FPG37',...%epileptic
    'FPG7','FPG25','FPG40',...%bad/noisy
    'AG5',...% bad channel, mirrors AG6 and negative corr with all others
    'DC03','DC04'....% Not real data
    'E','LSH','LLE','RSH','V1','V2','V3','V4','V5','V6','REF',...% Not real data
    'EDF Annotations','---(1)','---(2)','---(3)','---(4)','---(5)','---(6)',... % Not real data
    'EKG*'
    };
% shifted to ref_exclude:
%     'FPG19','FPG34',...% epileptic spread
%     'IHR20','IHR21','IHR30','IHR31',...% noisy
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1 1 1 1 1 1 2 2 2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
if numel(SBJ_vars.ch_lab.bad)~=numel(SBJ_vars.ch_lab.bad_code);error('bad ~= bad_code');end
SBJ_vars.ch_lab.eeg = {};
SBJ_vars.ch_lab.eog = {};
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
SBJ_vars.analysis_time = {{[62 1170]}};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 9;
SBJ_vars.artifact_params.hard_threshold_raw = 1000;%without FPG: 600;

SBJ_vars.artifact_params.std_limit_diff = 9;
SBJ_vars.artifact_params.hard_threshold_diff = 150;

%--------------------------------------
% Trials to Reject
%--------------------------------------
SBJ_vars.trial_reject_n = [193 244 251 307];
% BAD, and naive: SBJ_vars.trial_reject_ix = [146, 175, 181, 207, 209, 228, 232, 237, 242];
