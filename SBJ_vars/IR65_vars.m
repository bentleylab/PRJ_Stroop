%% IR65 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR65';
SBJ_vars.raw_file = {'2017120714_0015.besa'};
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
SBJ_vars.recon.wm_l       = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_wm_lh.mat'];
SBJ_vars.recon.wm_r       = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_wm_rh.mat'];
SBJ_vars.recon.infl_l     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_infl_lh.mat'];
SBJ_vars.recon.infl_r     = [SBJ_vars.dirs.recon 'Surfaces/' SBJ_vars.SBJ '_infl_rh.mat'];
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_f.mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_frv.mat'];
SBJ_vars.recon.elec_mni_s = [];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'RHH','RTH','RAC','ROF','LAM','LHH','LTH','LAC','LOF'};% RAM gone (all bad)
SBJ_vars.ch_lab.probe_type = {'seeg','seeg','seeg','seeg','seeg','seeg','seeg','seeg','seeg'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP','BP','BP','BP'};
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]); error('probes ~= type+ref');end;
SBJ_vars.ch_lab.ROI        = {'all',... %'RAC*','ROF*','LAC*','LOF*',...
                              '-RAC3-4'};% rejected because of epileptic spread in preproc
                              % BEWARE: RAC6-7, LAC1-2, and some other channels like ROF2-3 have small spread...
SBJ_vars.ch_lab.eeg_ROI    = {'CZ'};

%SBJ_vars.ch_lab.prefix = '';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'RHH1','RHH2','RHH3','RTH1','RTH2','RTH3','RTH4','RAM1','RAM2',...% epileptic
    'LHH4','LHH5','LAM1','LAM2','LAM3','LAM4',...% noisy, big ERP-like jumps
    'RAM3','RAM4','RAM5','RAM6','RAM7','RAM8','RAM9','RAM10','RTH5','RTH6','RHH4','RHH5',...%epileptic spread
    'ROF1','ROF2',...%epileptic spread
    'LAM10','RHH10',...% out of brain
    'AST*','PST*','AMG*','HH*','TH*',...% Not real data
    'C3-1','C4-1','CZ-1','OZ-1','FPZ','Z','LUE-1','LLE-1','RUE-1','RLE-1','EKG-1',...% Not real data
    'DC01','DC04','REF','EKG','GND','XREF'...% note real data
    };
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
if numel(SBJ_vars.ch_lab.bad)~=numel(SBJ_vars.ch_lab.bad_code);error('bad ~= bad_code');end
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
SBJ_vars.analysis_time = {{[145 490], [610 1290]}};

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
SBJ_vars.trial_reject_n = [];
