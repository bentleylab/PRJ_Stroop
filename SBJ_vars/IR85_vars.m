%% IR85 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip')); addpath(ft_dir); ft_defaults; end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR85';
SBJ_vars.raw_file = {'2019011311_0001.besa'};
SBJ_vars.block_name = {''};
SBJ_vars.low_srate  = [0];

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Stroop/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.import  = [SBJ_vars.dirs.SBJ '01_import/'];
SBJ_vars.dirs.preproc = [SBJ_vars.dirs.SBJ '02_preproc/'];
SBJ_vars.dirs.events  = [SBJ_vars.dirs.SBJ '03_events/'];
SBJ_vars.dirs.proc    = [SBJ_vars.dirs.SBJ '04_proc/'];
SBJ_vars.dirs.recon   = [SBJ_vars.dirs.SBJ '05_recon/'];
dirs_fields = fieldnames(SBJ_vars.dirs);
for field_ix = 1:numel(dirs_fields)
    if ~strcmp(dirs_fields{field_ix},'raw_filename') && ~exist(SBJ_vars.dirs.(dirs_fields{field_ix}),'dir')
        mkdir(SBJ_vars.dirs.(dirs_fields{field_ix}));
    end
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
SBJ_vars.recon.elec_mni_s = [];%[SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_s.mat'];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'LAM','LHH','LTH',...
                              'RAM','RHH','RTH','RAC','ROF','RPC','RMC','RIN','RIT'};
SBJ_vars.ch_lab.probe_type = {'seeg','seeg','seeg',...
                              'seeg','seeg','seeg','seeg','seeg','seeg','seeg','seeg','seeg'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP',...
                              'BP','BP','BP','BP','BP','BP','BP','BP','BP'};
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]); error('probes ~= type+ref');end;
SBJ_vars.ch_lab.nlx        = [0,0,0,0,0,0,0,0,0,0,0,0];
SBJ_vars.ch_lab.ROI        = {'all',...
                              '-RIN9-10','-RPC6-7'}; % excluded due to large artifacts
SBJ_vars.ch_lab.eeg_ROI    = {};

%SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
% RIT2-5 have some spike-like waves too that spread to RHH, RAM, RTH when large...
SBJ_vars.ch_lab.bad = {...
    'LHH8',... % epileptic
    'LAM7','LAM8','LTH8','LTH9','LTH10',... % epileptic spread; LTH10 maybe out too?
    'LTH5',... % loose channel
    'RIT10',... % HF noise, other RIT to some degree
    'Z',... % extra EEG channel?
    'EKG',... % EKG
    'Mark1','Mark2','DC03','DC04','E','REF','Events'... % not real data
    };
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1,2,2,2,2,2,1,2,0,0,0,0,0,0,0,0,0];
if numel(SBJ_vars.ch_lab.bad)~=numel(SBJ_vars.ch_lab.bad_code);error('bad ~= bad_code');end
SBJ_vars.ch_lab.eeg = {'FZ','CZ','OZ','C3','C4'};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
SBJ_vars.ch_lab.eog = {'LUC','LLC','RUC','RLC'};
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
SBJ_vars.analysis_time = {{[333 1814]}};

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
SBJ_vars.trial_reject_n = [16 176 275];
