%% IR82 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR82';
SBJ_vars.raw_file = {'IR82_stroop_raw_clinical.mat'};
SBJ_vars.block_name = {''};% there's a second block I don't want to process right now, so leaving blank here
SBJ_vars.low_srate  = [0];

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Stroop/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.nlx     = [SBJ_vars.dirs.raw 'SU_2018-09-28_15-17-47/'];
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
%hdr = ft_read_header(SBJ_vars.dirs.raw_filename);
%SBJ_vars.orig_n_ch = length(hdr.label);
%SBJ_vars.orig_n_samples = hdr.nSamples;
%SBJ_vars.orig_srate = hdr.Fs;
%clear hdr;

SBJ_vars.ch_lab.probes     = {'RAM','RHH','RTH','ROF','RAC','LHH'};
SBJ_vars.ch_lab.probe_type = {'seeg','seeg','seeg','seeg','seeg','seeg'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP'};
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]); error('probes ~= type+ref');end;
SBJ_vars.ch_lab.ROI        = {'all',...%'ROF*','FOA*'};
                              '-RHH7-8','-RAM7-8'}; % excluded due to large amplitude artifacts
SBJ_vars.ch_lab.eeg_ROI    = {'CZ'};
SBJ_vars.ch_lab.wires      = {'mram','mrhh','mrth','mrof','mrac'};
SBJ_vars.ch_lab.wire_type  = {'su','su','su','su','su'};
SBJ_vars.ch_lab.wire_ref   = {'','','','',''};
SBJ_vars.ch_lab.wire_ROI   = {'all'};

% SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
% SBJ_vars.ch_lab.suffix = '';    % after every channel except 'EDF Annotations'
% SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ch_lab.nlx_suffix = '_0007';    % after every channel except 'EDF Annotations'
SBJ_vars.ch_lab.nlx_nk_align = {'ROF3','ROF4'}; % tried RPC8,9 I think, maybe emodim: {'RIN4','RIN5'};
SBJ_vars.nlx_macro_inverted  = 1;

SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    'LHH1','LHH2','RHH1',... % epileptic
    'RAM1','RHH2','RHH3','RTH1','LHH3',.... % spread
    'E','DC01','DC02','DC03','DC04','Mark1','Mark2','REF','EKG','Events','GND'... % Junk
    };
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1,1,1,2,2,2,2,0,0,0,0,0,0,0,0,0,0,0];
SBJ_vars.ch_lab.eeg = {'FZ','CZ','OZ','C3','C4'};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
SBJ_vars.ch_lab.eog = {'LUC','LLC','RUC','RLC'};
SBJ_vars.ch_lab.photod = {'photo'};% was changed from original recording "mic1"
SBJ_vars.ch_lab.mic    = {'mic'};% was changed from original recording "spk1"
SBJ_vars.photo_inverted = 1;

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% R1: ~115 to ~1047
SBJ_vars.analysis_time = {{[40.0 1120.0]}};

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
