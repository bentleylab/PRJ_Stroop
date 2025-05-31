%% IR69 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR69';
SBJ_vars.raw_file = {'IR69_raw_R1.mat'};
SBJ_vars.block_name = {'_R1'};% there's a second block I don't want to process right now, so leaving blank here
SBJ_vars.low_srate  = [0, 0];

SBJ_vars.dirs.SBJ     = [root_dir 'PRJ_Stroop/data/' SBJ_vars.SBJ '/'];
SBJ_vars.dirs.raw     = [SBJ_vars.dirs.SBJ '00_raw/'];
SBJ_vars.dirs.SU      = [SBJ_vars.dirs.raw 'SU_2018-02-07_13-17-17/'];
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

SBJ_vars.ch_lab.probes     = {'RHH','RTH','ROF','RIN','FOA','FOP','LES'}; % 'RAM' was all bad
SBJ_vars.ch_lab.probe_type = {'seeg','seeg','seeg','seeg','seeg','seeg','seeg'};
SBJ_vars.ch_lab.ref_type   = {'BP','BP','BP','BP','BP','BP','BP'};
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]); error('probes ~= type+ref');end;
SBJ_vars.ch_lab.ROI        = {'all'};%'ROF*','FOA*'};
SBJ_vars.ch_lab.eeg_ROI    = {};

SBJ_vars.ch_lab.nlx          = [1,1,1,1,1,1,0,0];
SBJ_vars.ch_lab.wires        = {'mram','mrhh','mrth','mrof','mfoa'};
SBJ_vars.ch_lab.wire_type    = {'su','su','su','su','su'};
SBJ_vars.ch_lab.wire_ref     = {'','','','',''};
SBJ_vars.ch_lab.wire_ROI     = {'all'};
SBJ_vars.ch_lab.nlx_suffix   = '_0002';
SBJ_vars.ch_lab.nlx_nk_align = {'FOA3','FOA4'};
SBJ_vars.nlx_macro_inverted  = 1;

% SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
SBJ_vars.ch_lab.suffix = '_0010';    % after every channel except 'EDF Annotations'
% SBJ_vars.ch_lab.mislabel = {{'RLT12','FPG12'},{'IH;L8','IHL8'}};

SBJ_vars.ch_lab.ref_exclude = {}; %exclude from the CAR
SBJ_vars.ch_lab.bad = {...
    };
% emodim .bad:
%     'RAM*',...%crazy bad
%     'RHH1','RHH2','RHH3','RHH7',... % crazy bad
%     'RTH1','RTH2','RTH5','RTH6','RTH7','RTH8',...%crazy bad
%     'ROF3','ROF4','ROF8',...%crazy bad
%     'RIN4','RIN5','RIN6','RIN7','RIN10',...% crazy bad
%     'FOA7',... %crazy bad
%     'FOP1','FOP2','FOP5','FOP6','FOP7',...% crazy bad
%     'LES3','LES4','LES5',...%crazy bad
%     'RHH8','FOA8','FOP8','FOP9','FOP10','LES8','LES9','LES10',... % out of brain
%     'EKG',...% EKG
%     'Mark1','Mark2','REF',...% not real data
%     'DC01','DC02','DC03','DC04','E','GND','Events'...% not real data
SBJ_vars.ch_lab.eeg = {'C3','C4','CZ','FZ','OZ'};
SBJ_vars.ch_lab.eog = {'RUC','RLC','LLC','LUC'};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
SBJ_vars.ch_lab.photod = {'Photo1'};
SBJ_vars.ch_lab.mic    = {'Mic1'};
SBJ_vars.photo_inverted = 1;

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% R1: ~108 to ~1070
% R2: TBD...
SBJ_vars.analysis_time = {{[98.0 1080.0]}};

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
