%% IR48 Processing Variables
[root_dir, ft_dir] = fn_get_root_dir();
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'IR48';
SBJ_vars.raw_file = {'2016092915_0028.besa'};
SBJ_vars.block_name = {''};

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
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_r.mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_v.mat'];
SBJ_vars.recon.elec_mni_s = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_s.mat'];
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

SBJ_vars.ch_lab.probes     = {'FG','AF','MF'};
SBJ_vars.ch_lab.probe_type = {'ecog','ecog','ecog'};
SBJ_vars.ch_lab.ref_type   = {'CAR','CAR','CAR'};
SBJ_vars.ch_lab.ROI        = {'MF*','FG*','AF*'};
SBJ_vars.ch_lab.eeg_ROI    = {'FPZ'};

SBJ_vars.ch_lab.bad = {...
    'FG39','FG40','FG47','FG48','FG55','FG56','FG63','FG64',...%epileptic, Jack (via Jie) said FG47+48
    'AF12',...%noisy
    'DC01','DC03','E','REF','V1','V2','V3','V4','V5','V6','C10','OFG*','PST*','LAM*','LHH*', ...% Not real data
    'EKG'...
    };
%     'RSH','LSH','LLE',...% EOG?
%     'FG12',...% slow rhythm like FG43 but lower amplitude
%     'FG25',...% lots of drift, can maybe be saved with preprocessing?
SBJ_vars.ref_exclude = {...
    'FG27','FG28','FG29','FG37','FG38','FG45','FG46',...% strong mu rhythm that shouldn't be injected into other channels
    'FG12','FG43'...% large slow (~1 Hz) rhythms??? maybe artifact, so watch carefully!
    };
SBJ_vars.ch_lab.eeg = {'FPZ','FP1','FP2','OZ','T5','T6'};
SBJ_vars.ch_lab.FPZ_lap_ref = {'FP1','FP2'};
SBJ_vars.ch_lab.eog = {'RSH','LSH','LLE'}; % Janna says these channels are empty?
SBJ_vars.ch_lab.photod = {'DC02'};
SBJ_vars.ch_lab.mic    = {'DC04'};%burried in noise, so can't visually see responses in plot (can hear just barely)

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
SBJ_vars.notch_freqs = [60 120 180 240 300];
SBJ_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% task stats ~199.5s, ends ~1266s and goes to ~1280s
SBJ_vars.analysis_time = {{[189 1276]}};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 7;
SBJ_vars.artifact_params.hard_threshold_raw = 500;

SBJ_vars.artifact_params.std_limit_diff = 7;
SBJ_vars.artifact_params.hard_threshold_diff = 70;

%--------------------------------------
% Trials to Reject
%--------------------------------------
SBJ_vars.trial_reject_n = [27 165 192 215];
