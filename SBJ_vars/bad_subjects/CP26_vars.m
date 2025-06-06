%% CP26 Processing Variables
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
if isempty(strfind(path,'fieldtrip'))
    addpath(ft_dir);
    ft_defaults
end

%--------------------------------------
% Basics
%--------------------------------------
SBJ_vars.SBJ = 'CP26';
SBJ_vars.raw_file = {'CP26_stroop_raw.mat'};
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
SBJ_vars.recon.elec_pat   = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_acpc_fr.mat'];
SBJ_vars.recon.elec_mni_v = [SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_frv.mat'];
SBJ_vars.recon.elec_mni_s = [];%SBJ_vars.dirs.recon 'Electrodes/' SBJ_vars.SBJ '_elec_mni_s.mat'];
SBJ_vars.recon.fs_T1      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_T1.mgz'];
SBJ_vars.recon.fs_DK      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc+aseg.mgz'];
SBJ_vars.recon.fs_Dx      = [SBJ_vars.dirs.recon 'Scans/' SBJ_vars.SBJ '_fs_preop_aparc.a2009s+aseg.mgz'];

%--------------------------------------
% Channel Selection
%--------------------------------------
SBJ_vars.ch_lab.probes     = {'LIHA','LIHM','LIHP','RIHA','RIHM','RIHP','LOF','ROF','LG','RG'};%'A' is noise
SBJ_vars.ch_lab.probe_type = {'ecog','ecog','ecog','ecog','ecog','ecog','ecog','ecog','ecog','ecog'};
SBJ_vars.ch_lab.ref_type   = {'CARall','CARall','CARall','CARall','CARall','CARall','CARall','CARall','CARall','CARall'};
SBJ_vars.ch_lab.ref_group  = [1 1 1 2 2 2 3 4 3 4];
if ~all(numel(SBJ_vars.ch_lab.probes)==[numel(SBJ_vars.ch_lab.probe_type) numel(SBJ_vars.ch_lab.ref_type)]); error('probes ~= type+ref');end;
SBJ_vars.ch_lab.ROI        = {'all'};
SBJ_vars.ch_lab.eeg_ROI    = {};

%SBJ_vars.ch_lab.prefix = 'POL ';    % before every channel except 'EDF Annotations'
%SBJ_vars.ch_lab.suffix = '-Ref';    % after every channel except 'EDF Annotations'
SBJ_vars.ch_lab.mislabel = {...
    {'LTHP1','LIHP1'},{'LTHP2','LIHP2'},{'LTHP3','LIHP3'},...%,{'LTHP4','LIHP4'}}; (bad)
    {'LIH1','LIHM1'},{'LIH2','LIHM2'},{'LIH3','LIHM3'},{'LIH4','LIHM4'},...%otherwise LIH selects LIHA+LIHP
    {'RIH1','RIHM1'},{'RIH2','RIHM2'},{'RIH3','RIHM3'},{'RIH4','RIHM4'}...%otherwise RIH selects RIHA+RIHP
    };

SBJ_vars.ch_lab.ref_exclude = {
    'LOF3','LOF4',...%main leaders of network
    'LG9','LG17','LG26','LG28','LG20',...%anterior network of spike + slow
    'LG12','LG11','LG14',... % middle slowing network
    'LG8','LG16','LG24','LG32',...%posterior spiking network
    'ROF2','ROF4','RIHA1','RIHA2','RG33','RG34','RG36','RG37','RG41','RG43','RG44','RG59','RG60',...%anterior network
    'RG39','RG40','RG48','RG64','RG55','RG63'...%posterior network
    }; %exclude from the CAR
% old cleaning notes from others (first pass):
%     'RIH1','RIH3','ROF8',...%Christina said loose? I think they're fine, will likely keep
%     'LG24','ROF2','LOF3','LOF4','LIHA4','LIHP4','RIHA1','RIHA3','RIHA4','RG39','RG56','RG63','RG64'...%spikes?
SBJ_vars.ch_lab.bad = {...
    'LTHP4',...%spiking
    'LG25','RG35',...% most egregious spikers/slowers
    'LG1','LG2',...% HF noise
    'RG49',...% loose
    'RIHP1','RIHP2',... % empty channels, negative mirrors w/ perfect anticorrelation
    'LOF5','LOF6','LOF7','LOF8','ROF5','ROF6','ROF7','ROF8',...%on top of grid, no signal
    'A*','DC03','DC04','E','Mark1','Mark2','Events'...%not real data
    };
% bad_codes: 1 = toss (epileptic or bad); 2 = suspicious; 3 = out of brain; 0 = junk
SBJ_vars.ch_lab.bad_type = {'bad','sus','out'};
SBJ_vars.ch_lab.bad_code = [1 2 2 2 2 2 3 3 3 3 3 3 3 3 3 3 0 0 0 0 0 0 0];
if numel(SBJ_vars.ch_lab.bad)~=numel(SBJ_vars.ch_lab.bad_code);error('bad ~= bad_code');end
SBJ_vars.ch_lab.eeg = {};
% SBJ_vars.ch_lab.CZ_lap_ref = {};
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
SBJ_vars.analysis_time = {{[1 990]}};

%--------------------------------------
% Artifact Rejection Parameters
%--------------------------------------
SBJ_vars.artifact_params.std_limit_raw = 10;
SBJ_vars.artifact_params.hard_threshold_raw = 1000;

SBJ_vars.artifact_params.std_limit_diff = 10;
SBJ_vars.artifact_params.hard_threshold_diff = 100;

%--------------------------------------
% Trials to Reject
%--------------------------------------
SBJ_vars.trial_reject_n = [];
