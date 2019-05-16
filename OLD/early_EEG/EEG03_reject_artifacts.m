%% EEG03_reject_artifacts
SBJ = 'IR39';
pipeline_id = 'eeg_ft';
ecog_pipeline_id = 'main_ft';

%% Add paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load the data
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

eeg_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_cz_lap.mat');
eog_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_eog_filt.mat');
load(eeg_filename);
load(eog_filename);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_manual.mat'));

% Toss trials based on behavior and cleaning with Bob
trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,ecog_pipeline_id);

%% Select channels and events of interest
if strcmp(proc_vars.event_type,'stim')
    events = trial_info_clean.word_onset;
elseif strcmp(proc_vars.event_type,'resp')
    events = trial_info_clean.resp_onset;
else
    error(stract('ERROR: unknown event_type ',proc_vars.event_type));
end

% Convert trial_lim into samples
trial_lim = proc_vars.trial_lim_s*cz.fsample;

%% Toss specifically bad EEG/ROG epochs
bad_eeg_filename = [SBJ_vars.dirs.events SBJ '_eeg_bad_epochs_preclean.mat'];
load(bad_eeg_filename);

% Convert visually bad epochs from full time to analysis_time
%   NOTE: 1 keeps epochs that are only partially overlaping real data
%   (i.e., will be trimmed to edges of analysis_time)
eeg_bad_epochs = fn_convert_epochs_full2at(eeg_bad_epochs,SBJ_vars.analysis_time,...
                                    strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean.mat'),1);
                                
%% Cut data into trials
if strcmp(proc_vars.event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(proc_vars.event_type,'resp')
    events = trial_info.resp_onset;
else
    error(strcat('ERROR: unknown event_type ',proc_vars.event_type));
end
cz_trials = fn_ft_cut_trials_equal_len(cz,events,...
    trial_info.condition_n',proc_vars.trial_lim_s*cz.fsample);
eog_trials = fn_ft_cut_trials_equal_len(eog,events,...
    trial_info.condition_n',proc_vars.trial_lim_s*eog.fsample);

%% Visual Rejection Summary
% Run for CZ
cfg_reject = [];
cfg_reject.method = 'summary';mkdir Rana
cz_clean = ft_rejectvisual(cfg_reject,cz_trials);

% Run for EOG
eog_clean = ft_rejectvisual(cfg_reject,eog_trials);

%% Toss epochs that overlap with bad_epochs from Bob
% Bad epochs from initial visual
skip_epochs_ix = fn_find_trials_overlap_epochs(eeg_bad_epochs,1:size(cz.trial{1},2),events,trial_lim);

% Bad trials from ft_rejectvisual
bad_starts = [cz_clean.cfg.artfctdef.channel.artifact eog_clean.cfg.artfctdef.channel.artifact];
skip_ft_ix = [];
for bad_ix = 1:numel(bad_starts)
    skip_ft_ix = [skip_ft_ix find(bad_starts(bad_ix)==cz_trials.sampleinfo(:,1))];
end

% Check for compatibility between EEG and ECoG proc_vars
trial_info = trial_info_clean;
if ~strcmp(proc_vars.event_type,trial_info_clean.event_type)
    error('EEG and ECoG proc_vars have different event types!');
end
if any(trial_lim~=trial_info_clean.trial_lim)
    error('EEG and ECoG proc_vars have different trial_lim!');
end
if any(proc_vars.trial_lim_s~=trial_info_clean.trial_lim_s)
    error('EEG and ECoG proc_vars have different trial_lim_s!');
end
if proc_vars.RT_std_thresh~=trial_info_clean.RT_std_thresh
    error('EEG and ECoG proc_vars have different RT_std_thresh!');
end

% Document bad trials
trial_info.bad_trials.cz_eog = trial_info.trial_n(skip_ft_ix);
trial_info.bad_trials.eeg_vis = trial_info.trial_n(skip_epochs_ix);
trial_info.bad_trials.all = union(trial_info.bad_trials.all,[trial_info.trial_n(skip_epochs_ix) trial_info.trial_n(skip_ft_ix)]);

% Remove bad trials
trial_info.block_n(skip_eeg_ix)       = [];
trial_info.trial_n(skip_eeg_ix)       = [];
trial_info.word(skip_eeg_ix)          = [];
trial_info.color(skip_eeg_ix)         = [];
trial_info.trialtype(skip_eeg_ix)     = [];
trial_info.blocktype(skip_eeg_ix)     = [];
trial_info.response_time(skip_eeg_ix) = [];
trial_info.marker_time(skip_eeg_ix)   = [];
trial_info.onset_time(skip_eeg_ix)    = [];
trial_info.word_onset(skip_eeg_ix)    = [];
trial_info.resp_onset(skip_eeg_ix)    = [];
trial_info.condition_n(skip_eeg_ix)   = [];
trial_info.error(skip_eeg_ix)         = [];

%prevent formatting errors
trial_info.resp_onset = round(trial_info.resp_onset);
clear trial_info_clean

