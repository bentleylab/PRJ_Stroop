%% Save EEG to fiff
%   1. Load EEG bad_epochs and create trial_info for EEG analysis
%   2. Load EEG data and export to .fiff
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

SBJ = 'IR57';%'IR57';
pipeline_id = 'main_ft';

%% Load data
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

load([SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat']);
load(strcat(SBJ_vars.dirs.events,SBJ,'_eeg_bad_epochs_preclean.mat'));

load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
load(strcat(SBJ_vars.dirs.import,SBJ,'_eeg_',num2str(proc_vars.resample_freq),'hz.mat'));

%% Select events of interest
if strcmp(proc_vars.event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(proc_vars.event_type,'resp')
    events = trial_info.resp_onset;
else
    error(stract('ERROR: unknown event_type ',proc_vars.event_type));
end

% Convert trial_lim into samples
if ~isfield(proc_vars,'trial_lim_s')
    proc_vars.trial_lim_s = [-0.500 2.5];
end
trial_lim = proc_vars.trial_lim_s*data.fsample;

%% Toss bad epochs
% NOTE: skip_eeg will not have bob bad trials because those are already
%       gone from trial_info_final
eeg_bad_epochs = fn_convert_epochs_full2at(eeg_bad_epochs,SBJ_vars.analysis_time,...
                                    strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean.mat'),1);
skip_eeg = fn_find_trials_overlap_epochs(eeg_bad_epochs,1:size(data.trial{1},2),events,trial_lim);

%% Save trial_info_final_eeg
trial_info_final = trial_info;

% Document bad EEG trials
trial_info.bad_trials.eeg = skip_eeg;
trial_info.bad_trials.all = sort([trial_info.bad_trials.all; trial_info.bad_trials.eeg]);

% Remove bad trials
trial_rejected = ismember(trial_info.trial_n,trial_info.bad_trials.all);
trial_reject_ix = find(trial_rejected);
trial_info.block_n(trial_reject_ix) = [];
trial_info.trial_n(trial_reject_ix) = [];
trial_info.word(trial_reject_ix) = [];
trial_info.color(trial_reject_ix) = [];
trial_info.trialtype(trial_reject_ix) = [];
trial_info.blocktype(trial_reject_ix) = [];
trial_info.response_time(trial_reject_ix) = [];
trial_info.marker_time(trial_reject_ix) = [];
trial_info.onset_time(trial_reject_ix) = [];
trial_info.word_onset(trial_reject_ix) = [];
trial_info.resp_onset(trial_reject_ix) = [];
trial_info.condition_n(trial_reject_ix) = [];
trial_info.error(trial_reject_ix) = [];

% Get rid of non-matrix entries
trial_info = rmfield(trial_info,'word');
trial_info = rmfield(trial_info,'color');
trial_info = rmfield(trial_info,'trialtype');
trial_info = rmfield(trial_info,'blocktype');
trial_info = rmfield(trial_info,'condition_types');
trial_info = rmfield(trial_info,'sample_rate');
trial_info = rmfield(trial_info,'rt_window');
trial_info = rmfield(trial_info,'ignore_trials');
trial_info = rmfield(trial_info,'SBJ');
trial_info = rmfield(trial_info,'event_type');
trial_info = rmfield(trial_info,'trial_lim');
trial_info = rmfield(trial_info,'trial_lim_s');
trial_info = rmfield(trial_info,'RT_std_thresh');
trial_info = rmfield(trial_info,'bad_trials');

% Make sure no response times are in weird float format
trial_info.resp_onset = round(trial_info.resp_onset);

save(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final_eeg.mat'),'-v7.3','trial_info');

%% add fields tro FT structure
for e_ix = 1:numel(eeg.hdr.chantype)
    eeg.hdr.chantype{e_ix} = 'eeg';
    eeg.hdr.chanunit{e_ix} = '';
end

%% Save eeg and trial_info to .fiff
fiff_file  = [SBJ_vars.dirs.import SBJ '_eeg_' num2str(proc_vars.resample_freq) 'hz.fif'];
fieldtrip2fiff(fiff_file, eeg);
