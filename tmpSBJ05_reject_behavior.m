function trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,pipeline_id)
% Select data of interest and reject all bad trials and channels to get
% final clean dataset for processing
% Criteria:
%   1. Bob bad epoch
%   2. Bad trial (error, interruption, corrupt data, etc.)
%   3. Bad RT (no response, RT outlier, etc.)
% Inputs:
%   SBJ [str]- the dataset to process (e.g., 'IR54')
%   trial_info [struct]- structure with trial information
%       should be SBJ_trial_info_manual.mat from events_dir as saved by SBJ03_RT_manual_adjustments.m
%   event_type [str]- 'stim' or 'resp', determines the event to which trials are locked
%   trial_lim_s [2x1 float array]- boundaries in SECONDS of trials around events
%       trial_lim(1)- baseline, e.g., -0.5 would be 500 ms before event of interest
%       trial_lim(2)- post-event length, e.g., 2 would be 2000 ms after event of interest
%   RT_std_thresh [int]- # standarad deviation from mean for RT to be tossed as outlier
% Outputs:
%   trial_info_clean_behav [struct]- saves out final version after tossing all bad trials

% Directories
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% Load data
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_bob_bad_epochs_preclean.mat'));

% Parameters
if ~isfield(proc_vars,'RT_std_thresh')
    proc_vars.RT_std_thresh = 3;
end
if ~isfield(proc_vars,'trial_lim_s')
    proc_vars.trial_lim_s = [-0.500 2.5];
end

%% Select channels and events of interest
if strcmp(proc_vars.event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(proc_vars.event_type,'resp')
    events = trial_info.resp_onset;
else
    error(stract('ERROR: unknown event_type ',proc_vars.event_type));
end

% Convert trial_lim into samples
trial_lim = proc_vars.trial_lim_s*data.fsample;

%% Reject known artifacts
skip_rt1 = find(trial_info.resp_onset==-1); % this may be unnecessary, not sure if I make any -1 at some point...
skip_rt2 = find(isnan(trial_info.resp_onset));
skip_err = find(trial_info.error==1);
skip_bad = find(trial_info.error==-1);

% Convert visually bad epochs from full time to analysis_time
%   NOTE: 1 indicates trim epoch edges to include only valid points
bad_epochs = fn_convert_epochs_full2at(bad_epochs,SBJ_vars.analysis_time,...
                                    strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean.mat'),1);
% Toss epochs that overlap with bad_epochs from Bob
bob_bad_samples = [];
for bad_epoch_ix = 1:size(bad_epochs,1)
    bob_bad_samples = [bob_bad_samples bad_epochs(bad_epoch_ix,1):bad_epochs(bad_epoch_ix,2)];
end
trials_sample_idx = fn_epoch_cuts_datapad(1:size(data.trial{1},2),events,events,trial_lim);
skip_bob = [];
for epoch_ix = 1:size(trials_sample_idx,1)
    if ~isempty(intersect(trials_sample_idx(epoch_ix,:),bob_bad_samples))
        skip_bob = [skip_bob; epoch_ix];
    end
end

% Find RT outliers
RT_mean = nanmean(trial_info.response_time);
RT_std  = nanstd(trial_info.response_time);
skip_rt_outlier = find(abs(trial_info.response_time-RT_mean)>proc_vars.RT_std_thresh*RT_std);

% Compile all a priori bad trials
skip_trial_ix = unique([skip_bad; skip_rt1; skip_rt2; skip_err; skip_bob; skip_rt_outlier]);
ok_trial_ix = setdiff(1:length(trial_info.resp_onset),skip_trial_ix);

%% Compile Bad Trials
trial_info_clean = trial_info;
trial_info_clean.event_type = proc_vars.event_type;
trial_info_clean.trial_lim = trial_lim;
trial_info_clean.trial_lim_s = proc_vars.trial_lim_s;
trial_info_clean.RT_std_thresh = proc_vars.RT_std_thresh;

% Document bad trials
trial_info_clean.bad_trials.RT = trial_info.trial_n([skip_rt1 skip_rt2 skip_rt_outlier]);
trial_info_clean.bad_trials.bob = trial_info.trial_n(skip_bob);
trial_info_clean.bad_trials.error = trial_info.trial_n(skip_err);
trial_info_clean.bad_trials.bad = trial_info.trial_n(skip_bad);
trial_info_clean.bad_trials.all = trial_info.trial_n(skip_trial_ix);

% Remove bad trials
trial_info_clean.block_n(skip_trial_ix) = [];
trial_info_clean.trial_n(skip_trial_ix) = [];
trial_info_clean.word(skip_trial_ix) = [];
trial_info_clean.color(skip_trial_ix) = [];
trial_info_clean.trialtype(skip_trial_ix) = [];
trial_info_clean.blocktype(skip_trial_ix) = [];
trial_info_clean.response_time(skip_trial_ix) = [];
trial_info_clean.marker_time(skip_trial_ix) = [];
trial_info_clean.onset_time(skip_trial_ix) = [];
trial_info_clean.word_onset(skip_trial_ix) = [];
trial_info_clean.resp_onset(skip_trial_ix) = [];
trial_info_clean.condition_n(skip_trial_ix) = [];
trial_info_clean.error(skip_trial_ix) = [];

% Print results
fprintf('==============================================================================================\n');
fprintf('Num trials excluded for bad RT: %i\n',length(skip_rt1)+length(skip_rt2)+length(skip_rt_outlier));
fprintf('Num trials excluded for errors: %i\n',length(skip_err));
fprintf('Num trials excluded by Bob vis: %i\n',length(skip_bob));
fprintf('Num trials excluded for other : %i\n',length(skip_bad));
fprintf('TOTAL TRIALS EXCLUDED A PRIORI: %i\n',length(skip_trial_ix));
fprintf('TRIALS REMAINING: %i/%i\n',length(trial_info_clean.trial_n),length(trial_info.response_time));
fprintf('==============================================================================================\n');

end
