not_working_yet_function SBJ06_reject_artifacts(SBJ,trial_info,event_type,trial_lim_s,artifact_params)
% Select data of interest and reject all bad trials and channels to get
% final clean dataset for processing
% Criteria:
%   1. Bob bad epoch
%   2. Bad trial (error, interruption, corrupt data, etc.)
%   3. Bad RT (no response, RT outlier, etc.)
% Inputs:
%   SBJ [str]- the dataset to process (e.g., 'IR54')
%   trial_info [struct]- structure with trial information
%       should be trial_info_clean as output by SBJ05_reject_behav.m
%   event_type [str]- 'stim' or 'resp', determines the event to which trials are locked
%   trial_lim_s [2x1 float array]- boundaries in SECONDS of trials around events
%       trial_lim(1)- baseline, e.g., -0.5 would be 500 ms before event of interest
%       trial_lim(2)- post-event length, e.g., 2 would be 2000 ms after event of interest
%   RT_std_outlier [int]- # standarad deviation from mean for RT to be tossed as outlier
% Outputs:
%   trial_info [struct]- saves out final version after tossing all bad trials

% Parameters
if isempty(trial_lim_s)
    trial_lim = [-0.500 2.5];
end
if isempty(artifact_params)
    artifact_params.std_limit_raw = 3;
    artifact_params.hard_threshold_raw = 200;
    artifact_params.std_limit_dif = 6;
    artifact_params.hard_threshold_dif = 150;
end

% Directories
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
SBJ_dir     = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/'];
raw_dir     = [SBJ_dir '00_raw/'];
import_dir  = [SBJ_dir '01_import/'];
preproc_dir = [SBJ_dir '02_preproc/'];
events_dir = [SBJ_dir '03_events/'];
proc_dir = [SBJ_dir '04_proc/'];
if ~exist(proc_dir,'dir')
    mkdir(proc_dir);
end

%% Load data
load(strcat(preproc_dir,SBJ,'_proc_vars.mat'));
load(strcat(preproc_dir,SBJ,'_preproc_',proc_vars.line_filt_id,'.mat'));
load(strcat(events_dir,SBJ,'_bob_bad_epochs.mat'));

%% Select channels and events of interest
cfg = [];
cfg.channel = proc_vars.ch_lab.ROI;
data = ft_selectdata(cfg,data);

if strcmp(event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(event_type,'resp')
    events = trial_info.resp_onset;
else
    error(stract('ERROR: unknown event_type ',event_type));
end

% Convert trial_lim into samples
trial_lim = trial_lim_s*data.fsample;

%% Segment and Select Trials
trials = fn_ft_cut_trials_equal_len(data,events(ok_trial_ix),trial_info.condition_n(ok_trial_ix)',trial_lim);

% Pull good trials
cfg = [];
cfg.channel = proc_vars.ch_lab.ROI;
trials_ROI = ft_selectdata(cfg,trials);

% Compute Derivative
cfg = [];
cfg.derivative = 'yes';
trials_ROI_dif = ft_preprocessing(cfg,trials_ROI);

%% Visual Inspection
% Compute potential variance limits over trials and channels
[trial_mat,~] = fn_format_trials_ft2KLA(trials_ROI);
[trial_mat_dif,~] = fn_format_trials_ft2KLA(trials_ROI_dif);
std_mat = std(trial_mat,0,3);
std_mat_dif = std(trial_mat_dif,0,3);

ch_std_thresh = std(mean(std_mat,1))*std_thresh;
% figure
% subplot(2,1,1)
% scatter(1:,mean(var_mat,1));
% line([0 324],[std(mean(var_mat,1))*std_thresh std(mean(var_mat,1))*std_thresh],'Color','r','LineStyle','--')
% subplot(2,1,2)
% scatter(mean(var_mat,2),1:75);




%% Compile Bad Trials
fprintf('Num trials excluded for bad RT: %i\n',length(skip_rt1)+length(skip_rt2)+length(skip_rt_outlier));
fprintf('Num trials excluded for errors: %i\n',length(skip_err));
fprintf('Num trials excluded by Bob vis: %i\n',length(skip_bob));
fprintf('Num trials excluded for other : %i\n',length(skip_bad));
fprintf('TOTAL TRIALS EXCLUDED A PRIORI: %i\n',length(skip_trial_ix));
fprintf('TRIALS REMAINING: %i/%i\n',length(ok_trial_ix),length(trial_info.response_time));

trial_info.bad_trials.RT = trial_info.trial_n([skip_rt1 skip_rt2 skip_rt_outlier]);
trial_info.bad_trials.bob = trial_info.trial_n(skip_bob);
trial_info.bad_trials.error = trial_info.trial_n(skip_err);
trial_info.bad_trials.bad = trial_info.trial_n(skip_bad);
trial_info.bad_trials.all = trial_info.trial_n(skip_trial_ix);

% Remove bad trials
trial_info.block_n(skip_trial_ix) = [];
trial_info.trial_n(skip_trial_ix) = [];
trial_info.word(skip_trial_ix) = [];
trial_info.color(skip_trial_ix) = [];
trial_info.trialtype(skip_trial_ix) = [];
trial_info.blocktype(skip_trial_ix) = [];
trial_info.response_time(skip_trial_ix) = [];
trial_info.marker_time(skip_trial_ix) = [];
trial_info.onset_time(skip_trial_ix) = [];
trial_info.word_onset(skip_trial_ix) = [];
trial_info.resp_onset(skip_trial_ix) = [];
trial_info.condition_n(skip_trial_ix) = [];
trial_info.error(skip_trial_ix) = [];

%% Save clean output
%trial_info
%data? I guess not

end
