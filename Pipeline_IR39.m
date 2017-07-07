%% Preprocessing Pipeline: 'main_ft'
% This script should be run in sections. Functions/scripts with the SBJ##
% prefix can be run automatically, and all other sections should be
% manually editted for each dataset.
clear all; close all;
% Set Up Directories
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Processing Variables
SBJ = 'IR39';

pipeline_name = 'main_ft';

% Data Preprocessing
plot_psd      = '1by1';         % type of plot for channel PSDs
resample_freq = 1000;
notch_type    = 'bandstop';     % method for nothc filtering out line noise

% Behavioral Processing
rt_bounds = [0.3 2.0];          % bounds on a reasonable RT to be detected with KLA algorithm

% Trial Cut Parameteres
event_type    = 'stim';         % 'stim'/'resp': lock trial to these event
trial_lim_sec = [-0.25 2];      % data segments (in seconds) to grab around events
RT_std_thresh = 3;              % rejection threshold for RTs

% Varaince-Based Trial Rejection Parameters
var_std_warning_thresh = 3;

%% ========================================================================
%   Step 1- Load SBJ and Processing Variable Structures
%  ========================================================================
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% ======================================================================== 
%   Step 2- Quick Import and Processing for Data Cleaning/Inspection
%  ========================================================================
% FILE TOO BIG, RUNNING THIS VIA SGE

% block_prefix = '';
% SBJ00_cleaning_prep(SBJ,SBJ_vars.raw_file,plot_psd,block_prefix);

%% ========================================================================
%   Step 3- Import Data, Resample, and Save Individual Data Types
%  ========================================================================
SBJ01_import_data(SBJ,resample_freq);

%% ========================================================================
%   Step 4a- Manually Clean Photodiode Trace: Load & Plot
%  ========================================================================
% Load data
evnt_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt.mat');
load(evnt_filename);
[evnt, hdr] = fn_format_data_ft2KLA(evnt);

% Plot event channels
plot(linspace(0,hdr.length_in_seconds,hdr.n_samples), evnt);

%% ========================================================================
%   Step 4b- Manually Clean Photodiode Trace: Mark Sections to Correct
%  ========================================================================
% Create correction times and values in a separate file in ~/PRJ_Stroop/scripts/SBJ_evnt_clean/
SBJ_evnt_clean_filename = ['/home/knight/PRJ_Stroop/scripts/SBJ_evnt_clean/' SBJ '_evnt_clean_params.m'];
eval(SBJ_evnt_clean_filename);

%% ========================================================================
%   Step 4c- Manually Clean Photodiode Trace: Apply Corrections
%  ========================================================================
photod_ix = strmatch(preproc_vars.ch_lab.photod,hdr.channel_labels);
mic_ix = strmatch(preproc_vars.ch_lab.mic,hdr.channel_labels);
% Correct baseline shift
for shift_ix = 1:length(bsln_shift_times)
    epoch_idx = floor(bsln_shift_times{shift_ix}(1)*hdr.sample_rate):floor(bsln_shift_times{shift_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = evnt(photod_ix,epoch_idx) - bsln_shift_val(shift_ix);
end
% zero out drifts
for zero_ix = 1:length(bsln_times)
    epoch_idx = floor(bsln_times{zero_ix}(1)*hdr.sample_rate):floor(bsln_times{zero_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = bsln_val;
end

% level out stimulus periods
for stim_ix = 1:length(stim_times)
    epoch_idx = floor(stim_times{stim_ix}(1)*hdr.sample_rate):floor(stim_times{stim_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = stim_yval(stim_ix);
end

% Save corrected data
out_filename = [SBJ_vars.dirs.preproc SBJ '_evnt_clean.mat'];
save(out_filename, 'evnt', 'hdr', 'ignore_trials', 'photod_ix', 'mic_ix');

%% ========================================================================
%   Step 5- Parse Event Traces into Behavioral Data
%  ========================================================================
SBJ02_behav_parse(SBJ,rt_bounds,ignore_trials,1,1)
% Be sure to save the two figures coming from this function!
%   i.e., SBJ_photodiode_segmentation.fig & SBJ_events.fig

%% ========================================================================
%   Step 6- Manually Correct Reaction Times
%  ========================================================================
% Before running this function, open SBJ_events.fig
% run with default parameters (outlier_thresh=0.1s), not saving plot or trial_info
SBJ03_RT_manual_adjustments(SBJ,[],0,0)
% Check any big discrepancies between auto and manual RTs

% Run again after verifying accuracy of manual RTs, saving this time
SBJ03_RT_manual_adjustments(SBJ,0.15,1,1)

%% ========================================================================
%   Step 7- Preprocess Neural Data
%  ========================================================================
SBJ04_preproc(SBJ,notch_type)

%% ========================================================================
%   Step 8- Reject Bad Trials Based on Behavior and Bob
%  ========================================================================
% Load manually corrected trial_info
clear data trial_info
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_manual.mat'));

% Toss trials based on behavior and cleaning with Bob
trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,event_type,trial_lim_sec,RT_std_thresh);

%% ========================================================================
%   Step 9- Choose Thresholds for Variance-Based Trial Rejection
%  ========================================================================
% Load data for visualization
% load(strcat(preproc_dir,SBJ,'_proc_vars.mat'));
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_name,'.mat'));

% Select channels of interest
cfg = [];
cfg.channel = SBJ_vars.ch_lab.ROI;
data = ft_selectdata(cfg,data);

% Segment into trials
if strcmp(event_type,'stim')
    events = trial_info_clean.word_onset;
elseif strcmp(event_type,'resp')
    events = trial_info_clean.resp_onset;
else
    error(stract('ERROR: unknown event_type ',event_type));
end
trials = fn_ft_cut_trials_equal_len(data,events,...
    trial_info_clean.condition_n',trial_lim_sec*data.fsample);

% Compute Derivative
cfg = [];
cfg.derivative = 'yes';
trials_dif = ft_preprocessing(cfg,trials);

% Compute potential variance limits over trials and channels
[trial_mat,~] = fn_format_trials_ft2KLA(trials);
var_mat = std(trial_mat,0,3);
ch_var_mean = mean(var_mat,2);
ch_var_thresh = mean(ch_var_mean)+std(ch_var_mean)*var_std_warning_thresh;

trial_var_mean = mean(var_mat,1);
trial_var_thresh = mean(trial_var_mean)+std(trial_var_mean)*var_std_warning_thresh;

[trial_mat_dif,~] = fn_format_trials_ft2KLA(trials_dif);
var_mat_dif = std(trial_mat_dif,0,3);
ch_var_mean_dif = mean(var_mat_dif,2);
ch_var_dif_thresh = mean(ch_var_mean_dif)+std(ch_var_mean_dif)*var_std_warning_thresh;

trial_var_mean_dif = mean(var_mat_dif,1);
trial_var_dif_thresh = mean(trial_var_mean_dif)+std(trial_var_mean_dif)*var_std_warning_thresh;

% Report on potentially bad channels
fprintf('==============================================================================================\n');
fprintf('WARNING! Check channels over thresholds:\n');
fprintf('\tChannel Variance Names:');
disp(trials.label(abs(ch_var_mean) > ch_var_thresh));
fprintf('\n\tChannel Diff Variance Names:');
disp(trials.label(abs(ch_var_mean_dif) > ch_var_dif_thresh));
fprintf('\n\tTrial Variance Ix:');
disp(find(abs(trial_var_mean) > trial_var_thresh));
fprintf('\n\tTrial Diff Variance Ix:');
disp(find(abs(trial_var_mean_dif) > trial_var_dif_thresh));
fprintf('\n');
fprintf('==============================================================================================\n');

% Visualize data to set limits for variance-based rejection
cfg_reject = [];
cfg_reject.method = 'summary';
ft_rejectvisual(cfg_reject,trials);

% Visualize Derivative
ft_rejectvisual(cfg_reject,trials_dif);

% Choose thresholds based on plots above
artifact_params.std_limit_raw = 7;
artifact_params.hard_threshold_raw = 300; % trials go up to ~300, then outlier start

artifact_params.std_limit_diff = 7;
artifact_params.hard_threshold_diff = 100; % only one trial above 100 (I think)

% Comment here potential bad trials and channels:
%   Suspect Trials: 42, 118, 172***
%   Suspect Channels: LIN5-6, RIN5-6, ROF6-7

%% ========================================================================
%   Step 10- Reject Bad Trials Based on Variance
%  ========================================================================
% Re-load SBJ_vars after updating artifact field
clear SBJ_vars
eval(SBJ_vars_cmd);

% Run KLA artifact rejection based on robust variance estimates
% If too many/few trials are rejected, adjust artifact_params and rerun
plot_ch = {'worst',5};%ft_channelselection({'LPC*','LAC*','RIN*'},data.label);
report.hard_thresh = 2; % print total rejected and trial numbers
report.std_thresh  = 1; % print only total rejected
report.std_plot    = 1; % plot the std distribution and threshold
trial_info_final = SBJ06_reject_artifacts_KLA_report(trials,trial_info_clean,artifact_params,plot_ch,report);

% KLA Notes:
%   LIN5-6, 6-7, and 7-8 and RIN7-8 all end up with a lot of rejected trials
%   LAC1-2 has a few, but they might be interesting!
%   trial 202 has BIG fluctuations, good catch; 46 also not great

bad_samples = NaN([size(trial_info_final.bad_trials.var,1) 2]);
for t_ix = 1:size(bad_samples,1)
    bad_samples(t_ix,:) = trials.sampleinfo(find(trial_info_clean.trial_n==trial_info_final.bad_trials.var(t_ix)),:);
end
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'vertical';
cfg.artfctdef.visual.artifact = bad_samples;
ft_databrowser(cfg, trials);

%% ========================================================================
%   Step 10- Save the final trial_info
%  ========================================================================
clear trial_info
trial_info = trial_info_final;
save(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

