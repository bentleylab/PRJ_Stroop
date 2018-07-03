%% Preprocessing Pipeline
% This script should be run in sections. Functions/scripts with the SBJ##
% prefix can be run automatically, and all other sections should be
% manually editted for each dataset.
clear all; close all;

%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Step 0 - Processing Variables
% SBJ = 'IR';

pipeline_id = 'main_ft';
eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

%% ========================================================================
%   Step 1- Load SBJ and Processing Variable Structures
%  ========================================================================
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% ======================================================================== 
%   Step 2- Quick Import and Processing for Data Cleaning/Inspection
%  ========================================================================
% FILE TOO BIG, RUNNING THIS VIA SGE

% block_prefix = '';
% SBJ00_cleaning_prep(SBJ,SBJ_vars.raw_file,proc_vars.plot_psd,block_prefix);

%% ========================================================================
%   Step 3- Import Data, Resample, and Save Individual Data Types
%  ========================================================================
% FILE TOO BIG, RUNNING THIS VIA SGE
% SBJ01_import_data(SBJ,proc_vars.resample_freq);

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
SBJ_evnt_clean_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_evnt_clean/' SBJ '_evnt_clean_params.m'];
eval(SBJ_evnt_clean_cmd);

%% ========================================================================
%   Step 4c- Manually Clean Photodiode Trace: Apply Corrections
%  ========================================================================
photod_ix = strmatch(SBJ_vars.ch_lab.photod,hdr.channel_labels);
mic_ix = strmatch(SBJ_vars.ch_lab.mic,hdr.channel_labels);
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
if proc_vars.resample_freq~=1000
    error('ERROR!!! SBJ02_behav_parse assumes 1 kHz neural sampling rate!!!\n');
end
SBJ02_behav_parse_colinhoy(SBJ,proc_vars.rt_bounds,ignore_trials,1,1)
% Be sure to save the two figures coming from this function!
%   i.e., SBJ_photodiode_segmentation.fig & SBJ_events.fig

%% ========================================================================
%   Step 6- Create einfo based on ROIs
%  ========================================================================
% look at recon and create spreadsheet of general ROI, WM/GM, etc.
%   save that as tsv
fn_compile_einfo(SBJ,pipeline_id)

%% ========================================================================
%   Step 7- Manually Correct Reaction Times
%  ========================================================================
% Send to Rana, get her to create the basic excel sheet
% double check it!

% Before running this function, open SBJ_events.fig
% run with default parameters (outlier_thresh=0.1s), not saving plot or trial_info
outlier_thresh  = 0.15;
save_plot       = 0;
save_trial_info = 0;
SBJ03_RT_manual_adjustments(SBJ,outlier_thresh,save_plot,save_trial_info)
% Check any big discrepancies between auto and manual RTs
% Fix any RTs manually if necessary

% Run again after verifying accuracy of manual RTs, saving this time
save_plot       = 1;
save_trial_info = 1;
SBJ03_RT_manual_adjustments(SBJ,outlier_thresh,save_plot,save_trial_info)

%% ========================================================================
%   Step 8- Preprocess Neural Data
%  ========================================================================
SBJ04_preproc(SBJ,pipeline_id)

%% ========================================================================
%   Step 9- Reject Bad Trials Based on Behavior and Bob
%  ========================================================================
% Load manually corrected trial_info
clear data trial_info
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_manual.mat'));

% Toss trials based on behavior and cleaning with Bob
trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,pipeline_id);

%% ========================================================================
%   Step 10a- Prepare Variance Estimates for Variance-Based Trial Rejection
%  ========================================================================
% Load data for visualization
% load(strcat(preproc_dir,SBJ,'_proc_vars.mat'));
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));

% Select channels of interest
cfg = [];
cfg.channel = SBJ_vars.ch_lab.ROI;
data = ft_selectdata(cfg,data);

% Segment into trials
if strcmp(proc_vars.event_type,'stim')
    events = trial_info_clean.word_onset;
elseif strcmp(proc_vars.event_type,'resp')
    events = trial_info_clean.resp_onset;
else
    error(strcat('ERROR: unknown event_type ',proc_vars.event_type));
end
trials = fn_ft_cut_trials_equal_len(data,events,...
    trial_info_clean.condition_n',proc_vars.trial_lim_s*data.fsample);

% Compute Derivative
cfg = [];
cfg.derivative = 'yes';
trials_dif = ft_preprocessing(cfg,trials);

% Compute potential variance limits over trials and channels
[trial_mat,~] = fn_format_trials_ft2KLA(trials);
var_mat = std(trial_mat,0,3);
ch_var_mean = mean(var_mat,2);
ch_var_thresh = mean(ch_var_mean)+std(ch_var_mean)*proc_vars.var_std_warning_thresh;

trial_var_mean = mean(var_mat,1);
trial_var_thresh = mean(trial_var_mean)+std(trial_var_mean)*proc_vars.var_std_warning_thresh;

[trial_mat_dif,~] = fn_format_trials_ft2KLA(trials_dif);
var_mat_dif = std(trial_mat_dif,0,3);
ch_var_mean_dif = mean(var_mat_dif,2);
ch_var_dif_thresh = mean(ch_var_mean_dif)+std(ch_var_mean_dif)*proc_vars.var_std_warning_thresh;

trial_var_mean_dif = mean(var_mat_dif,1);
trial_var_dif_thresh = mean(trial_var_mean_dif)+std(trial_var_mean_dif)*proc_vars.var_std_warning_thresh;

% Report on potentially bad channels
bad_var_ch      = trials.label(abs(ch_var_mean) > ch_var_thresh);
bad_var_dif_ch  = trials.label(abs(ch_var_mean_dif) > ch_var_dif_thresh);
bad_var_trl     = trial_info_clean.trial_n(abs(trial_var_mean) > trial_var_thresh);
bad_var_dif_trl = trial_info_clean.trial_n(abs(trial_var_mean_dif) > trial_var_dif_thresh);
fprintf('==============================================================================================\n');
fprintf('Simple Variance Rejection:\n');
fprintf('\tChannel Variance Names: %s\n',bad_var_ch{:});
fprintf('\tChannel Diff Variance Names: %s\n',bad_var_dif_ch{:});
fprintf('\tTrial Variance: %i\n',bad_var_trl);
fprintf('\tTrial Diff Variance: %i\n',bad_var_dif_trl);
fprintf('==============================================================================================\n');

% Save results
var_rej_filename = [SBJ_vars.dirs.events SBJ '_variance_rejection_results.txt'];
r_file = fopen(var_rej_filename,'a');
fprintf(r_file,'===============================================================================================\n');
fprintf(r_file,'Simple Variance Rejection:\n');
fprintf(r_file,'Run Time: %s\n',datestr(datetime));
fprintf(r_file,'\tChannel Variance Names: %s\n',bad_var_ch{:});
fprintf(r_file,'\tChannel Diff Variance Names: %s\n',bad_var_dif_ch{:});
fprintf(r_file,'\tTrial Variance: %i\n',bad_var_trl);
fprintf(r_file,'\tTrial Diff Variance: %i\n',bad_var_dif_trl);
fprintf(r_file,'==============================================================================================\n');
fclose(r_file);

%% ========================================================================
%   Step 10b- Choose Thresholds for Variance-Based Trial Rejection
%  ========================================================================
% Visualize data to set limits for variance-based rejection
cfg_reject = [];
cfg_reject.method = 'summary';
ft_rejectvisual(cfg_reject,trials);

% Visualize Derivative
ft_rejectvisual(cfg_reject,trials_dif);

%% ========================================================================
%   Step 11a- Update Rejection parameters and electrodes based on variance
%  ========================================================================
% Comment in evernote note on bad trials and channels!
% Then the following variables should be written into SBJ_vars:

% Update SBJ_vars.ch_lab.var_rej field!

% % Choose thresholds based on plots above
% artifact_params.std_limit_raw = 7;
% artifact_params.hard_threshold_raw = 300; % based on maxabs()
% 
% artifact_params.std_limit_diff = 7;
% artifact_params.hard_threshold_diff = 100; % based on maxabs() for trials_dif

% Re-load SBJ_vars after updating artifact field
clear SBJ_vars
clear_cmd = ['clear ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(clear_cmd); %needed to delete cached version
eval(SBJ_vars_cmd);

% Reload data and re-select channels of interest after excluding bad ones
clear data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
cfg = [];
cfg.channel = SBJ_vars.ch_lab.ROI;
data = ft_selectdata(cfg,data);
trials = fn_ft_cut_trials_equal_len(data,events,...
    trial_info_clean.condition_n',proc_vars.trial_lim_s*data.fsample);

% Compute Derivative
cfg = [];
cfg.derivative = 'yes';
trials_dif = ft_preprocessing(cfg,trials);


%% ========================================================================
%   Step 11b- Automatically Reject Bad Trials Based on Variance
%  ========================================================================
% Run KLA artifact rejection based on robust variance estimates
% If too many/few trials are rejected, adjust artifact_params and rerun

% plot_chans is a cell array of channels to be plotted
%   {} - empty array indicates skip plotting
%   'worst' - says plot channels that had > 5 bad trials
%   full array - list of channel names to plot
% report [struct] - list of 0/1/2 flags for plotting and reporting options,
%   0 = no report; 1 = concise report; 2 = verbose report
%   .hard_thresh: 1=number of trials; 2=print rejected trial_n
%   .std_thresh: 1=number of trials; 2=print rejected trial_n
%   .std_plot: 1=plot std distribution with cut off
plot_ch = {'worst',5};%ft_channelselection({'LPC*','LAC*','RIN*'},data.label);
report.hard_thresh = 2; % print total rejected and trial numbers
report.std_thresh  = 1; % print only total rejected
report.std_plot    = 1; % plot the std distribution and threshold

trial_info_KLA_clean = SBJ06_reject_artifacts_KLA_report(trials,trial_info_clean,...
                                        SBJ_vars.artifact_params,plot_ch,report);

bad_samples = NaN([size(trial_info_KLA_clean.bad_trials.var,1) 2]);
for t_ix = 1:size(bad_samples,1)
    bad_samples(t_ix,:) = trials.sampleinfo(trial_info_clean.trial_n==trial_info_KLA_clean.bad_trials.var(t_ix),:);
end
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'vertical';
cfg.artfctdef.visual.artifact = bad_samples;
ft_databrowser(cfg, trials);

ft_databrowser(cfg, trials_dif);

%% ========================================================================
%   Step 12- Compile Variance-Based Trial Rejection and Save Results
%  ========================================================================
% Re-load SBJ_vars after updating artifact field
clear SBJ_vars
eval(clear_cmd); %needed to delete cached version
eval(SBJ_vars_cmd);

clear trial_info
trial_info = trial_info_clean;
% Document bad trials
trial_info.bad_trials.variance = SBJ_vars.trial_reject_n';
trial_info.bad_trials.all = sort([trial_info.bad_trials.all; trial_info.bad_trials.variance]);

% Remove bad trials
trial_rejected = ismember(trial_info.trial_n,SBJ_vars.trial_reject_n);
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

% Make sure no response times are in weird float format
trial_info.resp_onset = round(trial_info.resp_onset);

save(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

