%% IR57 Preprocessing Pipeline
% This script should be run in sections. Functions/scripts with the SBJ##
% prefix can be run automatically, and all other sections should be
% manually editted for each dataset.

SBJ      = 'IR57';
raw_file = '2017032319_0023.besa';

% Set Up Directories
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults
SBJ_dir     = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/'];
raw_dir     = [SBJ_dir '00_raw/'];
import_dir  = [SBJ_dir '01_import/'];
preproc_dir = [SBJ_dir '02_preproc/'];
events_dir   = [SBJ_dir '03_events/'];
proc_dir    = [SBJ_dir '04_proc/'];
if ~exist(import_dir,'dir')
    mkdir(import_dir);
end
if ~exist(preproc_dir,'dir')
    mkdir(preproc_dir);
end
if ~exist(events_dir,'dir')
    mkdir(events_dir);
end
if ~exist(proc_dir,'dir')
    mkdir(proc_dir);
end

%% ======================================================================== 
%   Step 1- Quick Import and Processing for Data Cleaning/Inspection
%  ========================================================================
plot_psd     = '1by1';
block_prefix = '';
SBJ00_cleaning_prep(SBJ,raw_file,plot_psd,block_prefix);

%% ========================================================================
%   Step 2- Build Preprocessing Variable Structure
%  ========================================================================
%--------------------------------------
% Basics
%--------------------------------------
preproc_vars.SBJ = SBJ;
preproc_vars.raw_filename = strcat(raw_dir,raw_file);
preproc_vars.resamp_freq = 1000;

%--------------------------------------
% Channel Selection
%--------------------------------------
load(strcat(import_dir,SBJ,'_raw_labels.mat'));
preproc_vars.orig_n_ch = length(raw_labels);

ch_labels.ch_lab.probes = {'RSM','RAC','ROF','RIN','RTI','RAM','RHH','RTH',...
                         'LSMA','LAC','LOF','LIN','LTI','LAM','LHH','LTH'};
preproc_vars.ch_lab.ROI    = {'RSM*','RAC*','ROF*','RIN*','RTI*',...
                              'LAM4-5','LAC*','LOF*'}; %LAM4 is inferior anterior insula

% Not necessary to write out all the neural channels like this in each SBJ
% ch_labels.nrl = {
%     'RSM1' 'RSM2' 'RSM3' 'RSM4' 'RSM5' 'RSM6' 'RSM7' 'RSM8' 'RSM9' 'RSM10'...
%     'RAC1' 'RAC2' 'RAC3' 'RAC4' 'RAC5' 'RAC6' 'RAC7' 'RAC8' 'RAC9' 'RAC10'...
%     'ROF1' 'ROF2' 'ROF3' 'ROF4' 'ROF5' 'ROF6' 'ROF7' 'ROF8' 'ROF9' 'ROF10'...
%     'RIN1' 'RIN2' 'RIN3' 'RIN4' 'RIN5' 'RIN6' 'RIN7' 'RIN8' 'RIN9' 'RIN10'...
%     'RTI1' 'RTI2' 'RTI3' 'RTI4' 'RTI5' 'RTI6' 'RTI7' 'RTI8' 'RTI9' 'RTI10'...
%     'RAM1' 'RAM2' 'RAM3' 'RAM4' 'RAM5' 'RAM6' 'RAM7' 'RAM8' 'RAM9' 'RAM10'...
%     'RHH1' 'RHH2' 'RHH3' 'RHH4' 'RHH5' 'RHH6' 'RHH7' 'RHH8' 'RHH9' 'RHH10'...
%     'RTH1' 'RTH2' 'RTH3' 'RTH4' 'RTH5' 'RTH6' 'RTH7' 'RTH8' 'RTH9' 'RTH10'...
%     'LSMA1' 'LSMA2' 'LSMA3' 'LSMA4' 'LSMA5' 'LSMA6' 'LSMA7' 'LSMA8' 'LSMA9' 'LSMA10'...
%     'LAC1' 'LAC2' 'LAC3' 'LAC4' 'LAC5' 'LAC6' 'LAC7' 'LAC8' 'LAC9' 'LAC10'...
%     'LOF1' 'LOF2' 'LOF3' 'LOF4' 'LOF5' 'LOF6' 'LOF7' 'LOF8' 'LOF9' 'LOF10'...
%     'LIN1' 'LIN2' 'LIN3' 'LIN4' 'LIN5' 'LIN6' 'LIN7' 'LIN8' 'LIN9' 'LIN10'...
%     'LTI1' 'LTI2' 'LTI3' 'LTI4' 'LTI5' 'LTI6' 'LTI7' 'LTI8' 'LTI9' 'LTI10'...
%     'LAM1' 'LAM2' 'LAM3' 'LAM4' 'LAM5' 'LAM6' 'LAM7' 'LAM8' 'LAM9' 'LAM10'...
%     'LHH1' 'LHH2' 'LHH3' 'LHH4' 'LHH5' 'LHH6' 'LHH7' 'LHH8' 'LHH9' 'LHH10'...
%     'LTH1' 'LTH2' 'LTH3' 'LTH4' 'LTH5' 'LTH6' 'LTH7' 'LTH8' 'LTH9' 'LTH10'...
%     };

preproc_vars.ch_lab.bad = {...
    'RHH1','RHH2','RHH3','RHH4','ROF1','RAM1','RTH1','RTH2',...%epileptic
    'LHH1','LHH2','LHH3','LHH4','LHH5','LHH6','LHH7','LHH8','LHH9','LHH10',...%epileptic
    'LTH1','LTH2','LTH3','LAM1','LAM2','LAM3',...%epileptic
    'RSM2','RSM3','RIN9','RIN10','LTI2','LTI3','RHH6','LHH1','LHH2','LHH4','LHH6',...%bad line noise
    'RSM8','RSM9','RAC10','RAM10','RTH9','RTH10','LOF10','LAM10',...%out of brain
    'RHH7',...%marked as bad in visual for some reason?
    'DC01','DC03',...% not real data
    'LUC','LLC','RUC','RLC','XREF','E',...% not real data
    'EKG'...
    };
preproc_vars.ch_lab.eeg = {'FPZ' 'CZ' 'OZ' 'C3' 'C4' 'Z' 'FP1' 'FP2' 'T3' 'T4' 'O1' 'O2'};
%                           not sure about 'Z' but PSD looks ok
preproc_vars.ch_lab.photod = {'DC02'};
preproc_vars.ch_lab.mic    = {'DC04'};

%--------------------------------------
% Line Noise Parameters
%--------------------------------------
% most have only regular harmonics with normal width, and 120 and 240 are weak
% RBT, RPIN, RSMA,LAC have an extra peak at 200
% RHH6 has really bad at all harmonics, like LUE and other nonsense chan
preproc_vars.notch_freqs = [60 120 180 240 300]; %200 shoudl come out in re-referencing
preproc_vars.bs_width    = 2;

%--------------------------------------
% Time Parameters
%--------------------------------------
% 9 blocks, starts ~90s with some response like voices before
% ends ~1027s
preproc_vars.analysis_time = {[70 1040]};

%--------------------------------------
% Save Preproc_vars
%--------------------------------------
out_filename = [import_dir SBJ '_preproc_vars.mat'];
save(out_filename,'preproc_vars');


%% ========================================================================
%   Step 3- Import Data, Resample, and Save Individual Data Types
%  ========================================================================
SBJ01_import_data(SBJ,preproc_vars.resamp_freq);


%% ========================================================================
%   Step 4a- Manually Clean Photodiode Trace: Load & Plot
%  ========================================================================
% Load data
evnt_filename = strcat(import_dir,SBJ,'_evnt.mat');
load(evnt_filename);
[evnt, hdr] = fn_format_data_ft2KLA(evnt);

% Plot event channels
plot(linspace(0,hdr.length_in_seconds,hdr.n_samples), evnt);

%% ========================================================================
%   Step 4b- Manually Clean Photodiode Trace: Mark Sections to Correct
%  ========================================================================
% Mark trials to ignore e.g., interruptions
ignore_trials = [1:5];  % see photod comment below

% Set zero/baseline during a block
bsln_val = 2710;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 22.75],...% first 4 trials are missing photodiode; onset of 5th is messed up (toss all 5)
    [107.0 109.0]...% extra blip of a trial at the end of B1
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = {};
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
% wavy stimulus shade on T22 and T31 in B5
stim_times = {[22.96 24.5]};
stim_yval = [30005];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

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
out_filename = [preproc_dir SBJ '_evnt_clean.mat'];
save(out_filename, 'evnt', 'hdr', 'ignore_trials', 'photod_ix', 'mic_ix');

%% ========================================================================
%   Step 5- Parse Event Traces into Behavioral Data
%  ========================================================================
SBJ02_behav_parse(SBJ,[0.300 2.000],ignore_trials,1,1)
% Be sure to save the two figures coming from this function!

%% ========================================================================
%   Step 6- Manually Correct Reaction Times
%  ========================================================================
% Before running this function, open SBJ_events.fig

% Run with default parameters (outlier_thresh=0.1s), not saving plot or trial_info
SBJ03_RT_manual_adjustments(SBJ,[],0,0)
% Check any big discrepancies between auto and manual RTs
% Removed B1T4 because laughing from mistake on B1T3

% Run again after verifying accuracy of manual RTs, saving this time
SBJ03_RT_manual_adjustments(SBJ,0.15,1,1)

%% ========================================================================
%   Step 7- Preprocess Neural Data
%  ========================================================================
SBJ04_preproc(SBJ,'bandstop')

%% ========================================================================
%   Step 8- Reject Bad Trials Based on Behavior and Bob
%  ========================================================================
% Load manually corrected trial_info
clear data trial_info
load(strcat(events_dir,SBJ,'_trial_info_manual.mat'));

% Toss trials based on behavior and cleaning with Bob
event_type = 'stim';
trial_lim_s = [-0.25 2];
RT_std_thresh = 3;
trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,event_type,trial_lim_s,RT_std_thresh);

%% ========================================================================
%   Step 9- Choose Thresholds for Variance-Based Trial Rejection
%  ========================================================================
warning_std_thresh = 3;

% Load data for visualization
load(strcat(preproc_dir,SBJ,'_proc_vars.mat'));
load(strcat(preproc_dir,SBJ,'_preproc_',proc_vars.line_filt_id,'.mat'));

% Select channels of interest
cfg = [];
cfg.channel = proc_vars.ch_lab.ROI;
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
    trial_info_clean.condition_n',trial_lim_s*data.fsample);

% Compute Derivative
cfg = [];
cfg.derivative = 'yes';
trials_dif = ft_preprocessing(cfg,trials);

% Compute potential variance limits over trials and channels
[trial_mat,~] = fn_format_trials_ft2KLA(trials);
var_mat = std(trial_mat,0,3);
ch_var_mean = mean(var_mat,2);
ch_var_thresh = mean(ch_var_mean)+std(ch_var_mean)*warning_std_thresh;

trial_var_mean = mean(var_mat,1);
trial_var_thresh = mean(trial_var_mean)+std(trial_var_mean)*warning_std_thresh;

[trial_mat_dif,~] = fn_format_trials_ft2KLA(trials_dif);
var_mat_dif = std(trial_mat_dif,0,3);
ch_var_mean_dif = mean(var_mat_dif,2);
ch_var_dif_thresh = mean(ch_var_mean_dif)+std(ch_var_mean_dif)*warning_std_thresh;

trial_var_mean_dif = mean(var_mat_dif,1);
trial_var_dif_thresh = mean(trial_var_mean_dif)+std(trial_var_mean_dif)*warning_std_thresh;

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
artifact_params.std_limit_raw = 8;
artifact_params.hard_threshold_raw = 700; % trials go up to ~600-700, then outlier start

artifact_params.std_limit_diff = 10;
artifact_params.hard_threshold_diff = 100; % only one trial above 100 (I think)

% Comment here potential bad trials and channels:
%   Suspect Trial Ix: 134 135 183 185 186 190 299 (same in my warning and ft plot)
%       Diff added: 270, 282, 299 (87,187,189?)
%   Suspect Channels: ROF7-8, LOF9-10, but especially LOF6-7*, LOF7-8*
%       Diff added: RAC8-9, LAC8-9

%% ========================================================================
%   Step 10- Reject Bad Trials Based on Variance
%  ========================================================================
% Run KLA artifact rejection based on robust variance estimates
% If too many/few trials are rejected, adjust artifact_params and rerun
plot_ch = {'worst',5};  %ft_channelselection({ROI_list},data.label); --> RAC, LAC
report.hard_thresh = 2; % print total rejected and trial numbers
report.std_thresh  = 1; % print only total rejected
report.std_plot    = 1; % plot the std distribution and threshold
trial_info_final = SBJ06_reject_artifacts_KLA_report(trials,trial_info_clean,artifact_params,plot_ch,report);

% KLA Notes:
%   Ran with STD thresh 7, lost 107 and 184 diff trials
%       LAC1-2 drove nearly all diff std rejections, but even std_lim=10 is
%       bad so I'll just live with it

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
save(strcat(events_dir,SBJ,'_trial_info_final.mat'),'trial_info');

