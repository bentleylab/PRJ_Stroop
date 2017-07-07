function A03_combine_blocks(subject_id, input_file_strs)
% Make input_file_strs the same order as the actual blocks
% Assumes all blocks have the same number of channels, sample rate, subject_id, etc. Only number of trials and points should be different

%% File paths
helper_function_dir_name = '/_TOOLBOXES';
data_dir_name = ['../' subject_id];
datain_dir_name = '01_import';
event_dir_name = '02_events';
output_dir_name = '03_fulldata';
if(~exist(fullfile(pwd,helper_function_dir_name)))
  fprintf('\nERROR: Helper function directory does not exist.\n');
  return;
end

% Import helper functions to Matlab path
addpath(genpath(fullfile(pwd,helper_function_dir_name)));

% Define path to input data and output data
data_path = fullfile(pwd,data_dir_name,datain_dir_name);
event_path = fullfile(pwd,data_dir_name,event_dir_name);
%   create directory if needed
output_path = fullfile(pwd,data_dir_name,output_dir_name);
if(exist(output_path, 'dir') == 0)
  mkdir(output_path);
end

fprintf('%s - Combining blocks:\n', subject_id);
for file_n = 1:length(input_file_strs)
  fprintf('\t%s\n', input_file_strs{file_n});
end

%% Loop through input files and combine them into one

% Get total number of trials and samples
total_samples = 0;
total_trials = 0;
fprintf('.');
for file_n = 1:length(input_file_strs)
  % Load ECoG data
  load(fullfile(data_path,input_file_strs{file_n})); clear data_evnt; clear data_ecog;
  total_samples = total_samples + header_ecog.n_samples;
  block_headers{file_n}.header_ecog = header_ecog;
  block_headers{file_n}.header_evnt = header_evnt;
  n_channels = header_ecog.n_channels;
  % Load trial data
  load(fullfile(event_path,input_file_strs{file_n}));
  total_trials = total_trials + length(trial_info.trial_n);
  block_headers{file_n}.trial_info = trial_info;
  block_headers{file_n}.input_file_name = input_file_strs{file_n};
end

% Initialize variables
all_data_ecog = zeros(n_channels,total_samples);
all_trial_info.sample_rate = trial_info.sample_rate;
all_trial_info.rt_window = trial_info.rt_window;
all_trial_info.ignore_trials = trial_info.ignore_trials;
all_trial_info.log_dir_struct = trial_info.log_dir_struct;
all_trial_info.subject_id = trial_info.subject_id;
all_trial_info.condition_types = trial_info.condition_types;
all_trial_info.condition_n = [];
all_trial_info.block_n = [];
all_trial_info.trial_n = [];
all_trial_info.response_time = [];
all_trial_info.marker_time = [];
all_trial_info.onset_time = [];
all_trial_info.word_onset = [];
all_trial_info.resp_onset = [];
all_trial_info.word = {};
all_trial_info.color = {};
all_trial_info.trialtype = {};
all_trial_info.blocktype = {};

% Loop through files
last_sample_n = 0;
last_block_n = 0;
last_trial_n = 0;
fprintf('.');
for file_n = 1:length(input_file_strs)
  % Load ECoG data
  load(fullfile(data_path,input_file_strs{file_n})); clear data_evnt;
  all_data_ecog(:,(last_sample_n+1):(last_sample_n+header_ecog.n_samples)) = data_ecog; clear data_ecog;
  % Load trial data
  load(fullfile(event_path,input_file_strs{file_n}));
  all_trial_info.block_n = [all_trial_info.block_n; trial_info.block_n+last_block_n];
  last_block_n = max(all_trial_info.block_n);
  all_trial_info.trial_n = [all_trial_info.trial_n; trial_info.trial_n+last_trial_n];
  last_trial_n = max(all_trial_info.trial_n);
  all_trial_info.condition_n = [all_trial_info.condition_n; trial_info.condition_n];
  all_trial_info.response_time = [all_trial_info.response_time trial_info.response_time];
  all_trial_info.marker_time = [all_trial_info.marker_time trial_info.marker_time];
  all_trial_info.onset_time = [all_trial_info.onset_time trial_info.onset_time];
  all_trial_info.word_onset = [all_trial_info.word_onset; trial_info.word_onset+last_sample_n];
  all_trial_info.resp_onset = [all_trial_info.resp_onset; trial_info.resp_onset+last_sample_n];
  all_trial_info.word = {all_trial_info.word{:}, trial_info.word{:}};
  all_trial_info.color = {all_trial_info.color{:}, trial_info.color{:}};
  all_trial_info.trialtype = {all_trial_info.trialtype{:}, trial_info.trialtype{:}};
  all_trial_info.blocktype = {all_trial_info.blocktype{:}, trial_info.blocktype{:}};
  % Increment sample_n
  last_sample_n = last_sample_n+header_ecog.n_samples;
end

% Fix header_ecog
header_ecog = rmfield(header_ecog, {'recording_startdate', 'recording_starttime','orig_samples_extracted','pid','rid'});
header_ecog.n_samples = total_samples;
header_ecog.length_in_seconds = header_ecog.n_samples/header_ecog.sample_rate;

% Rename variables
data_ecog = all_data_ecog; clear all_data_ecog;
trial_info = all_trial_info; clear all_trial_info;

%% Save output
fprintf('.');
save(fullfile(output_path,[subject_id '.mat']), 'subject_id', 'data_ecog', 'header_ecog', 'trial_info', 'block_headers');
fprintf('\n');

%% Delete inut files
for file_n = 1:length(input_file_strs)
  % Load ECoG data
  delete(fullfile(data_path,input_file_strs{file_n}));
end













































































