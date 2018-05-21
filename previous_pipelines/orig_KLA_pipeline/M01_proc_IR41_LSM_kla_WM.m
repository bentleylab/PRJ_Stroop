%% DATA_ID: IR41_LSM_kla_WM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run sections individually to process the dataset
close all
clear all

SBJ = 'IR41';
SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
raw_data_id = 'IR41_LSM_kla';
input_file = strcat(raw_data_id, '.mat');
input_filename = strcat(SBJ_dir,'01_import/',input_file);
proc_vars_filename = strcat(SBJ_dir,'01_import/',raw_data_id,'_proc_vars.mat');
preproc_id = strcat(raw_data_id,'_WM');
preproc_file = strcat(preproc_id,'.mat');
% output_dir = strcat(SBJ_dir,'03_proc/'); now specified in A03_comb blocks

% Filtering parameters
hp_cutoff = 0.5;

% Add toolboxes to path
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath('/home/knight/hoycw/Apps/fieldtrip-20160927/'); % read .besa
ft_defaults

%% Look at noise profile of all channels
% Load data
load(proc_vars_filename);
load(input_filename);

% Loop through channels to see noise profile
for channel_n = 1:header_ecog.n_channels
    [fft_data,freqs] = pwelch(data_ecog(channel_n,:),2048,0,2048,header_ecog.sample_rate);
    loglog(freqs,fft_data);
    xlim([1 350]);
    ax = gca;
    ax.XTick = [25 30 60 100 120 180 200 240 300 360 420];
    title(['Channel ' num2str(channel_n) ' [' header_ecog.channel_labels{channel_n} ']']);
    pause;
end
close all;
%% Set line filters based on noise profiles plotted above
% Looks like there is noise in 25-30 (???!), nothing at 60, then again at
% 180, 300 Hz...
line_noise_freqs = [60 120 180 300];
dummy_line_noise_freqs = []; %FIX ME!!! because cleanline2.m doesn't work

%% Importing all data blocks

% Re-import the data using the information we just noted
%  Also highpass filter which is slow
%   !!!CWH: will detrend (demean, only linear) if no high pass; higher polynomial fits???

A01_import_reref_filter_ft_besa_good_time(preproc_id, input_file, preproc_file, analysis_channels, event_channels, ...
  'SBJ_dir_name', SBJ_dir, ...
  'resample_rate', 1000, ...
  'reref_channels', reref_channels, ...
  'reref_name', 'WM', ...
  'reref_weights', reref_weights, ...
  'line_noise_freqs', dummy_line_noise_freqs, ...
  'analysis_time', analysis_time, ...
  'highpass_freq', hp_cutoff, ...
  'save_output', 1 ...
  );

%% use ft to remove line noise since cleanline2 failed...
clear header_ecog data_ecog
preproc_filename = strcat(SBJ_dir,'02_preproc/',preproc_file);
load(preproc_filename);
data_notch = [];
data_notch.label = header_ecog.channel_labels;
data_notch.fsample = header_ecog.sample_rate;
% NOTE: other functions will assume your data is in the trial field
data_notch.trial{1} = data_ecog; % data appears here by default before you do anything
data_notch.time{1} = linspace(0, size(data_ecog,2), size(data_ecog,2)); % ms time spacing

ft_line_noise_freqs = line_noise_freqs;
cfg = [];
cfg.dftfilter = 'yes';
cfg.dftfreq = ft_line_noise_freqs;
data_notch = ft_preprocessing(cfg, data_notch);

% Convert back to mat...
data_ecog = data_notch.trial{:};
header_ecog.line_noise_freqs = ft_line_noise_freqs;

%% Check event channels
% Look at event channels for Block 1
plot(linspace(0,header_evnt.length_in_seconds,header_evnt.n_samples), data_evnt);
% The photodiode data between Block 1 and Block 2 is bad

% Manually correct bad event traces
% Fix IR21: before B1, b/w B1/B2, first two trials of B2, very end
% photo_channel_n = 1;
yval = 6500; %this is the zero/baseline during a block
% before cutting to analysis_time: [0.0 161.0 261.0 421.0 751.0 757.0];
set_times = [0.0 21.0 122.0 140.0 240.0 280.0 606.0 616.0]; % in sec
data_evnt(photo_channel_n,1:floor(set_times(2)*header_evnt.sample_rate)) = yval;
data_evnt(photo_channel_n,floor(set_times(3)*header_evnt.sample_rate):floor(set_times(4)*header_evnt.sample_rate)) = yval;
data_evnt(photo_channel_n,floor(set_times(5)*header_evnt.sample_rate):floor(set_times(6)*header_evnt.sample_rate)) = yval;
data_evnt(photo_channel_n,floor(set_times(7)*header_evnt.sample_rate):floor(set_times(8)*header_evnt.sample_rate)) = yval;
save(preproc_filename, 'header_ecog', 'data_ecog', 'header_evnt', 'data_evnt', 'data_id');


%% Load events
% mic_channel_n = 2;
rt_window = [0.300 2.000]; % Possible response times from word onset in seconds
plot_it = 1; % Set to 0 if you don't want to see plots
save_it = 0; % 0 if yo don't want to write files
log_dir_struct.log_file_str = strcat(SBJ,'_strooptask_log.txt');
log_dir_struct.wav_file_dir = strcat(SBJ,'_strooptask_wavfiles');
% ingore_trials is for the log file; those should already be zeroed out in the photodiode trace
ignore_trials = []; % 
A02_parse_events_CWH(preproc_id, preproc_file, log_dir_struct, photo_channel_n,...
    mic_channel_n, rt_window, ignore_trials, plot_it, save_it);


%% Combine blocks
A03_combine_blocks_CWH(SBJ, preproc_id, {preproc_file});
% Check that trial types recovered by event triggers match those in the log file
check_trial_info(SBJ,data_id);

%% Spectral decompose
freq_band = [2 350];        % Range of frequency bands to analyze in Hz
frac_bandwidth = 0.25;      % Bandwith of each frequency.
                            %   Lower values give more precise phase estimation but worse temporal resolution.
                            %   Higher values give less precise phase estimation but better temporal resolution.
                            %   Don't go more than 0.35 or less than 0.15
freq_spacing_ratio = 0.5;  % Each center freq is spaced apart by this number multiplied by the lower freq
                            %   [ex if 0.1: 30 33 36 40 44 48 53 58 64 71 78 ...]
phase_cutoff_freq = 50;     % Don't save phase values for frequencies above this
% Call the function
A04_spectral_decompose_CWH(SBJ, preproc_id, freq_band, freq_spacing_ratio, frac_bandwidth, phase_cutoff_freq);

%% Artifact reject and define epochs
plot_channels = [1 2 3 4];
artifact_struct.std_limit_raw = 5.0;
artifact_struct.std_limit_diff = 6.0;
artifact_struct.hard_threshold_raw = 150;
artifact_struct.hard_threshold_diff = 300;
time_range = [-0.200 2.000];
A05_artifact_reject(preproc_id, time_range, artifact_struct, plot_channels);

% %% Plot ERPS
% n_channels = 1;
% epoch_time = [-0.1 1.0];
% for channel_n = 1:n_channels
%   C01_plot_erp(subject_id, channel_n, epoch_time);
%   pause;
%   close all
% end
% 
% C04_stack_erp_by_rt(subject_id, 1, epoch_time);
% 
% %% Plot spectral erps
% n_channels = 1;
% epoch_time = [-0.1 1.0];
% 
% % HG
% freq_range = [100 250];
% for channel_n = 1:n_channels
%   C02_plot_spectral_erp(subject_id, channel_n, epoch_time, freq_range);
%   pause;
%   close all
% end
% % Theta
% freq_range = [4 6];
% for channel_n = 1:n_channels
%   C02_plot_spectral_erp(subject_id, channel_n, epoch_time, freq_range);
%   pause;
%   close all
% end


