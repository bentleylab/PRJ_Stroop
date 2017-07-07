function A001_run_IR39WM_RAC(no_interaction)
% Set no_interaction to 1 to run script without user intervention
% KLA email:
%   The way to do it for a new patient is:
%   0 - Look at scans to determine which electrodes you want to use (for reference and for analysis)
%   1 - make a copy of the A00_run_S04IR26WM.m file and rename to the new patient.
%   2 - Don't run the whole script, just select sections and replace information as you go. I don't have anything in here to read besa files but there's a function in fieldtrip to do it that you could add in.

%% SUBJECT ID: IR39WM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SBJ = 'IR39';
SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
data_id = 'IR39_RAC';
input_file = strcat(data_id, '.mat');
input_filename = strcat(SBJ_dir,'01_import/',input_file);
proc_vars_filename = strcat(SBJ_dir,'01_import/',data_id,'_proc_vars.mat');
preproc_file = strcat(data_id,'_WM.mat');
output_dir = strcat(SBJ_dir,'03_prcs_data/');

% Filtering parameters
hp_cutoff = 0.5;

% Add toolboxes to path
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath(genpath('/home/knight/hoycw/Apps/fieldtrip-20160927/')); % read .besa


%% Look at noise profile of all channels
if(~no_interaction)
    % Load data
    load(proc_vars_filename);
    load(input_filename);

    % Loop through channels to see noise profile
    for channel_n = 1:header_ecog.n_channels
      [fft_data,freqs] = pwelch(data_ecog(channel_n,:),2048,0,2048,header_ecog.sample_rate);
      loglog(freqs,fft_data);
      xlim([1 350]);
      ax = gca;
      ax.XTick = [10 20 25 30 60 100 120 180 200 240 300 360 420];
      title(['Channel ' num2str(channel_n) ' [' header_ecog.channel_labels{channel_n} ']']);
      pause;
    end
end
close all;
% Again, low peak this time at ~20-25 (???!), but super clean (nothing at 60, 120, etc.),
% peaks at 180 (small) and clearly 300 Hz...
ft_line_noise_freqs = [60 120 180 300 ];

%% Importing all data blocks

% Re-import the data using the information we just noted
%  Also highpass filter which is slow
%   !!!CWH: will detrend (demean, only linear) if no high pass; higher polynomial fits???
line_noise_freqs = [];
A01_import_reref_filter_ft_besa(data_id, input_file, preproc_file, analysis_channels, event_channels, ...
  'SBJ_dir_name', SBJ_dir, ...
  'resample_rate', 1000, ...
  'reref_channels', reref_channels, ...
  'reref_name', 'WM', ...
  'reref_weights', reref_weights, ...
  'line_noise_freqs', line_noise_freqs, ...
  'analysis_time', analysis_time, ...
  'highpass_freq', hp_cutoff, ...
  'save_output', 1 ...
  );

%% use ft to remove line noise since cleanline2 failed...
clear header_ecog data_ecog
preproc_filename = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/02_preproc/',preproc_file);
load(preproc_filename);

% resample if you didn't run A01
% data_ecog = ft_preproc_resample(data_ecog, header_ecog.sample_rate, 1000, 'resample');
% header_ecog.sample_rate = 1000;

data_clean = [];
data_clean.label = header_ecog.channel_labels;
data_clean.fsample = header_ecog.sample_rate;
% NOTE: other functions will assume your data is in the trial field
data_clean.trial{1} = data_ecog; % data appears here by default before you do anything
data_clean.time{1} = linspace(0, size(data_ecog,2), size(data_ecog,2)); % ms time spacing

cfg = [];
cfg.dftfilter = 'yes';
cfg.dftfreq = ft_line_noise_freqs;
data_clean = ft_preprocessing(cfg, data_clean);

% Convert back to mat...
data_ecog = data_clean.trial{:};
header_ecog.line_noise_freqs = ft_line_noise_freqs;

%% Check event channels
if(~no_interaction)
    % Look at event channels for Block 1
    plot(linspace(0,header_evnt.length_in_seconds,header_evnt.n_samples), data_evnt);
    % The photodiode data between Block 1 and Block 2 is bad
end
% Manually correct bad event traces
% this means find places the photodiode is going crazy and set those tothe
% constant baseline (yval)
% A little blip after the first block. We need to fix that.
% Fix the space between Block 1 and Block 2
% photo_channel_n = 1;
yval = 400;
set_times = [0.0 150.0 1362.0 1370.0];
data_evnt(photo_channel_n,1:floor(set_times(2)*header_evnt.sample_rate)) = yval;
data_evnt(photo_channel_n,floor(set_times(3)*header_evnt.sample_rate):floor(set_times(4)*header_evnt.sample_rate)) = yval;
save(preproc_filename, 'header_ecog', 'data_ecog', 'header_evnt', 'data_evnt', 'data_id');


%% Load events
% mic_channel_n = 2;
rt_window = [0.300 2.000]; % Possible response times from word onset in seconds
plot_it = 1; % Set to 0 if you don't want to see plots
if(no_interaction)
  plot_it = 0; % Set to 0 if you don't want to see plots
end
log_dir_struct.log_file_str = 'IR39_strooptask_log.txt';
log_dir_struct.wav_file_dir = 'IR39_strooptask_wavfiles';
ignore_trials = [];
A02_parse_events_CWHedits(data_id, preproc_file, log_dir_struct, photo_channel_n, mic_channel_n, rt_window, ignore_trials, plot_it);


%% Combine blocks
A03_combine_blocks(data_id, {[data_id '.mat']});


%% Spectral decompose
freq_band = [2 350];        % Range of frequency bands to analyze in Hz
frac_bandwidth = 0.25;      % Bandwith of each frequency.
                            %   Lower values give more precise phase estimation but worse temporal resolution.
                            %   Higher values give less precise phase estimation but better temporal resolution.
                            %   Don't go more than 0.35 or less than 0.15
freq_spacing_ratio = 0.03;  % Each center freq is spaced apart by this number multiplied by the lower freq
                            %   [ex if 0.1: 30 33 36 40 44 48 53 58 64 71 78 ...]
phase_cutoff_freq = 50;     % Don't save phase values for frequencies above this
% Call the function
A04_spectral_decompose(data_id, freq_band, freq_spacing_ratio, frac_bandwidth, phase_cutoff_freq);

%% Artifact reject and define epochs
if(~no_interaction)
  plot_channels = [1 2 3 4];
else
  plot_channels = [];
end
artifact_struct.std_limit_raw = 5.0;
artifact_struct.std_limit_diff = 6.0;
artifact_struct.hard_threshold_raw = 150;
artifact_struct.hard_threshold_diff = 300;
time_range = [-0.100 1.000];
A05_artifact_reject(data_id, time_range, artifact_struct, plot_channels);

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


