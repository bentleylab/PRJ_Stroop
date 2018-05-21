%% DATA_ID: IR21_RC_WM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run sections individually to process the dataset
close all
clear all

plot_raw_PSDs = 0;

SBJ = 'IR32';
SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
raw_data_id = 'IR32_IH';
input_file = strcat(raw_data_id, '.mat');
input_filename = strcat(SBJ_dir,'01_import/',input_file);
proc_vars_filename = strcat(SBJ_dir,'01_import/',raw_data_id,'_proc_vars.mat');
preproc_id = strcat(raw_data_id,'_CA');
preproc_file = strcat(preproc_id,'.mat');
preproc_filename = strcat(SBJ_dir,'02_preproc/',preproc_file);
% output_dir = strcat(SBJ_dir,'03_proc/'); now specified in A03_comb blocks

% Filtering parameters
hp_cutoff = 0.5;

% Add toolboxes to path
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/'));
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath('/home/knight/hoycw/Apps/fieldtrip-20160927/'); % read .besa
ft_defaults

%% Look at noise profile of all channels
% Load data
load(proc_vars_filename);
load(input_filename);

% Loop through channels to see noise profile
if plot_raw_PSDs
    for channel_n = 1:header_ecog.n_channels
        [fft_data,freqs] = pwelch(data_ecog(channel_n,:),2048,0,2048,header_ecog.sample_rate);
        loglog(freqs,fft_data);
        xlim([1 350]);
        ax = gca;
        ax.XTick = [25 30 60 100 120 150 180 200 240 300 360 420];
        title(['Channel ' num2str(channel_n) ' [' header_ecog.channel_labels{channel_n} ']']);
        pause;
    end
    close all;
end
%% Set line filters based on noise profiles plotted above
% Big bumps at 60, 180, 300 Hz, some tiny at 120 but not really
line_noise_freqs = [60 120 180 240 300];
dummy_line_noise_freqs = []; %FIX ME!!! because cleanline2.m doesn't work

%% Eliminate bad channels
bad_chan = {'IHR20','IHR21','IHR30','IHR31'};
bad_chan_idx = [];
for bc_ix = 1:length(bad_chan)
    bad_chan_idx = [bad_chan_idx find(ismember(header_ecog.channel_labels,bad_chan{bc_ix}),1)];
end
good_chan_idx = setdiff(1:size(data_ecog,1),bad_chan_idx);
data_ecog = data_ecog(good_chan_idx,:);
header_ecog.bad_channel_labels = bad_chan;
header_ecog.bad_channel_n = header_ecog.orig_channel_n(bad_chan_idx);
header_ecog.channel_labels(bad_chan_idx) = [];
header_ecog.orig_channel_n(bad_chan_idx) = [];
header_ecog.n_channels = length(good_chan_idx);
car_exclude_chan = {'IHR18','IHR27','IHR28'};
bad_chan_idx = [];

%% Importing all data blocks
% I'm implementing common average referencing myself here
% A01_import_reref_filter_ft_besa_good_time(preproc_id, input_file, preproc_file, analysis_channels, event_channels, ...
%   'SBJ_dir_name', SBJ_dir, ...
%   'resample_rate', 1000, ...
%   'reref_channels', reref_channels, ...
%   'reref_name', 'WM', ...
%   'reref_weights', reref_weights, ...
%   'line_noise_freqs', dummy_line_noise_freqs, ...
%   'analysis_time', analysis_time, ...
%   'highpass_freq', hp_cutoff, ...
%   'save_output', 1 ...
%   );
orig_analysis_time = analysis_time;

%% Resample
resample_rate = 1000;
header_ecog.sample_rate = resample_rate;
if(header_ecog.sample_rate ~= header_ecog.original_sample_rate)
  fprintf('Resampling from %dHz to %dHz\n',header_ecog.original_sample_rate,header_ecog.sample_rate);
  resample_data = zeros(header_ecog.n_channels, ...
      ceil(header_ecog.n_samples*header_ecog.sample_rate/header_ecog.original_sample_rate));
  fprintf('.');
  for channel_n = 1:header_ecog.n_channels
    resample_data(channel_n,:) = resample(data_ecog(channel_n,:),header_ecog.sample_rate,...
        header_ecog.original_sample_rate,100); % 100 increases accuracy for resample()
    fprintf('.');
  end
  fprintf('\n');
  data_ecog = resample_data;
  clear resample_data;
  header_ecog.n_samples = size(data_ecog,2);
end

%% High Pass Filter
% fprintf('High pass filtering at %2.2fHz\n',hp_cutoff);
% %   First remove the temporal mean (detrend)
% data_ecog = detrend(data_ecog','constant')'; % This assumes data is [channels x samples]
% %   Then apply highpass filter
% hpfilt_order = ceil((header_ecog.sample_rate/hp_cutoff)*5);
% hpfilt_cutoff_norm = 2*hp_cutoff/header_ecog.sample_rate;
% hpfilt_coeff = fir1(hpfilt_order, hpfilt_cutoff_norm, 'high');
% hpfilt_group_delay = median(grpdelay(hpfilt_coeff));
% %   Filter the data (causal) and then shift signal to make it zero-phase
% data_ecog = filter(hpfilt_coeff,1,[data_ecog zeros(header_ecog.n_channels,hpfilt_group_delay)]')';  % This assumes data is [channels x samples]
% data_ecog = data_ecog(:,(hpfilt_group_delay+1):end);
% %data_ecog = X_eegfilt(data_ecog,header_ecog.sample_rate,hp_cutoff,[]); % Use EEGLAB (filtfilt which is slower)
% header_ecog.highpass = hp_cutoff;

%% Trim data to analysis_time
if(~isempty(analysis_time{1}))
    all_analysis_points_ecog = [];
    all_analysis_points_evnt = [];
    for time_grp = 1:length(analysis_time)
        all_analysis_points_ecog = [all_analysis_points_ecog ...
            (analysis_time{time_grp}(1)*header_ecog.sample_rate):(analysis_time{time_grp}(2)*header_ecog.sample_rate)];
%             analysis_limits{time_grp}(1):analysis_limits{time_grp}(2)];
        all_analysis_points_evnt = [all_analysis_points_evnt ...
            (analysis_time{time_grp}(1)*header_evnt.sample_rate):(analysis_time{time_grp}(2)*header_evnt.sample_rate)];
    end
    if sum(all_analysis_points_ecog<1) | sum(all_analysis_points_evnt<1)
        disp('===================================================');
        fprintf('WARNING: %i ecog points less than 0\n',sum(all_analysis_points_ecog<1));
        fprintf('WARNING: %i event points less than 0\n',sum(all_analysis_points_evnt<1));
        disp('===================================================');
    end
    all_analysis_points_ecog(all_analysis_points_ecog<1) = [];
    all_analysis_points_evnt(all_analysis_points_evnt<1) = [];
    if sum(all_analysis_points_ecog>header_ecog.n_samples) | sum(all_analysis_points_evnt>header_evnt.n_samples)
        disp('===================================================');
        fprintf('WARNING: %i ecog points longer than n_samples\n',sum(all_analysis_points_ecog>header_ecog.n_samples));
        fprintf('WARNING: %i event points longer than n_samples\n',sum(all_analysis_points_evnt>header_evnt.n_samples));
        disp('===================================================');
    end
    all_analysis_points_ecog(all_analysis_points_ecog>header_ecog.n_samples) = [];
    all_analysis_points_evnt(all_analysis_points_evnt>header_evnt.n_samples) = [];
    data_ecog = data_ecog(:,all_analysis_points_ecog);
    data_evnt = data_evnt(:,all_analysis_points_evnt);
    header_ecog.n_samples = size(data_ecog,2);
    header_evnt.n_samples = size(data_evnt,2);
    header_ecog.length_in_seconds = header_ecog.n_samples / header_ecog.sample_rate;
    header_evnt.length_in_seconds = header_evnt.n_samples / header_evnt.sample_rate;
end
header_ecog.analysis_time = orig_analysis_time;
header_evnt.analysis_time = orig_analysis_time;

%% Rereference
ref_ix = 1:size(data_ecog,1);
% Get rid of CAR Excluded Channels
%   IHR18, which has a big spike around 750s
%   IHR27,28 for some suspicious entrainment with epileptic rhythms (potential)
for lab_ix = 1:length(car_exclude_chan)
    ref_ix(find(ismember(header_ecog.channel_labels,car_exclude_chan{lab_ix}),1)) = [];
end
CAR = mean(data_ecog(ref_ix,:),1);

header_ecog.orig_channel_labels = header_ecog.channel_labels; % Save the original channel labels
header_ecog.channel_labels = [];
for ch_n = 1:size(data_ecog,1)
    data_ecog(ch_n,:) = data_ecog(ch_n,:) - CAR;
    header_ecog.channel_labels{ch_n} = [header_ecog.orig_channel_labels{ch_n} '-CA'];
end
header_ecog.reref_channels = header_ecog.orig_channel_n(ref_ix);
header_ecog.reref_weights = 'CAR-IHR18,27,28';
header_ecog.data_id = preproc_id;
header_evnt.data_id = preproc_id;
data_id = preproc_id;

%% use ft to remove line noise
data_notch = [];
data_notch.label = header_ecog.channel_labels;
data_notch.fsample = header_ecog.sample_rate;
data_notch.trial{1} = data_ecog;
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
plot(linspace(0,header_evnt.length_in_seconds,header_evnt.n_samples), data_evnt);
% 9 blocks, but nurse interrupts B6 after T4, then picks up and finishes block

% Manually correct bad event traces
% Fix IR32: 
% noise at front end, B1T1 might have some trouble due to noise at onset
% blip at end of B1 as usual
% B6T4 and T5 for interruption
yval_bsln = 633; %this is the zero/baseline during a block
bsln_set_times = {[0.0 13.0], [112.0 113.0], [707.5 787.0]}; % in sec
data_evnt(photo_channel_n,1:floor(bsln_set_times{1}(2)*header_evnt.sample_rate)) = yval_bsln;
data_evnt(photo_channel_n,floor(bsln_set_times{2}(1)*header_evnt.sample_rate):floor(bsln_set_times{2}(2)*header_evnt.sample_rate)) = yval_bsln;
data_evnt(photo_channel_n,floor(bsln_set_times{3}(1)*header_evnt.sample_rate):floor(bsln_set_times{3}(2)*header_evnt.sample_rate)) = yval_bsln;
save(preproc_filename, 'header_ecog', 'data_ecog', 'header_evnt', 'data_evnt', 'data_id');


%% Load events
% mic_channel_n = 2;
rt_window = [0.300 2.000]; % Possible response times from word onset in seconds
plot_it = 1; % Set to 0 if you don't want to see plots
save_it = 1; % 0 if you don't want to write files
log_dir_struct.log_file_str = strcat(SBJ,'_strooptask_log.txt');
log_dir_struct.wav_file_dir = strcat(SBJ,'_strooptask_wavfiles');
% ingore_trials is for the log file; those should already be zeroed out in the photodiode trace
ignore_trials = [220 221];
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
plot_channels = [1 2];
artifact_struct.std_limit_raw = 5.0;
artifact_struct.std_limit_diff = 6.0;
artifact_struct.hard_threshold_raw = 250;
artifact_struct.hard_threshold_diff = 300;
time_range = [-0.200 2.000]; 
A05_artifact_reject_CWH(preproc_id, time_range, artifact_struct, plot_channels);

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


