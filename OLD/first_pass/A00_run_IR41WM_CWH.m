function A00_run_IR41WM_CWH(no_interaction)
% Set no_interaction to 1 to run script without user intervention
% KLA email:
%   The way to do it for a new patient is:
%   0 - Look at scans to determine which electrodes you want to use (for reference and for analysis)
%   1 - make a copy of the A00_run_S04IR26WM.m file and rename to the new patient.
%   2 - Don't run the whole script, just select sections and replace information as you go. I don't have anything in here to read besa files but there's a function in fieldtrip to do it that you could add in.

%% SUBJECT ID: S04IR26WM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subject_id = 'IR41WM';
file_name = '2016051011_0001.besa';

% Add toolboxes to path
addpath(genpath(fullfile(pwd,'_TOOLBOXES')));
addpath(genpath('/home/knight/hoycw/Apps/fieldtrip-20160927/')); % read .besa

%% Print channel names and look at event channels
if(~no_interaction)
    
    % Print channel correspondence
%     print_channel_names(['./../' subject_id '/_EDF/' file_name]);
    cfg = [];
    cfg.dataset = raw_filename;
%     cfg.trialfun = 'trialfun_besa';
%     cfg = ft_definetrial(cfg);
    cfg.continuous = 'yes';
    cfg.channel = 'all';
    data_raw = ft_preprocessing(cfg);   % just load the data, don't process it
    
    % Load data and plot event channels
%     [header, data] = extract_nk_converted_edf(['./../' subject_id '/_EDF/' file_name], ...
%         [1:4], []);
%     figure('Position',[100 100 1200 1000])
%     subplot(4,1,1);
%     plot(linspace(0,header.length_in_seconds,header.n_samples),data(1,:)); title('1');
%     subplot(4,1,2);
%     plot(linspace(0,header.length_in_seconds,header.n_samples),data(2,:)); title('2');
%     subplot(4,1,3);
%     plot(linspace(0,header.length_in_seconds,header.n_samples),data(3,:)); title('3');
%     subplot(4,1,4);
%     plot(linspace(0,header.length_in_seconds,header.n_samples),data(4,:)); title('4');
end
% File is 9 blocks
%  Two other files, one with 1 block and the other with nothing

% Set channels
event_channels = [1 2]; % From all imported channels
photo_channel_n = 1; % Relative to event_channels
mic_channel_n = 2;   % Relative to event_channels

% Set valid times to use. Take out long pauses with no events.
%events start ~155 or 160s to ~384; ~580 to end (~1360?)
data_times = {[140 400],[560 1360]}; % In seconds. Leave at least 10 seconds on edge for edge artifacts

% Focusing on RSM and LSM since they're close to aMCC
% RSM1-3 are GM (esp 1!), use 4 and 5 as ref (4 is questionable GM/WM)
%   
% LSM1-3 are GM (esp 1!), use 4 as ref
possible_data_lab = {'LSM1', 'LSM2', 'LSM3'};
possible_data = [];
for lab_ix = 1:length(data_raw.label)
    for data_lab_ix = 1:length(possible_data_lab)
        if strcmp(data_raw.label(lab_ix),possible_data_lab(data_lab_ix))
            possible_data = [possible_data lab_ix];
        end
    end
end
possible_refs_lab = {'LSM4'};
possible_refs = [];
for lab_ix = 1:length(data_raw.label)
    for ref_lab_ix = 1:length(possible_data_lab)
        if strcmp(data_raw.label(lab_ix),possible_data_lab(ref_lab_ix))
            possible_refs = [possible_refs lab_ix];
        end
    end
end
   
%% Look at electrode locations and choose white matter ref electrodes
%    Look for clean refs
if(~no_interaction)
    possible_refs = [112:115]; % Got these from looking at the scans
    [header, data] = extract_nk_converted_edf(['./../' subject_id '/_EDF/' file_name], ...
        [possible_refs], [], 1000);
    for channel_n = 1:length(possible_refs)
        figure('Position',[100 500 1400 400]);
        for channel_n2 = 1:length(possible_refs)
            plot(linspace(0,header.length_in_seconds,header.n_samples),data(channel_n2,:), 'Color', [0.7 0.7 0.7]); hold on;
        end
        plot(linspace(0,header.length_in_seconds,header.n_samples),data(channel_n,:),'r');
        title(num2str(possible_refs(channel_n)));
        ylim([-1500 1500]);
        pause;
        close all;
    end
end
% Just choosing all 4 as WM because I don't have anatomy data yet.
wmref = [112:115];


%% Look at data channels
%    What should the thresholds be for artifact rejection? Any channels to throw out?
if(~no_interaction)
    possible_data = [110:117]; % Got these from looking at the scans
    [header, data] = extract_nk_converted_edf(['./../' subject_id '/_EDF/' file_name], ...
        [possible_data wmref], [], 1000);
    refdata = mean(data((length(possible_data)+1):end,:),1);
    for channel_n = 1:size(data,1)
        data(channel_n,:) = data(channel_n,:) - refdata;
    end
    for channel_n = 1:length(possible_data)
        figure('Position',[100 500 1400 400]);
        for channel_n2 = 1:length(possible_data)
            plot(linspace(0,header.length_in_seconds,header.n_samples),data(channel_n2,:), 'Color', [0.7 0.7 0.7]); hold on;
        end
        plot(linspace(0,header.length_in_seconds,header.n_samples),data(channel_n,:),'r');
        title(num2str(possible_data(channel_n)));
        ylim([-250 250]);
        pause;
        close all;
    end
end
% Looks okay
% Hard cutoffs should be -150 to 150. Hard cutoffs are data that is thrown out before calculating
%    standard deviation for throwing out atrifactual epochs.


%% Define the rereferencing scheme (WM ref)
%    Each cell is an array of channels to include in one channel of rereferenced data
%     reref_weights assigns a weight to each of the channels for rereferencing. Should add up to 1
reref_channels = {[110 wmref], [111 wmref], [112 wmref], [113 wmref], [114 wmref], [115 wmref], [116 wmref], [117 wmref] ... % LPC1
                 };
refwt = -1*ones(1,length(wmref))/length(wmref);
reref_weights  = {[1 refwt], [1 refwt], [1 refwt], [1 refwt], [1 refwt], [1 refwt], [1 refwt], [1 refwt] ...
                 };

% Analysis channels are all of the channels used in reref_channels
analysis_channels = [];
for chan_n = 1:length(reref_channels)
  analysis_channels = [analysis_channels reref_channels{chan_n}];
end
analysis_channels = unique(analysis_channels);


%% Look at noise profile of all channels
if(~no_interaction)
    % Load data
    save_file_name = [subject_id '.mat'];
    input_file_name = file_name;
    A01_import_reref_filter(subject_id, input_file_name, save_file_name, analysis_channels, event_channels, ...
      'data_dir_name', ['../' subject_id], ...
      'resample_rate', 1000, ...
      'reref_channels', reref_channels, ...
      'reref_weights', reref_weights ...
      );

    % Loop through channels to see noise profile
    load(['../' subject_id '/01_import/' save_file_name]);
    for channel_n = 1:header_ecog.n_channels
      [fft_data,freqs] = pwelch(data_ecog(channel_n,:),2048,0,2048,header_ecog.sample_rate);
      loglog(freqs,fft_data);
      xlim([1 350]);
      ax = gca;
      ax.XTick = [60 100 120 180 200 240 300 360 420];
      title(['Channel ' num2str(channel_n) ' [' header_ecog.channel_labels{channel_n} ']']);
      pause;
    end
end
close all;
% Looks like there is noise in 60, 120, 180, 300 Hz
line_noise_freqs = [60 120 180 300 ];



%% Importing all data blocks

% Re-import the data using the information we just noted
%  Also highpass filter which is slow
%   !!!CWH: will detrend (demean, only linear) if no high pass; higher polynomial fits???

save_file_name = [subject_id '.mat'];
input_file_name = file_name;
analysis_time = data_times;
hp_cutoff = 0.5;
A01_import_reref_filter(subject_id, input_file_name, save_file_name, analysis_channels, event_channels, ...
  'data_dir_name', ['../' subject_id], ...
  'resample_rate', 1000, ...
  'reref_channels', reref_channels, ...
  'reref_name', 'WM', ...
  'reref_weights', reref_weights, ...
  'line_noise_freqs', line_noise_freqs, ...
  'analysis_time', analysis_time, ...
  'highpass_freq', hp_cutoff ...
  );


%% Check event channels
if(~no_interaction)
    % Look at event channels for Block 1
    save_file_name = [subject_id '.mat'];
    load(['../' subject_id '/01_import/' save_file_name]);
    plot(linspace(0,header_evnt.length_in_seconds,header_evnt.n_samples), data_evnt);
    % The photodiode data between Block 1 and Block 2 is bad
end
% Manually correct bad event traces
% A little blip after the first block. We need to fix that.
% Fix the space between Block 1 and Block 2
yval = 10990;
set_times = [0.1 10.2 113.0 115.0];
save_file_name = [subject_id '.mat'];
load(['../' subject_id '/01_import/' save_file_name]);
data_evnt(photo_channel_n,floor(set_times(1)*header_evnt.sample_rate):floor(set_times(2)*header_evnt.sample_rate)) = yval;
data_evnt(photo_channel_n,floor(set_times(3)*header_evnt.sample_rate):floor(set_times(4)*header_evnt.sample_rate)) = yval;
save(['../' subject_id '/01_import/' save_file_name], 'header_ecog', 'data_ecog', 'header_evnt', 'data_evnt', 'subject_id');


%% Load events
rt_window = [0.300 2.000]; % Possible response times from word onset in seconds
plot_it = 1; % Set to 0 if you don't want to see plots
if(no_interaction)
  plot_it = 0; % Set to 0 if you don't want to see plots
end
save_file_name = [subject_id '.mat'];
log_dir_struct.log_file_str = '826_strooptask_log.txt';
log_dir_struct.log_file_dir = '2015_Aug_20_2042_826_strooptask';
log_dir_struct.wav_file_dir = '826_strooptask_wavfiles';
ignore_trials = [];
A02_parse_events(subject_id, save_file_name, log_dir_struct, photo_channel_n, mic_channel_n, rt_window, ignore_trials, plot_it);


%% Combine blocks
A03_combine_blocks(subject_id, {[subject_id '.mat']});


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
A04_spectral_decompose(subject_id, freq_band, freq_spacing_ratio, frac_bandwidth, phase_cutoff_freq);

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
A05_artifact_reject(subject_id, time_range, artifact_struct, plot_channels);

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


