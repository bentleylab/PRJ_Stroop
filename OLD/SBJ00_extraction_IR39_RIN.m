%% Extract all raw data with correct names
% This script should be run interactively, section by section
% Create a new version of this script for each subject to record important
% parameters and events
%   1) Load .edf or .besa format data
%   2) Get DC channel info, plot to get a feeling for the blocks and analysis times
%   3) Plot data to get parameters for automatic cleaning
%   4) Set up variable for preprocessing, save as .mat
%
% Differences from original A00_prep_raw_data:
%   - no referencing (left for separate script)
%   - saves out all data, so subsequent scripts don't need extraction
clear all; close all;

addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath('/home/knight/hoycw/Apps/fieldtrip-20160927/');
ft_defaults

SBJ = 'IR39';%input('Enter subject ID:  ','s');
raw_filename = '2016042320_0029.besa';
raw_data_id = 'IR39_RIN';

SBJ_dir = fullfile('/home/knight/hoycw/PRJ_Stroop/data/',SBJ);
raw_dir = fullfile(SBJ_dir, '00_raw/');
raw_filename = fullfile(raw_dir,raw_filename);

out_dir = fullfile(SBJ_dir,'01_import/');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

% Filtering parameters --> FOR PLOTTING PURPOSES ONLY! Not saved this way
filter_it   = 1;
notch_freqs = [60 120 180 240 300]; 

% basic plotting parameters
plot_it = 1;
cfg_plot = [];
cfg_plot.viewmode = 'vertical'; %trace per line
cfg_plot.continuous = 'yes'; % not trials
cfg_plot.plotlabels = 'yes'; %to know what to cut
cfg_plot.blocksize = 10; % good for bob, can adjust

%% Load the data
if strcmp(raw_filename(end-4:end),'.besa')
    cfg = [];
    cfg.dataset = raw_filename;
%     cfg.trialfun = 'trialfun_besa';
%     cfg = ft_definetrial(cfg);
    cfg.continuous = 'yes';
    cfg.channel = 'all';
    data_raw = ft_preprocessing(cfg);   % just load the data, don't process it
    % Define just the three fields needed for this script
    data_ecog = data_raw.trial{1};
    header_ecog.n_samples = size(data_ecog,2);
    header_ecog.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
    % Remove extra characters from channel labels
    for channel_n = 1:length(data_raw.label)
        data_raw.label{channel_n} = strrep(data_raw.label{channel_n},'POL','');
        data_raw.label{channel_n} = strrep(data_raw.label{channel_n},'Ref','');
    end
elseif strcmp(raw_filename(end-3:end),'.edf')
    [header_ecog, data_ecog] = extract_nk_converted_edf(raw_filename, [], []);
    % Remove extra characters from channel labels
    for channel_n = 1:header_ecog.n_channels
        header_ecog.channel_labels{channel_n} = strrep(header_ecog.channel_labels{channel_n},'POL','');
        header_ecog.channel_labels{channel_n} = strrep(header_ecog.channel_labels{channel_n},'Ref','');
    end
    % Put into ft format for plotting
    data_raw = [];
    data_raw.label = header_ecog.channel_labels;
    data_raw.fsample = header_ecog.sample_rate;
% NOTE: other functions will assume your data is in the trial field
    data_raw.trial{1} = data_ecog; %data appears here by default before you do anything
    data_raw.time{1} = linspace(0, header_ecog.length_in_seconds, header_ecog.n_samples); % ms time spacing
% elseif strcmp(raw_filename(end-3:end),'.mat')
%     load(raw_filename);
%     data_raw = [];
%     data_raw.label = elec_labels;
%     data_raw.fsample = 1000; % WARNING: might be false! shoudl automate detection
%NOTE: other functions will assume your data is in the trial field
%     data_raw.trial{1} = gdat_clean_filt; %data appears here by default before you do anything
%     data_raw.time{1} = linspace(0, size(gdat_clean_filt,2), size(gdat_clean_filt,2)); % ms time spacing
else
    error(strcat('Unknown raw data format: ',raw_filename));
end

%% Filter data for ease of viewing
if filter_it
    cfg_notch = [];
    cfg_notch.continuous = 'yes';
    cfg_notch.dftfilter = 'yes'; % line noise removal using discrete fourier transform
    cfg_notch.dftfreq = notch_freqs;
    cfg_notch.demean = 'yes';
    data_notch = ft_preprocessing(cfg_notch, data_raw);
else
    data_notch = data_raw;
end

%% Plot data: check good channels, refs, artifacts, etc.
if plot_it
    ft_databrowser(cfg_plot, data_notch); %pulls up your plotted data
end

% QUESTIONS:
%    What should the thresholds be for artifact rejection? Any channels to throw out?
% Hard cutoffs should be -150 to 150. Hard cutoffs are data that is thrown out before calculating
%    standard deviation for throwing out atrifactual epochs.


%% Create variable for KLA scripts
event_channels = [1 2]; % From all imported channels
photo_channel_n = 1; % Relative to event_channels
mic_channel_n = 2;   % Relative to event_channels

% Set valid times to use. Take out long pauses with no events.
% In seconds. Leave at least 10 seconds on edge for edge artifacts
% events start ~155 or 160s to ~384; ~580 to end (~1360?)
analysis_time = {[135 380], [535 1370]};

%% Data and reference selection
% !!! Cut this out!
% % Focusing on RAC (good aMCC)
% %   RAC2-6 area ll good GM, RAC4 looks like great aMCC
% %   RAC1 is good WM ref on L hemisphere, RAC7 is good on R (more superior)
% possible_data_lab = {'RAC2', 'RAC3', 'RAC4', 'RAC5', 'RAC6'};
% possible_refs_lab = {'RAC1', 'RAC7'}; %!!! careful with this weird surrounding ref...
% 
% data_ix = [];
% for lab_ix = 1:length(data_raw.label)
%     for data_lab_ix = 1:length(possible_data_lab)
%         if strcmp(data_raw.label(lab_ix),possible_data_lab(data_lab_ix))
%             data_ix = [data_ix lab_ix];
%         end
%     end
% end
% wm_refs_ix = [];
% for lab_ix = 1:length(data_raw.label)
%     for ref_lab_ix = 1:length(possible_refs_lab)
%         if strcmp(data_raw.label(lab_ix),possible_refs_lab(ref_lab_ix))
%             wm_refs_ix = [wm_refs_ix lab_ix];
%         end
%     end
% end
%    
% %% Define the rereferencing scheme (WM ref)
% % Plot the potential WM refs
% for channel_n = 1:length(wm_refs_ix)
%     figure%('Position',[100 500 1400 400]);
%     for channel_n2 = 1:length(wm_refs_ix)
%         plot(linspace(0,header_ecog.length_in_seconds,header_ecog.n_samples),...
%             data_ecog(wm_refs_ix(channel_n2),:), 'Color', [0.7 0.7 0.7]); hold on;
%     end
%     plot(linspace(0,header_ecog.length_in_seconds,header_ecog.n_samples),...
%         data_ecog(wm_refs_ix(channel_n),:),'r');
%     title(strcat(data_raw.label{wm_refs_ix(channel_n)},'-',num2str(wm_refs_ix(channel_n))));
% %     ylim([-1500 1500]);
%     pause;
%     close all;
% end
% 
% % Define the rereference weightings
% %    Each cell is an array of channels to include in one channel of rereferenced data
% %     reref_weights assigns a weight to each of the channels for rereferencing. Should add up to 1
% reref_channels = {};
% reref_weights = {};
% refwt = -1*ones(1,length(wm_refs_ix))/length(wm_refs_ix);
% for ch_n = 1:length(data_ix)
%     reref_channels = {reref_channels{:}, [data_ix(ch_n) wm_refs_ix]};
%     reref_weights  = {reref_weights{:}, [1 refwt]};
% end

% Analysis channels are all of the channels used in reref_channels
reref_channels = {};
analysis_channels = [];
for chan_n = 1:length(reref_channels)
  analysis_channels = [analysis_channels reref_channels{chan_n}];
end
analysis_channels = unique(analysis_channels);

%% Set up files to prepare for KLA scripts
if strcmp(raw_filename(end-4:end),'.besa')
    data_evnt = data_raw.trial{:}(event_channels,:);
    header_evnt.sample_rate = data_raw.fsample;
    header_evnt.orig_n_channels = size(data_raw.label,1);
    header_evnt.n_samples = size(data_ecog,2);
    header_evnt.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
    header_evnt.original_sample_rate = data_raw.fsample;
    header_evnt.channel_labels = {data_raw.label{event_channels}};
elseif strcmp(raw_filename(end-3:end),'.edf')
    data_evnt = data_ecog(event_channels,:);
    header_evnt = header_ecog;
    header_evnt.channel_labels = header_ecog.channel_labels(event_channels);
end
header_evnt.n_channels = size(event_channels,2);
header_evnt.orig_channel_n = event_channels;
header_evnt.raw_file = raw_filename;

if strcmp(raw_filename(end-4:end),'.besa')
    data_ecog = data_raw.trial{:}(analysis_channels,:); %only analysis_channels
    %   BUT! analysis_channels variable should be original chan #s
    header_ecog.sample_rate = data_raw.fsample;
    header_ecog.orig_n_channels = size(data_raw.label,1);
    header_ecog.n_samples = size(data_ecog,2);
    header_ecog.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
    header_ecog.original_sample_rate = data_raw.fsample;
    header_ecog.channel_labels = {data_raw.label{analysis_channels}};
elseif strcmp(raw_filename(end-3:end),'.edf')
    data_ecog = data_ecog(analysis_channels,:);
    header_ecog.channel_labels = header_ecog.channel_labels(analysis_channels);
end
header_ecog.n_channels = size(analysis_channels,2);
header_ecog.orig_channel_n = analysis_channels;
header_ecog.raw_file = raw_filename;

% Grab microphone recording separately to export and listen
mic_data = data_evnt(mic_channel_n,:);
%rescale to prevent clipping, add 0.05 fudge factor
mic_data_rescale = mic_data./(max(abs(mic_data))+0.05); 

%% Save out those files
data_out_filename = strcat(out_dir,raw_data_id,'.mat');
evnt_out_filename = strcat(out_dir,SBJ,'_events.mat');
proc_out_filename = strcat(out_dir,raw_data_id,'_proc_vars.mat');
mic_data_filename = strcat(out_dir,raw_data_id,'_mic_recording.wav');

save(data_out_filename, 'header_ecog', 'data_ecog');
save(evnt_out_filename, 'header_evnt', 'data_evnt');
save(proc_out_filename, 'analysis_time', 'analysis_channels', ...
    'event_channels', 'reref_channels', 'reref_weights', 'photo_channel_n', 'mic_channel_n');
audiowrite(mic_data_filename,mic_data_rescale,header_evnt.sample_rate);



% Call to extract_nk_...:
% subject_id, input_file_name, save_file_name, analysis_channels, event_channels, ...
%   'data_dir_name', ['../' subject_id], ...
%   'resample_rate', 1000, ...
%   'reref_channels', reref_channels, ...
%   'reref_name', 'WM', ...
%   'reref_weights', reref_weights, ...
%   'line_noise_freqs', line_noise_freqs, ...
%   'analysis_time', analysis_time, ...
%   'highpass_freq', hp_cutoff ...

% Contents of header that extract_nk_... produces:
%    header.recording_startdate  - [string], DD.MM.YY of original recording
%    header.recording_starttime  - [string], HH.MM.SS of original recording
%    header.orig_n_channels      - [integer], The number of channels originally recorded
%    header.n_channels           - [integer], The number of channels currently extracted
%    header.n_samples            - [integer], The number of samples per channel
%    header.length_in_seconds    - [float], The length of the extracted data in seconds
%    header.sample_rate          - [integer], The current sampling rate in Hz.
%    header.original_sample_rate - [integer], The recording sampling rate in Hz. Usually 5000.
%    header.channel_labels       - [cell array of strings], The labels of each extracted channel.
%    header.orig_channel_n       - [array of integers], The corresponding channel numbers of the
%                                  extracted data in the original data.

