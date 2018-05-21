% Plot data using ft_databrowser.
addpath(genpath('/home/knight/hoycw/Apps/fieldtrip-20160927/'));

SBJ = 'IR41';%input('Enter subject ID:  ','s');
raw_filename = '2016051011_0001.besa';
data_id = 'IR41_LSM';

SBJ_dir = fullfile('/home/knight/hoycw/PRJ_Stroop/data/',SBJ);
raw_dir = fullfile(SBJ_dir, '00_raw_data/');
raw_filename = fullfile(raw_dir,raw_filename);

out_dir = fullfile(SBJ_dir,'01_import/');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end


% Filtering parameters
filter_it   = 0;
notch_freqs = [60 120 180 300]; 

% basic plotting parameters
plot_it = 0;
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
elseif strcmp(raw_filename(end-3:end),'.edf')
    error('.edf file detected, implement kla code to extract');
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

% Filter data for ease of viewing
if filter_it
    cfg_notch = [];
    cfg_notch.continuous = 'yes';
    cfg_notch.dftfilter = 'yes'; % line noise removal using discrete fourier transform
    cfg_notch.dftfreq = notch_freqs;
    cfg_notch.demean = 'yes';
    data_notch = ft_preprocessing(cfg_notch, data_raw); %get to know this guy
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
%events start ~155 or 160s to ~384; ~580 to end (~1360?)
analysis_time = {[140 400],[560 1360]}; % In seconds. Leave at least 10 seconds on edge for edge artifacts

% Focusing on RSM and LSM since they're close to aMCC
% RSM1-3 are GM (esp 1!), use 4 and 5 as ref (4 is questionable GM/WM)
% LSM1-3 are GM (esp 1!), use 4 as ref
possible_data_lab = {'LSM1', 'LSM2', 'LSM3'};
possible_refs_lab = {'LSM4'}; % this is just based on anatomy,
%   actually LSM4 looks very similar and will ikely destroy my ERPs

data_ix = [];
for lab_ix = 1:length(data_raw.label)
    for data_lab_ix = 1:length(possible_data_lab)
        if strcmp(data_raw.label(lab_ix),possible_data_lab(data_lab_ix))
            data_ix = [data_ix lab_ix];
        end
    end
end
wm_refs_ix = [];
for lab_ix = 1:length(data_raw.label)
    for ref_lab_ix = 1:length(possible_refs_lab)
        if strcmp(data_raw.label(lab_ix),possible_refs_lab(ref_lab_ix))
            wm_refs_ix = [wm_refs_ix lab_ix];
        end
    end
end
   
%% Define the rereferencing scheme (WM ref)
%    Each cell is an array of channels to include in one channel of rereferenced data
%     reref_weights assigns a weight to each of the channels for rereferencing. Should add up to 1
reref_channels = {[45 48], [46 48], [47 48]};
refwt = -1*ones(1,length(wm_refs_ix))/length(wm_refs_ix);
reref_weights  = {[1 refwt], [1 refwt], [1 refwt]};

% Analysis channels are all of the channels used in reref_channels
analysis_channels = [];
for chan_n = 1:length(reref_channels)
  analysis_channels = [analysis_channels reref_channels{chan_n}];
end
analysis_channels = unique(analysis_channels);

%% Save files to prepare for KLA scripts
data_ecog = data_raw.trial{:}(analysis_channels,:); %only analysis_channels
%   BUT! analysis_channels variable should be original chan #s
header_ecog.sample_rate = data_raw.fsample;
header_ecog.orig_n_channels = size(data_raw.label,1);
header_ecog.n_channels = size(analysis_channels,2);
header_ecog.n_samples = size(data_ecog,2);
header_ecog.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
header_ecog.original_sample_rate = data_raw.fsample;
header_ecog.channel_labels = {data_raw.label{analysis_channels}};
header_ecog.orig_channel_n = analysis_channels;

data_evnt = data_raw.trial{:}(event_channels,:);
header_evnt.sample_rate = data_raw.fsample;
header_evnt.orig_n_channels = size(data_raw.label,1);
header_evnt.n_channels = size(event_channels,2);
header_evnt.n_samples = size(data_ecog,2);
header_evnt.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
header_evnt.original_sample_rate = data_raw.fsample;
header_evnt.channel_labels = {data_raw.label{event_channels}};
header_evnt.orig_channel_n = event_channels;

data_out_filename = strcat(out_dir,data_id,'.mat');
save(data_out_filename, 'header_ecog', 'data_ecog', 'header_evnt', 'data_evnt');
proc_out_filename = strcat(out_dir,data_id,'_proc_vars.mat');
save(proc_out_filename, 'analysis_time', 'analysis_channels', ...
    'event_channels', 'reref_channels', 'reref_weights', 'photo_channel_n', 'mic_channel_n');

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

