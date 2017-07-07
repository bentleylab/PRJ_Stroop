%% Plot preprocessed data for cleaning with Bob
%   1) Load preproc data (downsample, highpass, trim, notch)
%   2) Plot data to identify bad channels and epochs
clear all; %close all;

addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath('/home/knight/hoycw/Apps/fieldtrip-20160927/');
ft_defaults

SBJ_files = {...
    {'IR32','2015121512_0002.edf'},...
    {'ST46','E15-473_0022.edf'},...
    {'IR27','2015082615_0003.edf'},...
    {'IR37','2016031117_0051.besa'}...
%     {'IR26','2015082020_0002.edf vs. 5082020_0002.edf???'},...
%     {'IR21','2015063014_0002.edf'},...
%     {'IR31','2015111117_0005.besa'},...
%     {'IR35','2016021711_0006.besa'},...
%     {'IR39','2016042320_0029.besa'},...
%     {'IR48','2016092915_0028.besa'},...
    };

% Filtering parameters --> FOR PLOTTING PURPOSES ONLY! Not saved this way
save_it       = 0;
resample_it   = 1;
resample_freq = 1000;
filter_it     = 1;
notch_freqs   = [60 100 120 180 200 240 300];

% basic plotting parameters
plot_it = 1;                    %if 1, assumes you're cleaning with bob and saves that way
cfg_plot = [];
cfg_plot.viewmode = 'vertical'; %trace per line
cfg_plot.continuous = 'yes';    % not trials
cfg_plot.plotlabels = 'yes';    %to know what to cut
cfg_plot.blocksize = 10;        % good for bob, can adjust

%% Processing Loop
for SBJ_ix = 1:length(SBJ_files)
    SBJ = SBJ_files{SBJ_ix}{1};
    SBJ_dir = fullfile('/home/knight/hoycw/PRJ_Stroop/data/',SBJ);
    raw_dir = fullfile(SBJ_dir, '00_raw/');
    raw_filename = fullfile(raw_dir,SBJ_files{SBJ_ix}{2});
    fprintf('============== Processing %s ==============\n',SBJ);
    
    out_dir = fullfile(SBJ_dir,'02_preproc/');
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    
    %% Load the data
%     preclean_filename = strcat(out_dir,SBJ,'_preclean.mat');
%     if exist(preclean_filename,'file')
%         load(preclean_filename);
%         data_notch = [];
%         data_notch.label = header_ecog.channel_labels;
%         data_notch.fsample = header_ecog.sample_rate;
%         data_notch.trial{1} = data_ecog;
%         data_notch.time{1} = linspace(0, size(data_ecog,2), size(data_ecog,2)); % ms time spacing
    if strcmp(raw_filename(end-4:end),'.besa')
        cfg = [];
        cfg.dataset = raw_filename;
        cfg.continuous = 'yes';
        cfg.channel = 'all';
        data_raw = ft_preprocessing(cfg);   % just load the data, don't process it
        % Define just the three fields needed for this script
        data_ecog = data_raw.trial{1};
        header_ecog.n_samples = size(data_ecog,2);
        header_ecog.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
        orig_srate = data_raw.fsample;
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
        orig_srate = header_ecog.sample_rate;
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
    
    %% Resample data to speed things up
    fprintf('============== Resampling %s ==============\n',SBJ);
    if (resample_it) && (data_raw.fsample ~= resample_freq)
        cfg_resamp = [];
        cfg_resamp.resamplefs = resample_freq;
        cfg_resamp.detrend = 'yes';
        data_resamp = ft_resampledata(cfg_resamp, data_raw);
    else
        data_resamp = data_raw;
    end
    clear data_raw
    %% Filter data for ease of viewing
    fprintf('============== Notching %s ==============\n',SBJ);
    if (filter_it)% && (data_resamp.line_noise ~= notch_freqs)
        cfg_notch = [];
        cfg_notch.continuous = 'yes';
        cfg_notch.dftfilter = 'yes'; % line noise removal using discrete fourier transform
        cfg_notch.dftfreq = notch_freqs;
        cfg_notch.demean = 'yes';
        data_notch = ft_preprocessing(cfg_notch, data_resamp);
    else
        data_notch = data_resamp;
    end
    clear data_resamp
    %% Plot data: check good channels, refs, artifacts, etc.
    if plot_it
        plot_out = ft_databrowser(cfg_plot, data_notch); %pulls up your plotted data
    else
        plot_out = [];
    end
    
    %% Set up files to prepare for KLA scripts
%     analysis_channels = 1:length(data_notch.label);
%     if strcmp(raw_filename(end-4:end),'.besa')
%         data_ecog = data_notch.trial{:}(analysis_channels,:); %only analysis_channels
%         %   BUT! analysis_channels variable should be original chan #s
%         header_ecog.sample_rate = data_notch.fsample;
%         header_ecog.orig_n_channels = size(data_notch.label,1);
%         header_ecog.n_samples = size(data_ecog,2);
%         header_ecog.length_in_seconds = size(data_notch.trial{:},2)/data_notch.fsample;
%         header_ecog.original_sample_rate = data_notch.fsample;
%         header_ecog.channel_labels = {data_notch.label{analysis_channels}};
%     elseif strcmp(raw_filename(end-3:end),'.edf')
%         data_ecog = data_ecog(analysis_channels,:);
%         header_ecog.channel_labels = header_ecog.channel_labels(analysis_channels);
%     end
%     header_ecog.n_channels = size(analysis_channels,2);
%     header_ecog.orig_channel_n = analysis_channels;
%     header_ecog.raw_file = raw_filename;
%     header_ecog.line_noise_freqs = notch_freqs;
%     if resample_it == 1
%         header_ecog.original_sample_rate = orig_srate;
%         header_ecog.sample_rate = resample_freq;
%     end
%     
    %% Save out those files
    if plot_it == 1
        out_filename = strcat(out_dir,SBJ,'_bob_clean.mat');
    else
        out_filename = strcat(out_dir,SBJ,'_preclean.mat');
    end
    if save_it == 1
        fprintf('============== Saving %s ==============\n',out_filename);
        save(out_filename, 'plot_out', 'data_notch','-v7.3');%'header_ecog', 'data_ecog','-v7.3');%, 'header_evnt', 'data_evnt');
    end
    clear data_ecog header_ecog data_notch
end

%% OLD STUFF
% % Grab microphone recording separately to export and listen
% mic_data = data_evnt(mic_channel_n,:);
% %rescale to prevent clipping, add 0.05 fudge factor
% mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);

% if strcmp(raw_filename(end-4:end),'.besa')
%     data_evnt = data_raw.trial{:}(event_channels,:);
%     header_evnt.sample_rate = data_raw.fsample;
%     header_evnt.orig_n_channels = size(data_raw.label,1);
%     header_evnt.n_samples = size(data_ecog,2);
%     header_evnt.length_in_seconds = size(data_raw.trial{:},2)/data_raw.fsample;
%     header_evnt.original_sample_rate = data_raw.fsample;
%     header_evnt.channel_labels = {data_raw.label{event_channels}};
% elseif strcmp(raw_filename(end-3:end),'.edf')
%     data_evnt = data_ecog(event_channels,:);
%     header_evnt = header_ecog;
%     header_evnt.channel_labels = header_ecog.channel_labels(event_channels);
% end
% header_evnt.n_channels = size(event_channels,2);
% header_evnt.orig_channel_n = event_channels;
% header_evnt.raw_file = raw_filename;

% data_out_filename = strcat(out_dir,raw_data_id,'.mat');
% proc_out_filename = strcat(out_dir,raw_data_id,'_proc_vars.mat');
% mic_data_filename = strcat(out_dir,raw_data_id,'_mic_recording.wav');

% save(proc_out_filename, 'analysis_time', 'analysis_channels', ...
%     'event_channels', 'reref_channels', 'reref_weights', 'photo_channel_n', 'mic_channel_n');
% audiowrite(mic_data_filename,mic_data_rescale,header_evnt.sample_rate);



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

