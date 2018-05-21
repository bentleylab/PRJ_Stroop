%% Preprocessing Script for IR39
%   Bipolar reference, resample, notch, detrend, demean 

close all
clear all
% Add toolboxes to path
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/'));
addpath('/home/knight/hoycw/Apps/fieldtrip-20160927/'); % read .besa
ft_defaults

%% Script Variables
SBJ         = 'IR39';
raw_data_id = 'IR39_RAC';

% Logic Flags
resample_it = 1;
reref_it    = 1;
filter_it   = 1;        % also does demean and detrend
compare_it  = 0;        % keep raw data in plotting for comparison (also apply same resamp/filt)
plot_it     = 0;

% FieldTrip Processing Variables
resample_freq = 1000;
detrend       = 'no';       % minimalist view
demean        = 'no';       % see what the data looks like without these, then add if necessary
notch         = 'yes';
highpass      = 'no';
hp_freq       = 0.5;
lowpass       = 'no';
lp_freq       = 250;

SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
input_file = strcat(raw_data_id, '_data.mat');
input_filename = strcat(SBJ_dir,'01_import/',input_file);
extract_vars_filename = strcat(SBJ_dir,'01_import/',raw_data_id,'_extract_vars.mat');
bad_epoch_filename = strcat(SBJ_dir,'02_preproc/',raw_data_id,'_bad_epochs.mat');
preproc_id = strcat(raw_data_id,'_BP');     % NOTE the BP for bipolar
preproc_filename = strcat(SBJ_dir,'02_preproc/',preproc_id,'_preproc.mat');

%% Look at noise profile of all channels
% Load data
load(extract_vars_filename);
load(input_filename);

%% Loop through channels to see noise profile
% Check for really noisy channels that should be tossed
fn_plot_PSD_compare(data_ecog, header_ecog.channel_labels, header_ecog.sample_rate, 1);

%% Set line filters based on noise profiles plotted above
% Serious bump in PSD ~20 Hz!!!, nothing at 60, only blip at 300 Hz
notch_freqs = [60 120 180 240 300];

%% Resample data to speed things up
% Cut to analysis_time blocks
[data_ecog_cut, header_ecog] = fn_cut_analysis_time_blocks(data_ecog,header_ecog,analysis_time);
header_ecog.analysis_time = analysis_time;

% Move into FieldTrip
data_ft = [];
data_ft.fsample = header_ecog.sample_rate;
data_ft.label = header_ecog.channel_labels;
data_ft.time{1} = linspace(0, size(data_ecog_cut,2), size(data_ecog_cut,2));
data_ft.trial{1} = data_ecog_cut;

% Resample
fprintf('============== Resampling %s ==============\n',SBJ);
if (resample_it) && (data_ft.fsample ~= resample_freq)
    cfg_resamp = [];
    cfg_resamp.resamplefs = resample_freq;
    cfg_resamp.debug = 'saveonerror';
    data_resamp = ft_resampledata(cfg_resamp, data_ft);
else
    data_resamp = data_ft;
end
header_ecog.n_samples = size(data_resamp.trial{1},2);
header_ecog.sample_rate = data_resamp.fsample;

%% Re-reference the data
% Use a bipolar scheme
fprintf('============== Rereferencing %s ==============\n',SBJ);
if reref_it
    data_reref = [];
    data_reref.fsample = data_resamp.fsample;
    data_reref.time{1} = data_resamp.time{1};
    data_reref.trial{1} = NaN([size(data_resamp.trial{1},1)-1,size(data_resamp.trial{1},2)]);
    reref_labels = {};
    for ix = 1:size(data_reref.trial{1},1)
        % Medial-Lateral
        data_reref.trial{1}(ix,:) = data_resamp.trial{1}(ix,:)-data_resamp.trial{1}(ix+1,:);
        data_reref.label{ix} = strcat(data_resamp.label{ix},'-',data_resamp.label{ix+1});
    end
else
    data_reref = data_resamp;
end
header_ecog.orig_channel_labels = data_resamp.label;
header_ecog.channel_labels = data_reref.label;
header_ecog.n_channels = length(data_reref.label);
header_ecog.rereference_type = 'bipolar_medial-lateral';

% Concatenate original and rereferenced data to compare
if compare_it
    data_compare = [];
    data_compare.label = [data_resamp.label data_reref.label];
    data_compare.trial{1} = vertcat(data_resamp.trial{1}, data_reref.trial{1});
    final_chan = size(data_resamp.trial{1},1)+1:size(data_compare.trial{1},1);
else
    data_compare = data_reref;
    final_chan = 1:size(data_compare.trial{1},1);
end

%% Filter data
fprintf('============== Notching %s ==============\n',SBJ);
if (filter_it)% && (data_resamp.line_noise ~= notch_freqs)
    cfg_notch = [];
    cfg_notch.continuous = 'yes';
    cfg_notch.detrend    = detrend;
    cfg_notch.demean     = demean;
    cfg_notch.dftfilter  = notch; % line noise removal using discrete fourier transform
    cfg_notch.dftfreq    = notch_freqs;
    cfg_notch.lpfilter   = lowpass;
    if strcmp(lowpass,'yes')
        cfg_notch.lpfreq     = lp_freq;
        header_ecog.lp_freq  = lp_freq;
    end
    cfg_notch.hpfilter   = highpass;
    if strcmp(highpass,'yes')
        cfg_notch.hpfreq     = hp_freq;
        header_ecog.hp_freq  = hp_freq;
    end
    data_notch = ft_preprocessing(cfg_notch, data_compare);
else
    data_notch = data_compare;
end
header_ecog.detrend  = detrend;
header_ecog.demean   = demean;
header_ecog.notch    = notch;
header_ecog.lowpass  = lowpass;
header_ecog.highpass = highpass;
header_ecog.notch_freqs = notch_freqs;

%% Plot both side by side
cfg_plot = [];
cfg_plot.viewmode = 'vertical'; %trace per line
cfg_plot.continuous = 'yes';    % not trials
cfg_plot.plotlabels = 'yes';    %to know what to cut
cfg_plot.blocksize = 10;        % good for bob, can adjust
if exist(bad_epoch_filename)
    load(bad_epoch_filename);
    cfg_plot.artfctdef.visual.artifact = bad_epochs;
end

if plot_it
    ft_databrowser(cfg_plot, data_notch)
end

%% Save preprocessed data
data_preproc = data_notch;
data_preproc.trial{1} = data_preproc.trial{1}(final_chan,:);
data_preproc.label = data_preproc.label(final_chan);

% Save
save(preproc_filename, 'header_ecog', 'data_preproc', 'preproc_id');
