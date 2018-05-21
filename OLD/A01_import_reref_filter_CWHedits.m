function A01_import_reref_filter_CWHedits(subject_id, data_file_str, save_file_str, analysis_channels, event_channels, varargin)
% CWH Edits:
%   1) Added warnign when analysis time points are out of range of recordings
% CONSIDER:
%   - polynomial detrending???

% This function imports ECoG data from UC Irvine's NK system and saves it as a .mat file
%   - ECoG data is trimmed, resampled, highpass filtered, rereferenced, and cleaned of line noise
%   - A particular file structure is assumed:
%       Functions that this function depend on should be in a subdirectory of the pwd defined
%        by [helper_function_dir_name] - default is pwd/_TOOLBOXES/
%       Input data should be a subdirectory of a directory at the level of pwd
%        by [data_dir_name]/[input_dir_name] - default is pwd/../data/00_edf/
%       This structure canbe changed in the File Paths section of this code
%
% Outputs:
%  none - A file [subject_id].mat is saved in the output directory inside the data directory
% 
% Inputs:
%  Required
%   subject_id    - [string] Identifier for the dataset
%   data_file_str - [string] The input file name including extension (no directory)
%   save_file_str - [string] The output file name including extension (no directory)
%   analysis_channels - [int array] The ECoG channels from the original dataset to use
%   event_channels - [int array] The analog input channels with event information
%  Optional
%       These should be input as name-value pairs.
%   'resample_rate'   - [number] Resample the ECoG data to this. Default is keep original sample rate
%   'highpass_freq'   - [number] Highpass the data at this cutoff (in Hz). Default is no highpass filter
%   'analysis_time'   - [cell array of 2 element number arrays] Ex: {[20.0 400.0]}
%                        Only pull data in this time range (in seconds). Default is to keep all data
%   'reref_channels'  - [cell aray of number arrays] Ex: {[1 2],[2 3],[3 4]}
%                        Each element of the cell array will be a new rereferenced channel.
%                        The numbers in each element array are the channels to use for rereferencing.
%                        The values should use the original channel numbers.
%   'reref_weights'   - [cell array of number arrays] Ex: {[-1 1],[-1 1],[-1 1]}
%                        These weights will be applied to each channel specified in reref_channels for
%                        rereferencing. Standard practice would be for the values in each cell element
%                        to add up to zero. An example for bipolar rereferencing would be [-1 1]
%   'line_noise_freqs'- [number array] Array of frequencies where there might be line noise.
%                        This should include harmonics. It is best to look at the data first before
%                        settling on frequencies to use. Typical input would be [60 120 180 240 300 360]
% 
% Example function call:
% A01_import_reref_filter('S1','infile.edf',[5:50],[1 2],'resample_rate',1024,'highpass_freq',1.0,'analysis_time',{[15 450]});
% 
% Version: 2014-08-27
% Author:  Kris Anderson, UC Berkeley




%% Defaults
default_resample_rate    = -1;
default_highpass_freq    = -1;
default_analysis_time    = {[]};
default_reref_channels   = {};
default_reref_weights    = {};
default_line_noise_freqs = [];
default_data_dir_name    = '/data';
default_reref_name       = '';

%% Parse inputs
% A01_import_reref_filter(subject_id, data_file_str, ecog_resample_rate, highpass_cutoff, analysis_time, 
% analysis_channels, event_channels, reref_channels, reref_weights, line_noise_freqs)
ip = inputParser; ip.CaseSensitive = false;
addOptional(ip,'resample_rate',    default_resample_rate,    @isnumeric);
addOptional(ip,'highpass_freq',    default_highpass_freq,    @isnumeric);
addOptional(ip,'analysis_time',    default_analysis_time,    @iscell);
addOptional(ip,'reref_channels',   default_reref_channels,   @iscell);
addOptional(ip,'reref_weights',    default_reref_weights,    @iscell);
addOptional(ip,'line_noise_freqs', default_line_noise_freqs, @isnumeric);
addOptional(ip,'data_dir_name',    default_data_dir_name,    @ischar);
addOptional(ip,'reref_name',       default_reref_name,       @ischar);
parse(ip,varargin{:});
resample_rate    = ip.Results.resample_rate;
highpass_freq    = ip.Results.highpass_freq;
analysis_time    = ip.Results.analysis_time;
reref_channels   = ip.Results.reref_channels;
reref_weights    = ip.Results.reref_weights;
line_noise_freqs = ip.Results.line_noise_freqs;
data_dir_name    = ip.Results.data_dir_name;
reref_name       = ip.Results.reref_name;

%% File paths
helper_function_dir_name = '/_TOOLBOXES';      % Ex: pwd/_TOOLBOXES/
% data_dir_name is defined in input or is set as pwd/data by default
input_dir_name = '_EDF';                       % Ex: pwd/../data/00_edf/
output_dir_name = '01_import';                 % Ex: pwd/../data/01_import/
if(~exist(fullfile(pwd,helper_function_dir_name)))
  fprintf('\nERROR: Helper function directory does not exist.\n');
  return;
end
if(~exist(fullfile(pwd,data_dir_name,input_dir_name,data_file_str)))
  fprintf('\nERROR: Cannot find data at [%s].\n', fullfile(pwd,data_dir_name,input_dir_name,data_file_str));
  return;
end

% Import helper functions to Matlab path
addpath(genpath(fullfile(pwd,helper_function_dir_name)));

% Define path to input data and output data
data_path = fullfile(pwd,data_dir_name,input_dir_name);
%   create directory if needed
output_path = fullfile(pwd,data_dir_name,output_dir_name);
if(exist(output_path, 'dir') == 0)
  mkdir(output_path);
end

%% Import data from .edf file
% Convert analysis_time from seconds to samples
analysis_points = {};
analysis_limits = [];
orig_analysis_time = analysis_time;
if(~isempty(analysis_time{1}))
  [temp_header, ~] = extract_nk_converted_edf(fullfile(data_path,data_file_str), 1:2, [1 10]); % Quickly load data to get header info
  orig_sample_rate = temp_header.sample_rate;
  analysis_min = 99999999999;
  analysis_max = 1;
  for time_grp = 1:length(analysis_time)
      analysis_points{time_grp}(:) = analysis_time{time_grp}(:) * orig_sample_rate;
      if(min(analysis_points{time_grp}(:)) < analysis_min)
          analysis_min = min(analysis_points{time_grp}(:));
      end
      if(max(analysis_points{time_grp}(:)) > analysis_max)
          analysis_max = max(analysis_points{time_grp}(:));
      end
  end
  analysis_limits = [max(analysis_min,1) analysis_max];
  % Subtract the beginning sample time from analysis_times so they will be correct later
  for time_grp = 1:length(analysis_time)
      analysis_time{time_grp}(:) = analysis_time{time_grp}(:) - (analysis_min/orig_sample_rate);
  end
end

%     Import ECoG data and resample
[header_ecog, data_ecog] = extract_nk_converted_edf(fullfile(data_path,data_file_str), analysis_channels, analysis_limits, resample_rate);
%     Import event data and keep original sample rate
[header_evnt, data_evnt] = extract_nk_converted_edf(fullfile(data_path,data_file_str), event_channels, analysis_limits);

%% Remove extra characters from channel labels
for channel_n = 1:header_ecog.n_channels
  header_ecog.channel_labels{channel_n} = strrep(header_ecog.channel_labels{channel_n},'POL','');
  header_ecog.channel_labels{channel_n} = strrep(header_ecog.channel_labels{channel_n},'Ref','');
end

%% Highpass filter ecog data
if(highpass_freq > 0)
  fprintf('High pass filtering at %2.2fHz\n',highpass_freq);
  %   First remove the temporal mean (detrend)
  data_ecog = detrend(data_ecog','constant')'; % This assumes data is [channels x samples]
  %   Then apply highpass filter
  hpfilt_order = ceil((header_ecog.sample_rate/highpass_freq)*5);
  hpfilt_cutoff_norm = 2*highpass_freq/header_ecog.sample_rate;
  hpfilt_coeff = fir1(hpfilt_order, hpfilt_cutoff_norm, 'high');
  hpfilt_group_delay = median(grpdelay(hpfilt_coeff));
  %   Filter the data (causal) and then shift signal to make it zero-phase
  data_ecog = filter(hpfilt_coeff,1,[data_ecog zeros(header_ecog.n_channels,hpfilt_group_delay)]')';  % This assumes data is [channels x samples]
  data_ecog = data_ecog(:,(hpfilt_group_delay+1):end);
  %data_ecog = X_eegfilt(data_ecog,header_ecog.sample_rate,highpass_freq,[]); % Use EEGLAB (filtfilt which is slower)
  header_ecog.highpass = highpass_freq;
else
  %   Remove the temporal mean (detrend)
  data_ecog = detrend(data_ecog','constant')'; % This assumes data is [channels x samples]
  header_ecog.highpass = 0;
end

%% Pull out sections of data to analyze
if(~isempty(analysis_time{1}))
    all_analysis_points_ecog = [];
    all_analysis_points_evnt = [];
    for time_grp = 1:length(analysis_time)
        all_analysis_points_ecog = [all_analysis_points_ecog ...
            (analysis_time{time_grp}(1)*header_ecog.sample_rate):(analysis_time{time_grp}(2)*header_ecog.sample_rate)];
        all_analysis_points_evnt = [all_analysis_points_evnt ...
            (analysis_time{time_grp}(1)*header_evnt.sample_rate):(analysis_time{time_grp}(2)*header_evnt.sample_rate)];
    end
    all_analysis_points_ecog(all_analysis_points_ecog<1) = [];
    all_analysis_points_evnt(all_analysis_points_evnt<1) = [];
    disp('===================================================');
    disp('WARNING: some analysis points less than 0');
    disp('===================================================');
    all_analysis_points_ecog(all_analysis_points_ecog>header_ecog.n_samples) = [];
    all_analysis_points_evnt(all_analysis_points_evnt>header_evnt.n_samples) = [];
    disp('===================================================');
    disp('WARNING: some analysis points longer than recording');
    disp('===================================================');
    data_ecog = data_ecog(:,all_analysis_points_ecog);
    data_evnt = data_evnt(:,all_analysis_points_evnt);
    header_ecog.n_samples = size(data_ecog,2);
    header_evnt.n_samples = size(data_evnt,2);
    header_ecog.length_in_seconds = header_ecog.n_samples / header_ecog.sample_rate;
    header_evnt.length_in_seconds = header_evnt.n_samples / header_evnt.sample_rate;
end
header_ecog.analysis_time = orig_analysis_time;
header_evnt.analysis_time = orig_analysis_time;

%% Rereference ecog data
% Check that length of channels and weights cell arrays are the same
if(length(reref_channels)~=length(reref_weights))
  error('\nDifferent number of elements in reref_channels and reref_weights\n');
  return;
end
% Check that all channels in reref_channels are in analysis_channels
flat_reref_channels = [];
for channel_n = 1:length(reref_channels)
  flat_reref_channels = [flat_reref_channels reref_channels{channel_n}];
end
if(min(ismember(flat_reref_channels,analysis_channels)) == 0)
  error('\nCheck that reref_channels only contains channels in analysis_channels\n');
  return;
end
% Rereference if necessary
if(~isempty(reref_channels))
    % @ is creating a function handle (map_fun) that takes x as an arg and performs the find
  map_fun = @(x) find(header_ecog.orig_channel_n==x); % function to map original channel number to extracted channel number
  data_ecog_reref = zeros(length(reref_channels),header_ecog.n_samples); % Matrix to store the rereferenced data
  header_ecog.orig_channel_labels = header_ecog.channel_labels; % Save the original channel labels
  header_ecog.channel_labels = [];
  for channel_n = 1:length(reref_channels)
    % Map the original channel numbers to the analysis channel numbers
    mapped_reref_channels = arrayfun(map_fun,reref_channels{channel_n});
    % Rereference the data
    for ref_channel_n = 1:length(mapped_reref_channels)
        %CWH: Why does he multiply the weight by the ref_channel_n??? won't that increment each time and give later channels a higher weighting???
      data_ecog_reref(channel_n,:) = data_ecog_reref(channel_n,:) + ...
        reref_weights{channel_n}(ref_channel_n) * data_ecog(mapped_reref_channels(ref_channel_n),:); % Add weighted values
    end
    % Set the channel label to indicate the reference channels and weights
    header_ecog.channel_labels{channel_n} = '';
    if ~isempty(reref_name)
        header_ecog.channel_labels{channel_n} = [header_ecog.orig_channel_labels{mapped_reref_channels(1)} '-' ...
            reref_name];
    else
        if(length(mapped_reref_channels) == 2 && min(reref_weights{channel_n}) == -1 && max(reref_weights{channel_n}) == 1)
          % Special case, bipolar reference
          if(reref_weights{channel_n}(1) == 1)
            header_ecog.channel_labels{channel_n} = [header_ecog.orig_channel_labels{mapped_reref_channels(1)} '-' ...
              header_ecog.orig_channel_labels{mapped_reref_channels(2)}];
          else
              %CWH: why name it 2-1? What cases is that appropriate???
            header_ecog.channel_labels{channel_n} = [header_ecog.orig_channel_labels{mapped_reref_channels(2)} '-' ...
              header_ecog.orig_channel_labels{mapped_reref_channels(1)}];
          end
        else
          % Other reference cases
          for ref_channel_n = 1:length(mapped_reref_channels)
            header_ecog.channel_labels{channel_n} = [header_ecog.channel_labels{channel_n} ...
              '[(' num2str(reref_weights{channel_n}(ref_channel_n),'%1.2f') ')' header_ecog.orig_channel_labels{mapped_reref_channels(ref_channel_n)} ']'];
          end
        end
    end
    %fprintf('%d: %s\n', channel_n,  header_ecog.channel_labels{channel_n});
  end
  data_ecog = data_ecog_reref; clear data_ecog_reref;
  % Populate the header
  header_ecog.n_channels = size(data_ecog,1);
  header_ecog.reref_channels = reref_channels;
  header_ecog.reref_weights = reref_weights;
else
  header_ecog.orig_channel_labels = header_ecog.channel_labels;
  header_ecog.reref_channels = {};
  header_ecog.reref_weights = {};
end

%% Reduce line noise using cleanline (subtracts fitted sinusoid)
% cleanline2 is the same as cleanline from Tim Mullen but it takes
%  raw data instead of an EEGLAB dataset.
% See the comments in cleanline for an explanation of all of these
%  parameters. You might want to change them to get the best results.
if(~isempty(line_noise_freqs))
  fprintf('Removing line noise\n');
  data_ecog2 = cleanline2(data_ecog, header_ecog.sample_rate, ...
    'LineFrequencies',  line_noise_freqs, ...   % default is 60, 120
    'ScanForLines',     1, ...                  % finds exact line freq around given value
    'LineAlpha',        0.01, ...               % default = 0.01
    'Bandwidth',        3, ...                  % default = 1
    'ChanCompIndices',  1:header_ecog.n_channels, ...
    'SlidingWinLength', 4.0, ...                % in sec, default = 4
    'SlidingWinStep',   3.0, ...                % in sec, default = 4 (no overlap)
    'SmoothingFactor',  100, ...                % default=100; 1=linear, Inf=no smooth
    'PaddingFactor',    1, ...                  % default = 2; CWH good with this
    'ComputeSpectralPower', 0, ...              % default = 1/True; might be nice
    'PlotFigures', 1, ...                       % maybe alternative to CompSpecPow
    'VerboseOutput',    0);
  header_ecog.line_noise_freqs = line_noise_freqs;
else
  header_ecog.line_noise_freqs = [];
end
                   
%% Save output data
header_ecog.subject_id = subject_id; header_evnt.subject_id = subject_id;
fprintf('Saving [%s]\n',fullfile(output_path,save_file_str));
save(fullfile(output_path,save_file_str), 'header_ecog', 'data_ecog', 'header_evnt', 'data_evnt', 'subject_id');

















