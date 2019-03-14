function [header, data] = extract_nk_converted_edf_CWH(file_path, channels_to_save, samples_to_save, resample_rate)
% Colin W. Hoy adapted this function for personal use (May 22, 2016)
% Changes:
%   1) Personal file paths (e.g., add EDFData.m directory to path)
%   2) Saves the .mat file (without running a separate script)
%
% This function converts .edf files created with the .besa to .edf conversion manager (V01-04_20131112)
% from Nihon Kohden into a Matlab data structure
% 
% INPUTS:
%  file_path           - [string], The path to the file to be extracted
%  channels_to_save    - [integer array], Channel numbers that you want to be extracted
%                                          input [] for all channels.
%  samples_to_save     - [2 element integer array], Sample range from original data to extract.
%                                          input [] for all samples.
%  [resample_rate]     - Optional, [integer], Data will be resampled using Matlab's resample() function
% 
% OUTPUTS:
%  data                          - [channels x samples], The ECoG data
%  header                        - [structure], Information about the ECoG data
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
% 
% Caution: Extract only the channels that you need. Certain versions of Matlab cannot save files larger than 2GB
%   (assuming you are going to save the output). Check your Matlab preferences [General->MAT-Files] to see what
%   version of .mat files will be saved.
% 
% Caution 2: The data is scaled correctly, but the DC offsets might be off between channels.
% 
% Caution 3: This function uses Matlab's resample() function from the signal processing toolbox.
% 
% External function necessary for this script:
%  EDFData.m, created by J.B.Wagenaar, Copyright 2013 Trustees of the University of Pennsylvania
%     - This script uses a modified version of EDFData.m. Specifically the order of the matrix output
%        from getUnscaledData and getData are changed from [pnts x chans] to [chans x pnts]
% 
% Kris Anderson, UC Berkeley

if(length(samples_to_save) ~= 2 && ~isempty(samples_to_save))
  error('samples_to_save should be two element vector [start end] with start > end or []\n');
end
if(~isempty(samples_to_save) && samples_to_save(2) < samples_to_save(1))
  error('samples_to_save should be two element vector [start end] with start > end or []\n');
end
if(~isempty(samples_to_save) && min(samples_to_save) < 1)
  error('samples_to_save should be two element vector [start end] with start > end or []\n');
end
if(~isempty(channels_to_save) && min(channels_to_save) < 1)
  error('channels_to_save should array of channel numbers to save or []\n');
end

%% SET UP DIRECTORIES
disp(' ');
SBJ = input('Enter subject ID:  ','s');
SBJ_dir = ['/home/knight/hoycw/Data_Extraction/Data/' SBJ '/'];
raw_dir = [SBJ_dir 'raw_data/'];
if(exist(raw_dir, 'dir') == 0)
    mkdir(raw_dir);
end
clean_dir = [SBJ_dir 'clean/'];
if(exist(clean_dir, 'dir') == 0)
    mkdir(clean_dir);
end

%% EXTRACT EDF AND CONVERT TO MAT
% Call the function to read the .edf file
addpath('/home/knight/hoycw/Data_Extraction/Scripts/extractEDF/');
fprintf('Reading [%s]\n',file_path);
edf_object = EDFData(file_path);

% Pull out the useful header information
header.recording_startdate = num2str(edf_object.startDate);
header.recording_starttime = num2str(edf_object.startTime);
header.orig_n_channels = edf_object.ns;
if(isempty(channels_to_save))
  % If input is [] for channels_to_save, save all channels
  channels_to_save = 1:header.orig_n_channels;
end
header.n_channels = numel(channels_to_save);
if(header.n_channels > header.orig_n_channels)
  error(['Requesting more channels than exist in the data. File contains ' num2str(header.orig_n_channels) ' channels.']);
end
header.n_samples = edf_object.records * edf_object.samples(1);
if(isempty(samples_to_save))
  samples_to_save = [1 header.n_samples];
end
header.orig_samples_extracted = samples_to_save;
if(~isempty(samples_to_save))
  if(samples_to_save(2) > header.n_samples || samples_to_save(1) < 1)
    error(['samples_to_save is out of range. Max sample number is ' num2str(header.n_samples)]);
  end
end
if(range(edf_object.samples) ~= 0)
  error('Header information shows different number of samples in each channel.\nThis function does not work with that type of file.');
end
header.length_in_seconds = edf_object.records * edf_object.duration;
header.original_sample_rate = edf_object.sf(1);
if(range(edf_object.sf) ~= 0)
  error('Header information shows different sample rates in each channel.\nThis function does not work with that type of file.');
end
header.pid = edf_object.patientID;
header.rid = edf_object.recordID;
for channel_n = 1:header.n_channels
  header.channel_labels{channel_n} = edf_object.label{channels_to_save(channel_n)};
  header.orig_channel_n(channel_n) = channels_to_save(channel_n);
end

% Check if resampling is required
if(nargin < 4)
  resample_rate = header.original_sample_rate;
end
if(resample_rate <= 0)
  error('resample_rate must be a positive number');
end
header.sample_rate = resample_rate;

% Extract data
fprintf('Collecting requested channels\n');
% Extract specified samples
data = edf_object.getData([samples_to_save(1) samples_to_save(2)],channels_to_save);
header.n_samples = size(data,2);
header.length_in_seconds = header.n_samples / header.original_sample_rate;

% Resample
if(header.sample_rate ~= header.original_sample_rate)
  fprintf('Resampling from %dHz to %dHz\n',header.original_sample_rate,header.sample_rate);
  resample_data = zeros(header.n_channels, ceil(header.n_samples*header.sample_rate/header.original_sample_rate));
  fprintf('.');
  for channel_n = 1:header.n_channels
    resample_data(channel_n,:) = resample(data(channel_n,:),header.sample_rate,header.original_sample_rate,100); % 100 increases accuracy for resample()
    fprintf('.');
  end
  fprintf('\n');
  data = resample_data;
  clear resample_data;
  header.n_samples = size(data,2);
end

%% SAVE MAT

