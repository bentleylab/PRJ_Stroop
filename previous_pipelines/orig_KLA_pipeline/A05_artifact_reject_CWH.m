function A05_artifact_reject_CWH(data_id, time_range, artifact_struct, plot_chans)
% time_range is [start,end] times in sec, where start and end define an
%   epoch around the event of interest (e.g., -10:10 around word onset)
% plot_chans is index into the data_ecog and header_ecog arrays
% 

if(nargin < 4)
  plot_chans = [];
end

data_file_str = [data_id '.mat'];
u_score_pos = strfind(data_id,'_');
SBJ = data_id(1:u_score_pos(1)-1);

%% File paths
helper_function_dir_name = '/home/knight/hoycw/PRJ_Stroop/scripts/_TOOLBOXES/';
data_dir_name = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
input_dir_name = '04_proc';
output_dir_name = '06_epochs';
if(~exist(fullfile(helper_function_dir_name)))
  fprintf('\nERROR: Helper function directory does not exist.\n');
  return;
end
if(~exist(fullfile(data_dir_name,input_dir_name,data_file_str)))
  fprintf('\nERROR: Cannot find data at [%s].\n', fullfile(data_dir_name,input_dir_name,data_file_str));
  return
end

% Import helper functions to Matlab path
addpath(genpath(fullfile(helper_function_dir_name)));
% Define path to input data and output data
data_path = fullfile(data_dir_name,input_dir_name);
%   create directory if needed
output_path = fullfile(data_dir_name,output_dir_name);
if(exist(output_path, 'dir') == 0)
  mkdir(output_path);
end

%% Load input data
fprintf('Loading %s\n',fullfile(data_path,data_file_str));
load(fullfile(data_path,data_file_str));
n_channels = header_ecog.n_channels;
s_rate = header_ecog.sample_rate;

plot_struct.plot_chans = plot_chans;
plot_struct.plot_chan_names = header_ecog.channel_labels(plot_chans);

% Convert time_range from times to points
epoch_range = round(time_range(1)*s_rate):round(time_range(2)*s_rate);

%% Epoch data
n_epochs = length(trial_info.trial_n);
n_samples = length(epoch_range);
event_onsets = trial_info.word_onset;
data_epoched = zeros(n_channels,n_epochs,n_samples);
for epoch_n = 1:n_epochs
  data_epoched(:,epoch_n,:) = data_ecog(:,event_onsets(epoch_n)+epoch_range);
end
clear data_ecog;

%% Artifact rejection
% Reject artifacts from raw data
ok_epochs1 = artifact_reject(data_epoched, artifact_struct.std_limit_raw, artifact_struct.hard_threshold_raw, s_rate, plot_struct);
% Reject artifacts from diff data (looking for fast changes)
ok_epochs2 = artifact_reject(diff(data_epoched,1,3), artifact_struct.std_limit_diff, artifact_struct.hard_threshold_diff, s_rate, plot_struct);
% Get intersection of the two arrays
ok_epochs = intersect(ok_epochs1,ok_epochs2);
fprintf('\tTOTAL %d OUT OF %d EPOCHS RETAINED\n', length(ok_epochs), n_epochs);

%% Save output
artifact_struct.time_range = time_range;
time_range_str = [num2str(floor(time_range(1)*1000),'%+d') '.' num2str(floor(time_range(2)*1000),'%+d')];
save(fullfile(output_path,[data_id '_' time_range_str '.mat']), 'header_ecog', ...
  'ok_epochs', 'event_onsets', 'trial_info', 'artifact_struct', 'data_id');


%% Artifact rejection function
function [ok_epochs] = artifact_reject(current_data, std_limit, hard_threshold, s_rate, plot_struct)

n_channels = size(current_data,1);
n_epochs = size(current_data,2);

% Remove epochs that exceed a certain threshold
ok_epochs = 1:n_epochs;
threshold_epochs = squeeze(max(max(abs(current_data),[],3),[],1)) > hard_threshold;
ok_epochs(threshold_epochs) = [];

% Remove temporal mean
for channel_n = 1:n_channels
  for epoch_n = 1:n_epochs
    current_data(channel_n,epoch_n,:) = detrend(squeeze(current_data(channel_n,epoch_n,:)),'constant');
  end
end
% Keep removing epochs until all are below standard deviation threshold
%   Initialize variables for loop
n_bad_epochs_last = -1;
n_bad_epochs = 0;
exceeded_channels_last1 = -ones(n_channels,1);
exceeded_channels_last2 = zeros(n_channels,1);
%   Start checking EEG Standard Deviations
[eeg_std] = get_std_of_data(current_data(:,ok_epochs,:), [], s_rate);
std_last = eeg_std(:,1,:); clear eeg_std;
iter_n = 0;
while( n_bad_epochs > n_bad_epochs_last )
  iter_n = iter_n+1;
  n_bad_epochs_last = n_bad_epochs;
  % Get standard deviation of data (excluding bad epochs)
  [eeg_std] = get_std_of_data(current_data(:,ok_epochs,:), [], s_rate); clear ok_epochs; clear bad_epochs;
  std_data = eeg_std(:,1,:); clear eeg_std;
  % Replace standard deviation on channels whose number of bad epochs has not increased
  for channel_n = 1:n_channels
    if(exceeded_channels_last1(channel_n) <= exceeded_channels_last2(channel_n))
      % Number of bad epochs for channel has not increased, use old std
      std_data(channel_n,1,:) = std_last(channel_n,1,:);
    end
  end
  % Find bad epochs
  [ok_epochs, bad_epochs, exceeded_channels] = find_exceeded_std(current_data, std_data, std_limit, s_rate);
  n_bad_epochs = numel(bad_epochs);
  std_last = std_data;
  exceeded_channels_last2 = exceeded_channels_last1;
  exceeded_channels_last1 = exceeded_channels;
  
  fprintf('\tITERATION %2.0f\tOK: %3.0f\tREJECTED: %3.0f\tPCT REJ: %3.0f\n', iter_n, numel(ok_epochs), n_bad_epochs, (n_bad_epochs/n_epochs)*100);
end
fprintf('\tOK: %3.0f\tREJECTED: %3.0f\tPCT REJ: %3.0f\n', numel(ok_epochs), numel(bad_epochs), (n_bad_epochs/n_epochs)*100);

% Plot results
if(~isempty(plot_struct.plot_chans))
  plot_results(current_data, std_data, std_limit, s_rate, ok_epochs, plot_struct, exceeded_channels);
end

%% Find epochs with data that exceeds std limit function
function [ok_epochs, bad_epochs, exceeded_channels] = find_exceeded_std(eeg_data, std_data, std_limit, s_rate)
% Find data that exceeds threshold

% Get the absolute value of the difference between each sample point and the ensemble mean
%   Ensemble mean will be the robust mean (using a tukey window on the sorted values for each time point)
%   This is less affected by outliers
n_channels = size(eeg_data,1);
n_epochs = size(eeg_data,2);
n_samples = size(eeg_data,3);
ensemble_mean = zeros(n_channels, 1, n_samples);
weight_window = tukeywin(n_epochs,0.40)'; % 60% is square windowed and 20% on each side is tapered
weight_window = weight_window ./ sum(weight_window);
for channel_n = 1:n_channels
  for sample_n = 1:n_samples
    sorted_data = sort(eeg_data(channel_n,:,sample_n));
    ensemble_mean(channel_n,1,sample_n) = weighted_mean(sorted_data,weight_window);
  end
end
% Smooth the ensemble mean by first doing a moving maximum and then a moving average
movmax_window = ones(floor(s_rate*0.050),1); % Window for moving max is about 50ms
smooth_length = (length(movmax_window)*2)/n_samples; % This is a proportion compared to data length. Twice as long as movmax_window.
pad_length = length(movmax_window); % The number of samples to add to each side before smoothing to mitigate edge effects
for channel_n = 1:n_channels
  % Pad the data on each side;
  curr_data = padarray(squeeze(ensemble_mean(channel_n,1,:)),pad_length,'symmetric','both');
  curr_data = smooth(curr_data,smooth_length,'moving');
  curr_data = imdilate(curr_data,movmax_window);
  curr_data = smooth(curr_data,smooth_length,'moving');
  ensemble_mean(channel_n,1,:) = curr_data((pad_length+1):(end-pad_length));
end
abs_eeg_data = abs(eeg_data-repmat(ensemble_mean,[1 size(eeg_data,2) 1]));

exceeded_eeg_data = zeros(size(abs_eeg_data));

n_channels = size(eeg_data,1);
for channel_n = 1:n_channels
  for epoch_n = 1:n_epochs
    exceeded_eeg_data(channel_n,epoch_n,:) = abs_eeg_data(channel_n,epoch_n,:) > (std_data(channel_n,:,:)*std_limit);
  end
end
all_epochs = max(max(exceeded_eeg_data,[],3),[],1);
n_bad_epochs = 0;
n_ok_epochs = 0;
ok_epochs = [];
bad_epochs = [];
for epoch_n = 1:n_epochs
  if(all_epochs(epoch_n)>0)
    n_bad_epochs=n_bad_epochs+1;
    bad_epochs(n_bad_epochs)=epoch_n;
  else
    n_ok_epochs=n_ok_epochs+1;
    ok_epochs(n_ok_epochs)=epoch_n;
  end
end

% Determine channels that exceeded threshold
exceeded_channels = zeros(n_channels,1);
for channel_n = 1:n_channels
  exceeded_channels(channel_n) = sum(max(exceeded_eeg_data(channel_n,:,:),[],3));
end

%% Get standard deviation function
function [eeg_std] = get_std_of_data(eeg_data, bad_channels, s_rate)

n_channels = size(eeg_data,1);
n_samples = size(eeg_data,3);

% Calculate the ensemble standard deviation
eeg_std = std(eeg_data,[],2);

% Smooth the standard deviation by first doing a moving maximum and then a moving average
movmax_window = ones(floor(s_rate*0.050),1); % Window for moving max is about 50ms
smooth_length = (length(movmax_window)*2)+1; % This is a proportion compared to data length. Twice as long as movmax_window.
pad_length = length(movmax_window); % The number of samples to add to each side before smoothing to mitigate edge effects
for channel_n = 1:n_channels
  % Pad the data on each side;
  curr_data = padarray(squeeze(eeg_std(channel_n,1,:)),pad_length,'symmetric','both');
  curr_data = smooth(curr_data,smooth_length,'moving');
  curr_data = imdilate(curr_data,movmax_window);
  curr_data = smooth(curr_data,smooth_length,'moving');
  eeg_std(channel_n,1,:) = curr_data((pad_length+1):(end-pad_length));
end

% Make standard deviation in 'bad_channels' equal a large number
if(~isempty(bad_channels))
  eeg_std(bad_channels,1,:) = 1000;
end

%% Plot results function
function plot_results(eeg_data, std_data, std_limit, s_rate, ok_epochs, plot_struct, exceeded_channels)

plot_chans = plot_struct.plot_chans;

figure;%('Position', [20 50 1200 900]);

% Scale the data by dividing by the largest standard deviation in each channel
max_amp = squeeze(max(std_data*std_limit,[],3));
min_amp = squeeze(min(std_data*-std_limit,[],3));
scale_amp = max(max_amp,min_amp);
eeg_data =  eeg_data ./ repmat(scale_amp,[1 size(eeg_data,2) size(eeg_data,3)]);
std_data =  std_data ./ repmat(scale_amp,[1 size(std_data,2) size(std_data,3)]);

%   Ensemble mean will be the robust mean (using a tukey window on the sorted values for each time point)
%   This is less affected by outliers
n_channels = size(eeg_data,1);
n_epochs = size(eeg_data,2);
n_samples = size(eeg_data,3);
ensemble_mean = zeros(n_channels, 1, n_samples);
weight_window = tukeywin(n_epochs,0.40)'; % 60% is square windowed and 20% on each side is tapered
weight_window = weight_window ./ sum(weight_window);
for channel_n = 1:n_channels
  for sample_n = 1:n_samples
    sorted_data = sort(eeg_data(channel_n,:,sample_n));
    ensemble_mean(channel_n,1,sample_n) = weighted_mean(sorted_data,weight_window);
  end
end
% Smooth the ensemble mean by first doing a moving maximum and then a moving average
movmax_window = ones(floor(s_rate*0.050),1); % Window for moving max is about 50ms
smooth_length = (length(movmax_window)*2)/n_samples; % This is a proportion compared to data length. Twice as long as movmax_window.
pad_length = length(movmax_window); % The number of samples to add to each side before smoothing to mitigate edge effects
for channel_n = 1:n_channels
  % Pad the data on each side;
  curr_data = padarray(squeeze(ensemble_mean(channel_n,1,:)),pad_length,'symmetric','both');
  curr_data = smooth(curr_data,smooth_length,'moving');
  curr_data = imdilate(curr_data,movmax_window);
  curr_data = smooth(curr_data,smooth_length,'moving');
  ensemble_mean(channel_n,1,:) = curr_data((pad_length+1):(end-pad_length));
end
mean_eeg_data = ensemble_mean;

% Before
subplot(1,2,1); hold on;
set(gca, 'Position', [0.07 0.16 0.45 0.83]);
last_max = 0;
for chan_n = 1:numel(plot_chans)
  y_origin(chan_n) = last_max - min(min(min(std_data(plot_chans(chan_n),:,:)*-std_limit))) + 0.25;
  curr_std_data = std_data(plot_chans(chan_n),:,:)*std_limit;
  plot(squeeze(eeg_data(plot_chans(chan_n),:,:))'+y_origin(chan_n));
  plot(squeeze(mean_eeg_data(plot_chans(chan_n),:,:)+curr_std_data)+y_origin(chan_n), 'r','LineWidth', 2);
  plot(squeeze(mean_eeg_data(plot_chans(chan_n),:,:)-curr_std_data)+y_origin(chan_n),'r','LineWidth', 2);
  last_max = y_origin(chan_n) + max(max(max(std_data(plot_chans(chan_n),:,:)*std_limit)));
end
ylim([-0.2 last_max+0.2]);
xlim([0 size(eeg_data,3)]);
set(gca,'YTick', y_origin, 'YTickLabel', plot_struct.plot_chan_names);
set(gca,'XTick', []);
set(gca,'FontSize',5);

% After
subplot(1,2,2); hold on;
set(gca, 'Position', [0.53 0.16 0.45 0.83]);
for chan_n = 1:numel(plot_chans)
  plot(squeeze(eeg_data(plot_chans(chan_n),ok_epochs,:))'+y_origin(chan_n));
end
ylim([-0.2 last_max+0.2]);
xlim([0 size(eeg_data,3)]);
set(gca,'YTick', []);
set(gca,'XTick', []);

% Rejected epochs per channel
if length(exceeded_channels) > 1
  subplot(50,1,50); hold on;
  set(gca, 'Position', [0.02 0.03 0.96 0.12]);
  bar(exceeded_channels);
  xlim([0 length(exceeded_channels)]);
  set(gca,'XTick', [1 5:5:(length(exceeded_channels)-1) length(exceeded_channels)]);
end

hold off;
pause;
close all;


function y = weighted_mean(x,w,dim)
%WMEAN   Weighted Average or mean value.
%   For vectors, WMEAN(X,W) is the weighted mean value of the elements in X
%   using non-negative weights W. For matrices, WMEAN(X,W) is a row vector 
%   containing the weighted mean value of each column.  For N-D arrays, 
%   WMEAN(X,W) is the weighted mean value of the elements along the first 
%   non-singleton dimension of X.
%
%   Each element of X requires a corresponding weight, and hence the size 
%   of W must match that of X.
%
%   WMEAN(X,W,DIM) takes the weighted mean along the dimension DIM of X. 
%
%   Class support for inputs X and W:
%      float: double, single
%
%   Example:
%       x = rand(5,2);
%       w = rand(5,2);
%       wmean(x,w)
%
% Copyright (c) 2009, Adam Auton
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


if nargin<2
    error('Not enough input arguments.');
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if nargin==2, 
  % Determine which dimension SUM will use
  dim = find(size(x)~=1, 1 );
  if isempty(dim), dim = 1; end
end

y = sum(w.*x,dim)./sum(w,dim);















































































