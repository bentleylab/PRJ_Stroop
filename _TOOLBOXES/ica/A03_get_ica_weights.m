function A03_get_ica_weights(subject_id, bad_blocks, bad_channels)
% This function computes the ICA weight matrix for all of the .set files for a subject.
%  First it combines all of the blocks into one big EEG file and then it calls runica().
%  The output is only the weight matrix, etc and not the EEG data or ICA activations.
%  The output file can be used in a later function to remove ICA artifacts or whatever.
% 
% This function will take a while to run, and runica is pretty verbose.


ica_params.max_channel_n = 70;
ica_params.momentum = 0.05;
ica_params.maxsteps = 768;
ica_params.stoplearning = 0.000000005;
ica_params.edgecutoff = 10;           % In seconds. This will throw out the edge data from each block

valid_events_bin = {'10110000', ... % VIS L CUE
                    '10010000', ... % VIS R CUE
                    '01110000', ... % TAC L CUE
                    '01010000', ... % TAC R CUE
                    '10101000', ... % VIS L SAMPLE
                    '10001000', ... % VIS R SAMPLE
                    '01101000', ... % TAC L SAMPLE
                    '01001000', ... % TAC R SAMPLE
                    '10100100', ... % VIS L PROBE
                    '10000100', ... % VIS R PROBE
                    '01100100', ... % TAC L PROBE
                    '01000100', ... % TAC R PROBE
                    '10100010', ... % VIS L BLINK
                    '10000010', ... % VIS R BLINK
                    '01100010', ... % TAC L BLINK
                    '01000010' ...  % TAC R BLINK
                    };
valid_event_distance = 1.5; % In seconds
data_max_value = 250; % Data is ignored if it's greater than this amount
max_bad_channels = 0; % Ignore data if more channels than this are bad at a certain time point

max_eeg_channel_n = 64; % used for determining bad data

bad_channels = sort(bad_channels);
bad_blocks = sort(bad_blocks);

% Determine valid events
for event_n = 1:length(valid_events_bin)
  valid_events(event_n) = bin2dec(valid_events_bin{event_n}); %#ok<AGROW>
end

% Set directory to pwd and make save directory
input_directory = [pwd '/../' subject_id '/' 'A01_SETFILES/'];
save_directory = [input_directory '../A03_ICAWEIGHTS/'];
if(exist(save_directory, 'dir') == 0)
  mkdir(save_directory);
end

% Get file array
file_structure = eeglabFileCheck('set',input_directory);
n_files = length(file_structure);
if(isempty(file_structure))
  fprintf('\n\nERROR: No set files found in [%s]\n\n', input_directory);
  return;
end

% Remove bad blocks from file array
if(numel(bad_blocks) > 0)
  % Search through file_structure array for the bad block
  found_bad_blocks = 0;
  for file_n = n_files:-1:1
    for bad_block_n = 1:numel(bad_blocks)
      if(~isempty(strfind(file_structure(file_n).name, ['_B' num2str(bad_blocks(bad_block_n))])))
        file_structure(file_n) = [];
        found_bad_blocks = found_bad_blocks+1;
        break;
      end
    end
  end
  if(found_bad_blocks ~= numel(bad_blocks))
    fprintf('\n\nERROR: Bad block not found. Looking for %2.0f, found %2.0f\n\n', found_bad_blocks, numel(bad_blocks));
    return;
  end
  fprintf('\tRemoved %2.0f bad blocks\n', found_bad_blocks);
  n_files = length(file_structure);
end

% Determine size of all EEG blocks
samples_per_block = zeros(1,n_files);
for file_n = 1:n_files
  EEG = pop_loadset('filename', file_structure(file_n).name, 'filepath', input_directory);
  samples_per_block(file_n) = size(EEG.data(:,(ica_params.edgecutoff*EEG.srate):(end-ica_params.edgecutoff*EEG.srate)),2);
  clear EEG;
end

% Load a base data set and preallocate big array
EEG = pop_loadset('filename', file_structure(1).name, 'filepath', input_directory);
EEG.data = zeros(EEG.nbchan, sum(samples_per_block));
EEG.pnts = sum(samples_per_block);

% Combine all blocks into one big dataset
last_point = 0;
for file_n = 1:n_files
  tempEEG = pop_loadset('filename', file_structure(file_n).name, 'filepath', input_directory);
  EEG.data(:,(last_point+1):(last_point+samples_per_block(file_n))) = tempEEG.data(:,(ica_params.edgecutoff*EEG.srate):(end-ica_params.edgecutoff*EEG.srate));
  last_point = last_point+samples_per_block(file_n);
  clear tempEEG;
end

% Save reference channel number
ica_params.ref_channel = EEG.ref;

% Remove bad channels from data
fprintf('Removing bad channels for ICA\t[ ');
for channel_n = 1:numel(bad_channels)
  fprintf('%d ', bad_channels(channel_n));
end
fprintf(']\n');
ica_params.all_channels = 1:ica_params.max_channel_n;
ica_params.ica_channels = ica_params.all_channels; % Set ica_channels to all channels
ica_params.ica_channels(sort(unique([bad_channels ica_params.ref_channel]))) = [];           % Remove bad_channels from ica_channels
eeg_channels = ica_params.ica_channels(ica_params.ica_channels<=max_eeg_channel_n); % Set eeg_channels to all ica_channels less than the largest eeg channel
eeg_channels = 1:numel(eeg_channels);                                               % Since EEG channels are always first (1:64 usually), and bad channels have been pulled, make eeg_channels just the length of ica_channels less than maximum eeg channel number
EEG.data = EEG.data(ica_params.ica_channels,:);
EEG.nbchan = size(EEG.data,1);
EEG.chanlocs = EEG.chanlocs(ica_params.ica_channels);

fprintf('Original data length: %7.0f points\n', EEG.pnts); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Find data that is within a certain time from an event
good_times = zeros(EEG.pnts,1);
for event_n = 1:length(EEG.event)
  if(max(EEG.event(event_n).type == valid_events))
    curr_lat = EEG.event(event_n).latency;
    good_times(ceil(curr_lat-(valid_event_distance*EEG.srate)):ceil(curr_lat+(valid_event_distance*EEG.srate))) = 1;
  end
end

fprintf('After non-event data removed: %7.0f\n', sum(good_times)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Find really bad data (greater than a certain number of microvolts)
for point_n = 1:EEG.pnts
  n_bad_channels = 0;
  for channel_n = eeg_channels % Only EEG channels
    if(abs(EEG.data(channel_n,point_n)) > data_max_value)
      n_bad_channels = n_bad_channels+1;
    end
  end
  if(n_bad_channels > max_bad_channels)
    good_times(point_n) = 0;
  end
end
 
% Pull out data that is near a valid event and doesn't exceed threshold
new_data = zeros(EEG.nbchan,sum(good_times));
last_point = 0;
for point_n = 1:EEG.pnts
  if(good_times(point_n) > 0)
    last_point = last_point+1;
    new_data(:,last_point) = EEG.data(:,point_n);
  end
end
EEG.data = new_data;
EEG.pnts = size(EEG.data,2);
clear new_data;
 
fprintf('After extreme artifact data removed: %7.0f\n', EEG.pnts); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine data rank (for ICA)
ica_params.data_rank = rank(double(EEG.data));
if(ica_params.data_rank < size(EEG.data,1))
  fprintf('\n\nERROR: Data rank is %d, which is less than the number of (non-bad) channels %d. PCA needs to reduce data size first\n\n', ica_params.data_rank, size(EEG.data,1));
  return;
end

% Calculate ICA weights
function_start_time = clock;
[ica_data.weights, ica_data.sphere, ica_data.compvars, ica_data.bias, ica_data.signs] = ...
  runica(EEG.data, ...
  'extended', 1, 'stop', ica_params.stoplearning, 'maxsteps', ica_params.maxsteps, ...
  'momentum', ica_params.momentum);
function_stop_time = clock;
fprintf('%s\tComputed ICA weights\t%2.0f of %2.0f\tElapsed time: %4.0f seconds\n', subject_id, file_n, length(file_structure), etime(function_stop_time, function_start_time));


% Fix weights matrix, etc to add back in bad_channels
removed_channels = sort(unique([bad_channels ica_params.ref_channel]));
used_channels = 1:ica_params.max_channel_n;
used_channels(removed_channels) = [];
for channel_n = 1:numel(removed_channels)
  current_channel = removed_channels(channel_n);
  ica_data.signs(end+1) = 1;
  ica_data.bias(end+1) = 0;
  ica_data.compvars(end+1) = min(ica_data.compvars)/1.1;
  % Weights matrix
  ica_data.weights = [ica_data.weights(:,1:(removed_channels(channel_n)-1)) zeros(size(ica_data.weights,1),1) ica_data.weights(:,removed_channels(channel_n):end)];
  ica_data.weights(end+1,:) = 0;
  ica_data.weights(end,current_channel) = 1;
  % Sphere matrix (put 1 in diagonal)
  ica_data.sphere = [ica_data.sphere(1:(removed_channels(channel_n)-1),:) ; zeros(1,size(ica_data.sphere,2)) ; ica_data.sphere(removed_channels(channel_n):end,:)];
  ica_data.sphere = [ica_data.sphere(:,1:(removed_channels(channel_n)-1)) zeros(size(ica_data.sphere,1),1) ica_data.sphere(:,removed_channels(channel_n):end)];
  ica_data.sphere(current_channel,current_channel) = 1;
end


block_string_indices = strfind(file_structure(file_n).name, '_B');
save_file_name = [file_structure(file_n).name(1:(block_string_indices(1)-1)) '.mat'];
save([save_directory save_file_name], 'ica_data', 'ica_params', 'bad_blocks', 'bad_channels', 'subject_id', 'good_times');

