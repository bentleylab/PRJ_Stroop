function [] = A02_find_bad_channels(subject_id)
% Just inspect EEG data





% Set directory to pwd and make save directory
input_directory = [pwd '/../' subject_id '/' 'A01_SETFILES/'];

% Get file array
eegStruct = eeglabFileCheck('set',input_directory);
if(isempty(eegStruct))
  fprintf('\n\nERROR: No set files found in [%s]\n\n', input_directory);
  return;
end

for file_n = length(eegStruct):-1:1
  
  % Load data
  EEG = pop_loadset('filename', eegStruct(file_n).name, 'filepath', input_directory);
  fprintf('%s\tFile %2.0f of %2.0f\n', subject_id, file_n, length(eegStruct));
  
  % Plot data
  plot_seconds = 15;
  amplitude_range = 40;
  n_chans_to_display = 45;
  eegplot(EEG.data, 'srate', EEG.srate, 'winlength', plot_seconds, 'spacing', amplitude_range, 'dispchans', n_chans_to_display);
  
  pause;
  close all;
end



