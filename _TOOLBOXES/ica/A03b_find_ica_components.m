function A03b_find_ica_components(subject_id, comps_to_plot)
% 
% This function is for visual inspections of ICA component topographies and activations.
% 
% Use it to determine which components to remove for artifact correction.
% 
% Remember to be conservative and not take out brain signals. You will be interpolating channels
%  again after this, so if a component is fixed on a frontal channel but looks like it has brain
%  activity, don't choose it for deletion.
% 




n_comps_to_plot = 20;

if(nargin < 2)
  comps_to_plot = 1:n_comps_to_plot;
end

% Set directory to pwd and make save directory
input_directory   = [pwd '/../' subject_id '/' 'A01_SETFILES/'];
weights_directory = [pwd '/../' subject_id '/' 'A03_ICAWEIGHTS/'];

% Get file array
file_structure = eeglabFileCheck('set',input_directory);
n_files = length(file_structure);
if(isempty(file_structure))
  fprintf('\n\nERROR: No set files found in [%s]\n\n', input_directory);
  return;
end

% Load ICA weights for subject
weights_file_structure = eeglabFileCheck('mat',weights_directory);
if(length(weights_file_structure) == 1)
  load([weights_directory weights_file_structure(1).name]);
else
  fprintf('ERROR: There should only be one .mat file in [%s]', weights_directory);
  return;
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
fprintf('Combining blocks');
for file_n = 1:n_files
  tempEEG = pop_loadset('filename', file_structure(file_n).name, 'filepath', input_directory);
  EEG.data(:,(last_point+1):(last_point+samples_per_block(file_n))) = tempEEG.data(:,(ica_params.edgecutoff*EEG.srate):(end-ica_params.edgecutoff*EEG.srate));
  last_point = last_point+samples_per_block(file_n);
  clear tempEEG;
  fprintf('.');
end
fprintf('\n');

% Add ICA information to EEG structure
EEG.icasphere   = ica_data.sphere;
EEG.icaweights  = ica_data.weights;
EEG.icachansind = 1:EEG.nbchan;
EEG.icawinv     = pinv(EEG.icaweights*EEG.icasphere);

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

% Plot component information. This is the same for all blocks
% pop_selectcomps(EEG, 1:n_comps_to_plot ); % This function was incorrectly displaying the components
for comp_n = comps_to_plot
  pop_prop( EEG, 0, comp_n, [], {'freqrange', [1 50]} );
  pause;
  close all;
end

% Plot topoplots of each component
%figure('Position', [25 700 1400 250]);
pop_topoplot(EEG, 0, comps_to_plot , [subject_id ' Block ' num2str(file_n) ' of ' num2str(n_files)] , ...
  [2 10] , 0, 'electrodes', 'off', 'colorbar', 'off');%, 'figureposition', [25 800 1400 250]);
  
% Plot data
plot_seconds = 10;
n_chans_to_display = 20;
tmpEEGdata = eeg_getdatact(EEG, 'component', (1:(size(EEG.icaweights,1)-numel(sort(unique([bad_channels ica_params.ref_channel]))))));
amplitude_range = ( max(tmpEEGdata(1,(20*EEG.srate):(end-20*EEG.srate)))-min(tmpEEGdata(1,(20*EEG.srate):(end-20*EEG.srate))) ) / 2;
eegplot(tmpEEGdata, 'srate', EEG.srate, 'winlength', plot_seconds, 'spacing', amplitude_range, 'dispchans', n_chans_to_display);
clear tempEEGdata;
clear EEG;

pause;
close gcf;

