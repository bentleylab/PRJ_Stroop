function A04_remove_ica_components(subject_id, comps_to_remove, plot_results)






if(nargin < 3)
  plot_results = 0;
end

% Set directory to pwd and make save directory
input_directory   = [pwd '/../' subject_id '/' 'A01_SETFILES/'];
weights_directory = [pwd '/../' subject_id '/' 'A03_ICAWEIGHTS/'];
save_directory = [input_directory '../A04_CORRECTED/'];
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

% Load ICA weights for subject
weightsFileStruct = eeglabFileCheck('mat',weights_directory);
if(length(weightsFileStruct) == 1)
  load([weights_directory weightsFileStruct(1).name]);
else
  fprintf('ERROR: There should only be one .mat file in [%s]', weights_directory);
  return;
end

%% If there are no components to remove, just copy the data to the new directory
% If no bad channels, just copy data into new directory
if(length(comps_to_remove) < 1)
  for file_n = 1:n_files
    fprintf('\t%s\tNo components to remove. Copying data.\t%2.0f of %2.0f\n', subject_id, file_n, length(file_structure));
    % Load data
    EEG = pop_loadset('filename', file_structure(file_n).name, 'filepath', input_directory);
    % Save data
    EEG = pop_saveset(EEG, 'filename', file_structure(file_n).name, 'filepath', save_directory);
    % **********
    clear EEG;
  end
  return
end

%% Remove ICA components
for file_n = 1:n_files
  
  fprintf('\t%s\t%2.0f of %2.0f\n', subject_id, file_n, n_files);
  
  % Load data
  EEG = pop_loadset('filename', file_structure(file_n).name, 'filepath', input_directory);
  EEG.icasphere   = ica_data.sphere;
  EEG.icaweights  = ica_data.weights;
  EEG.icachansind = 1:EEG.nbchan;
  EEG.icawinv     = pinv(EEG.icaweights*EEG.icasphere);
  
  % Plot results
  if(plot_results)
    component_keep = setdiff(1:size(EEG.icaweights,1), comps_to_remove);
    compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
    compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);
    eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Black = channel before rejection; red = after rejection', ...
      'winlength', 15, 'data2', compproj, 'dispchans', 30);
    pause;
    close all;
  end
  
  % Remove components
  EEG = pop_subcomp( EEG, comps_to_remove, 0);
  
  % Update EEG.history
  EEG.history = [EEG.history (sprintf(['\n' 'Removed ICA Components ' num2str(comps_to_remove) '\n']))];
  
  % Save data
  pop_saveset(EEG, 'filename', file_structure(file_n).name, 'filepath', save_directory);
  % **********
  
  clear EEG;
  
  
end






