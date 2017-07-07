function [data, header] = fn_KLA_rereference(data, header, reref_channels,  reref_weights)
%% Rereference with KLA's function
%   'data'            - [n_channels, n_samples] double array of the data
%   'header'          - [struct of strings, cells, lists, etc.] describing the data
%   'reref_channels'  - [cell aray of number arrays] Ex: {[1 2],[2 3],[3 4]}
%                        Each element of the cell array will be a new rereferenced channel.
%                        The numbers in each element array are the channels to use for rereferencing.
%                        The values should use the original channel numbers.
%   'reref_weights'   - [cell array of number arrays] Ex: {[-1 1],[-1 1],[-1 1]}
%                        These weights will be applied to each channel specified in reref_channels for
%                        rereferencing. Standard practice would be for the values in each cell element
%                        to add up to zero. An example for bipolar rereferencing would be [-1 1]


%% Rereference if necessary
if(~isempty(reref_channels))
  % @ is creating a function handle (map_fun) that takes x as an arg and performs the find
  map_fun = @(x) find(header.orig_channel_n==x); % function to map original channel number to extracted channel number
  data_reref = zeros(length(reref_channels),header.n_samples); % Matrix to store the rereferenced data
  header.orig_channel_labels = header.channel_labels; % Save the original channel labels
  header.channel_labels = [];
  for channel_n = 1:length(reref_channels)
    % Map the original channel numbers to the analysis channel numbers
    mapped_reref_channels = arrayfun(map_fun,reref_channels{channel_n}); %!!!CWH make sure this is indexed properly
    % Rereference the data
    for ref_channel_n = 1:length(mapped_reref_channels)
      data_reref(channel_n,:) = data_reref(channel_n,:) + ...
        reref_weights{channel_n}(ref_channel_n) * data(mapped_reref_channels(ref_channel_n),:); % Add weighted values
    end
    % Set the channel label to indicate the reference channels and weights
    header.channel_labels{channel_n} = '';
    if ~isempty(reref_name)
        header.channel_labels{channel_n} = [header.orig_channel_labels{mapped_reref_channels(1)} '-' ...
            reref_name];
    else
        if(length(mapped_reref_channels) == 2 && min(reref_weights{channel_n}) == -1 && max(reref_weights{channel_n}) == 1)
          % Special case, bipolar reference
          if(reref_weights{channel_n}(1) == 1)
            header.channel_labels{channel_n} = [header.orig_channel_labels{mapped_reref_channels(1)} '-' ...
              header.orig_channel_labels{mapped_reref_channels(2)}];
          else
              %CWH: why name it 2-1? What cases is that appropriate???
            header.channel_labels{channel_n} = [header.orig_channel_labels{mapped_reref_channels(2)} '-' ...
              header.orig_channel_labels{mapped_reref_channels(1)}];
          end
        else
          % Other reference cases
          for ref_channel_n = 1:length(mapped_reref_channels)
            header.channel_labels{channel_n} = [header.channel_labels{channel_n} ...
              '[(' num2str(reref_weights{channel_n}(ref_channel_n),'%1.2f') ')' header.orig_channel_labels{mapped_reref_channels(ref_channel_n)} ']'];
          end
        end
    end
    %fprintf('%d: %s\n', channel_n,  header.channel_labels{channel_n});
  end
  data = data_reref; clear data_reref;
  % Populate the header
  header.n_channels = size(data,1);
  header.reref_channels = reref_channels;
  header.reref_weights = reref_weights;
else
  header.orig_channel_labels = header.channel_labels;
  header.reref_channels = {};
  header.reref_weights = {};
end
