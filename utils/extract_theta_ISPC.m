function ISPC_cell_array = extract_theta_ISPC(zISPC, frex)
    % Function to extract ISPC values for each MPFC-LPFC channel pair,
    % averaged over 4-8 Hz frequencies.
    %
    % Inputs:
    %   zISPC - 5-D array of ISPC values with dimensions:
    %           trials x MPFC_channels x LPFC_channels x frequencies x time
    %   frex  - Frequency vector corresponding to the frequencies dimension in zISPC
    %
    % Output:
    %   ISPC_cell_array - Cell array where each cell contains a [trials x time]
    %                     matrix for each MPFC-LPFC channel pair
    
    % Find indices of frequencies between 4 and 8 Hz
    freq_idx = frex >= 4 & frex <= 8;
    
    % Average over the selected frequencies
    zISPC_avg_freq = mean(zISPC(:, :, :, freq_idx, :), 4); 
    
    % Get the number of MPFC and LPFC channels
    [~, nMPFC_channels, nLPFC_channels, ~, ~] = size(zISPC_avg_freq);
    
    % Initialize the cell array
    ISPC_cell_array = cell(nMPFC_channels * nLPFC_channels, 1);
    cell_index = 1;
    
    % Loop over each MPFC-LPFC channel pair
    for m = 1:nMPFC_channels
        for l = 1:nLPFC_channels
            % Extract the trials x time matrix for the channel pair
            % Squeeze to remove singleton dimensions
            ISPC_cell_array{cell_index} = squeeze(zISPC_avg_freq(:, m, l, :, :));
            cell_index = cell_index + 1;
        end
    end
end
