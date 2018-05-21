function epochs_norm = fn_bsln_epochs(epochs, events, bsln_lim, bsln_type)
% Takes one time series over the whole experiment and cuts into epochs
% INPUTS:
%   epochs [2D array] - n_epochs by time array of cut epochs
%   buff_lim [int, int]- 2 int array of data indices for [start, end] of baseline period
%   bsln_type [str]   - type of baseline to implement
%       'zscore' = subtract mean and divide by SD
%       'demean' = subtract mean
% OUTPUTS:
%   bsln_epochs [2D array] - same epochs but baseline corrected

% Cutoff signals by time range
epochs_norm = NaN(size(epochs));
for k = 1:size(epochs,1)
    trial = epochs(k,:);
    switch bsln_type
        case 'zscore'
            bsln_mean = nanmean(trial(events(k)-bsln_lim(1):events(k)+bsln_lim(2)));
            bsln_std  = nanstd(trial(events(k)-bsln_lim(1):events(k)+bsln_lim(2)));
            epochs_norm(k,:) = (trial-bsln_mean)/bsln_std;
        case 'demean'
            bsln_mean = nanmean(trial(events(k)-bsln_lim(1):events(k)+bsln_lim(2)));
            epochs_norm(k,:) = trial-bsln_mean;
        case 'none'
            epochs_norm(k,:) = trial;
        otherwise
            error(strcat('unknown bsln_type: ',bsln_type));
    end
end

end




