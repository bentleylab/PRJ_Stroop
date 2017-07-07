function epochs_norm = fn_ft_bsln_epochs(epochs, bsln_times, bsln_type)
% Takes one time series over the whole experiment and cuts into epochs
% INPUTS:
%   epochs [2D array] - n_epochs by time array of cut epochs
%   bsln_times [int, int]- [start, end] of baseline period in SECONDS
%       will be matched against epochs.time field
%   bsln_type [str]   - type of baseline to implement
%       'zscore' = subtract mean and divide by SD
%       'demean' = subtract mean
% OUTPUTS:
%   epochs_norm [2D array] - same epochs but baseline corrected

epochs_norm = epochs;

% Cutoff signals by time range
for k = 1:numel(epochs.trial)
    % Find baseline time points
    bsln_ix = [find(epochs.time{k}==bsln_times(1)) find(epochs.time{k}==bsln_times(2))];
    
    switch bsln_type
        case 'zscore'
            bsln_mean = nanmean(epochs.trial{k}(:,bsln_ix(1):bsln_ix(2)),2);
            bsln_std  = nanstd(epochs.trial{k}(:,bsln_ix(1):bsln_ix(2)),[],2);
            epochs_norm.trial{k} = (epochs.trial{k}-repmat(bsln_mean,1,size(epochs.trial{1},2)))/bsln_std;
        case 'demean'
            bsln_mean = nanmean(epochs.trial{k}(:,bsln_ix(1):bsln_ix(2)),2);
            epochs_norm.trial{k} = epochs.trial{k}-repmat(bsln_mean,1,size(epochs.trial{1},2));
        case 'none'
            epochs_norm.trial{k} = epochs.trial{k};
        otherwise
            error(strcat('unknown bsln_type: ',bsln_type));
    end
end

epochs_norm.cfg.previous.bsln_type = bsln_type;
epochs_norm.cfg.previous.baselinewindow = bsln_times;

end




