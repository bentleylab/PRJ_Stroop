function epochs_norm = fn_epoch_bsln_combo(signal, starts, ends, buff_lim, epoch_type, bsln_events, bsln_lim, bsln_type)
% Takes one time series over the whole experiment and cuts into epochs
% Also baselines that data based on any other arbitrary epoch
% INPUTS:
%   signal [1D array] - one time series, to be cut into epochs
%   starts [1D array] - list of data indices at which epochs start
%   ends [1D array]   - list of data indices at which epochs end
%   buff_lim [int, int]- 2 int array of # data points to include [before, after] the event
%   epoch_type [str] - type of epoch cuts ('datapad' or 'nanpad')
%   bsln_lim [int, int]- 2 int array of data indices for [start, end] of baseline period
%   bsln_type [str]   - type of baseline to implement
%       'zscore' = subtract mean and divide by SD
%       'demean' = subtract mean
%       'none'   = no baseline correction
% OUTPUTS:
%   epochs [2D array] - an n_epochs by time array of cut epochs

% Cut Trials and Baseline Epochs
if strcmp(epoch_type,'nanpad')
    epochs = fn_epoch_cuts_nanpad(signal,starts,ends,buff_lim);
    bslns  = fn_epoch_cuts_nanpad(signal,bsln_events,bsln_events,bsln_lim);
elseif strcmp(epoch_type,'datapad')
    epochs = fn_epoch_cuts_datapad(signal,starts,ends,buff_lim);
    bslns  = fn_epoch_cuts_datapad(signal,bsln_events,bsln_events,bsln_lim);
else
    error('Bad epoch_type');
end


% Apply Baseline
epochs_norm = NaN(size(epochs));
for k = 1:size(epochs,1)
    trial = epochs(k,:);
    bsln = bslns(k,:);
    switch bsln_type
        case 'zscore'
            epochs_norm(k,:) = (trial-nanmean(bsln))/nanstd(bsln);
        case 'demean'
            epochs_norm(k,:) = trial-nanmean(bsln);
        case 'none'
            epochs_norm(k,:) = trial;
        otherwise
            error(strcat('unknown bsln_type: ',bsln_type));
    end
end

end




