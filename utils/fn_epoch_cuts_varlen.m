function epochs = fn_epoch_cuts_varlen(signal, starts, ends, buff_lim)
% Takes one time series over the whole experiment and cuts into epochs
% INPUTS:
%   signal [1D array] - one time series, to be cut into epochs
%   starts [1D array] - list of data indices at which epochs start
%   ends [1D array]   - list of data indices at which epochs end
%   buff_lim [int, int]- 2 int array of # data points to include [before, after] the event
% OUTPUTS:
%   epochs [cell array] - an n_epochs cell array, each cell has a n_channels
%                           by time array cut to that specific epoch length

% Check starts and ends are same length
if length(starts) ~= length(ends)
    error('starts and ends variables are different lengths');
end
% Cutoff signals by time range
% epochs = NaN(length(starts),length(signal(starts(1)-buff_lim(1):ends(1)+buff_lim(2))));
for k = 1:length(starts)
    epochs{k} = signal(starts(k)-buff_lim(1):ends(k)+buff_lim(2));
end

end




