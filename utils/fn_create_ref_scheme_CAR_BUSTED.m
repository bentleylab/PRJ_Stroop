THIS DOESN"T WORK!!!
function [weights] = fn_create_ref_scheme_CAR(labels,CAR_exclude)
%% Create a re-reference scheme using common average reference from a list of labels
% INPUTS:
%   labels [cell array of str] - names of the channels
%   CAR_exclude [cell array of str] - names of channels in labels that
%       should not be included in the average that is subtracted from all channels
%
% OUTPUTS:
%   weights = [length(labels),length(labels)] sized matrix of weights to combine elecs

% Find valid channels to average
valid_CAR_ch = setdiff(labels,CAR_exclude);

% Create montage
CAR_weights = zeros(size(labels));
for lab_ix = 1:numel(labels)
    % If it's not excluded, give this a negative weight proportional to # good channels
    if isempty(strmatch(labels{lab_ix},CAR_exclude))
        CAR_weights(lab_ix) = -1/numel(valid_CAR_ch);
    end
end

% Combine weights for each channel into a matrix with 1 on the diagonal
weights = repmat(CAR_weights,[1 numel(labels)]);
for lab_ix = 1:numel(labels)
    weights(lab_ix,lab_ix) = 1;
end

end
