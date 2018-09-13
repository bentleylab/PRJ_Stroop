function [new_lab] = fn_thin_labels(labels,thinner)
% Removes repeated labels.
% INPUTS:
%   labels [cell array strs] - labels to be thinned
%   thinner [int] - amount to thin the labels
%       0: not at all
%       1: immediately repeating labels
%       2: empty until a new label

last_lab = 'placeHolder';
new_lab = labels;
for l = 1:numel(labels)
    cur_lab = labels{l};
    if thinner && strcmp(cur_lab,last_lab)
        new_lab{l} = ' ';
        if thinner==2
            last_lab = ' ';
        end
    else
        last_lab = cur_lab;
    end
end

end
