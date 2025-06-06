function [labels, colors, line_styles] = fn_group_label_styles(model_id)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'CNI', 'CI', 'PC'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible labels and their colors
conditions  = {'RT','pCNI','CNI','PC','PCi'};
cond_colors = {[77 175 74]./256, ...
               [228 26 28]./256, [228 26 28]./256, ...
               [55 126 184]./256, [55 126 184]./256};
% Taken from: colorbrewer2.org, qualitative, 5-class Set1
%   currently: green, red, blue
%   purple for later: [152 78  163]
%   orange for later: [255 127 0]

%% Convert model_id into set of conditions
switch model_id
%     case 'RT_CNI_PC'
%         labels = {'RT', 'CNI', 'PC'};
    case 'rRT_CNI_PC'
        labels = {'CNI', 'PC'};
    case 'crRT_CNI_PC'
        labels = {'CNI', 'PC'};
    case 'pCNI_PC'
        labels = {'pCNI', 'PC'};
    case 'CNI_PC'
        labels = {'CNI', 'PC'};
    case 'CNI_p_PC'
        labels = {'CNI', 'pCNI', 'PC'};
    case 'PC'
        labels = {'PC'};
    case 'PCi'
        labels = {'PC'};
    case 'CNI'
        labels = {'CNI'};
    case 'RT'
        labels = {'RT'};
    case 'actv'
        labels = {};
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',labels{:}));
end

% Assign colors and line styles
colors = cell(size(labels));
line_styles = cell(size(labels));
for cond_ix = 1:numel(labels)
    colors{cond_ix} = cond_colors{strcmp(conditions,labels{cond_ix})};
    line_styles{cond_ix} = '-';
end

end

