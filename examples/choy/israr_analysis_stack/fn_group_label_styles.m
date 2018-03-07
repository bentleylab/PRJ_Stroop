function [labels, colors, line_styles] = fn_group_label_styles(model_id)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'CNI', 'CI', 'pcon'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

%% List of possible labels and their colors
conditions  = {'RT','CNI','pcon'};
cond_colors = {[77 175 74]./256, [228 26 28]./256, [55 126 184]./256};
% Taken from: colorbrewer2.org, qualitative, 5-class Set1
%   currently: green, red, blue
%   purple for later: [152 78  163]
%   orange for later: [255 127 0]

%% Convert model_id into set of conditions
switch model_id
    case 'RT_CNI_pcon'
        labels = {'RT', 'CNI', 'pcon'};
    case 'rRT_CNI_pcon'
        labels = {'CNI', 'pcon'};
    case 'corrRT_CNI_pcon'
        labels = {'CNI', 'pcon'};
    case 'RT'
        labels = {'RT'};
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',labels{:}));
end

% Assign colors and line styles
colors = {};
line_styles = {};
for cond_ix = 1:numel(labels)
    colors{cond_ix} = cond_colors{strmatch(labels{cond_ix},conditions,'exact')};
    line_styles{cond_ix} = '-';
end

end

%% Extra conditions
%     case 'CI'
%         labels = {'con', 'inc'};
%         colors = {[55,126,184]./256, [228,26,28]./256};
%         line_styles = {'-', '-'};    % colors for cond_lab plotting
%     case 'pcon'
%         labels = {'mcon', 'same', 'minc'};
%         colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
%         line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
%     case 'pcon_CI'
%         labels = {'con_mcon', 'con_minc', 'inc_mcon', 'inc_minc'};
%         colors = {[0 0 0], [0 0 0], [228,26,28]./256, [228,26,28]./256};    % colors for cond_lab plotting
%         line_styles = {'-', '--', '-','--'};    % colors for cond_lab plotting
%     case 'conseq'
%         cond_id = 'conseq';

% else
%     cond_id = 'cst';
%     for c_ix = 1:length(cond_lab)
%         cond_id = [cond_id fn_convert_condition_lab2num(cond_lab{c_ix})];
%     end
%     cond_colors = manual_cond_colors;
%     if length(cond_lab)~=length(cond_colors)
%         error('Mismatched condition labels and colors');
%     end
% end
