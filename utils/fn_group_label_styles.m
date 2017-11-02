function [labels, colors, line_styles] = fn_group_label_styles(model_id)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'CNI', 'CI', 'pcon'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

% if length(cond_lab) == 1
switch model_id
    case 'RT_CNI_pcon'
        labels = {'RT', 'CNI', 'pcon'};
        colors = {[127,201,127]./256, [190,174,212]./256, [253,192,134]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'rRT_CNI_pcon'
        labels = {'CNI', 'pcon'};
        colors = {[190,174,212]./256, [253,192,134]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
    case 'corrRT_CNI_pcon'
        labels = {'CNI', 'pcon'};
        colors = {[190,174,212]./256, [253,192,134]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
    case 'RT'
        labels = {'RT'};
        colors = {[127,201,127]./256};
        line_styles = {'-'};
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
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',labels{:}));
end
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

end
