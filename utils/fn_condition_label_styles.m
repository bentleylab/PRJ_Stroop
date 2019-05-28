function [labels, colors, line_styles] = fn_condition_label_styles(factor_name)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'CNI', 'CI', 'pcon'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

% if length(cond_lab) == 1
switch factor_name
    case 'CNI'
        labels = {'con', 'neu', 'inc'};
        colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'CI'
        labels = {'con', 'inc'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
    case 'pcon'
        labels = {'mcon', 'same', 'minc'};
        colors = {[11 255 200]./256, [255 200 111]./256, [200 111 255]./256};
        %colors = {[87 255 241]./256, [255 241 87]./256, [241 87 255]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'pcon_CI'
        labels = {'con_mcon', 'con_minc', 'inc_mcon', 'inc_minc'};
        colors = {[0 0 0], [0 0 0], [228,26,28]./256, [228,26,28]./256};    % colors for cond_lab plotting
        line_styles = {'-', '--', '-','--'};    % colors for cond_lab plotting
    case 'CSE'
        labels = {'cI','iI'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'actv'
        labels = {'actv'};
        colors = {'k'};
        line_styles = {'-'};
    otherwise
        error(strcat('Only one, unrecognized condition offered: ',factor_name));
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
