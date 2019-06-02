function [labels, colors, line_styles] = fn_condition_label_styles(factor_name)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'CNI', 'CI', 'PC'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

switch factor_name
    case 'CNI'
        labels = {'C', 'N', 'I'};
        colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'pCNI'
        labels = {'pC', 'pN', 'pI'};
        colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'CI'
        labels = {'C', 'I'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
    case 'PC'
        labels = {'MC', 'EQ', 'MI'};
        colors = {[11 255 200]./256, [255 200 111]./256, [200 111 255]./256};
        %colors = {[87 255 241]./256, [255 241 87]./256, [241 87 255]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'PC_CI'
        labels = {'C_MC', 'C_MI', 'I_MC', 'I_MI'};
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

end
