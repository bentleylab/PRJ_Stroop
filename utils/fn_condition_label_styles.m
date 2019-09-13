function [labels, colors, line_styles] = fn_condition_label_styles(factor_name)
%% Converts the name of a set of conditions into labels, plotting colors/styles
% condition_name: [str] 'CNI', 'CI', 'PC'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3

switch factor_name
    case 'CNI'
        labels = {'C', 'N', 'I'};
        colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'pCI'
        labels = {'pC', 'pI'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
    case 'pCNI'
        labels = {'pC', 'pN', 'pI'};
        colors = {[55,126,184]./256, [0 0 0], [228,26,28]./256};
        line_styles = {'-', '-', '-'};    % colors for cond_lab plotting
    case 'SDR'
        labels = {'S', 'D', 'R'};   % light blue, pink/red, gold (vivid)
        colors = {[3 209 255]./256, [255 3 108]./256, [255 205 3]./256};
        line_styles = {'-', '-', '-'};
    case 'ConCSPC'
        labels = {'Con', 'CS', 'PC'};   % gold, teal, pink (light but sharp)
        colors = {[255 200 111]./256, [11 255 200]./256, [200 111 255]./256};
        line_styles = {'-', '-', '-'};
    case 'CSPC'
        labels = {'CS', 'PC'};
        colors = {[11 255 200]./256, [200 111 255]./256};
        line_styles = {'-', '-'};
    case 'CI'
        labels = {'C', 'I'};
        colors = {[55,126,184]./256, [228,26,28]./256};
        line_styles = {'-', '-'};    % colors for cond_lab plotting
    case 'PC'
        labels = {'MC', 'EQ', 'MI'};
        colors = {[27 158 119]./256, [217 95 2]./256, [117 112 179]./256};
        %colors = {[87 255 241]./256, [255 241 87]./256, [241 87 255]./256};    % too light
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
