function condition_num_str = fn_convert_condition_lab2num(label)
% Converts condition string labels to number index
% Inputs:
%   label [str] - label in either KLA or CWH format
%       CWH format = ${trial_type}_${block_type}
%       KLA format = ${block_type}-${trial_type}
% Outputs:
%   condition_num [char] - char list of numbers for trial_info condition numbering
%       e.g., mcon-con = '3'; con = '123'

% Convert to CWH style label if necessary
if (isempty(strfind(label,'_'))) && (length(label) > 3)
    label = fn_convert_condition_label_format(label);
end
switch label
    case 'con'
        condition_num_str = [num2str(1) num2str(2) num2str(3)];
    case 'neu'
        condition_num_str = [num2str(4) num2str(5) num2str(6)];
    case 'inc'
        condition_num_str = [num2str(7) num2str(8) num2str(9)];
    case 'mcon'
        condition_num_str = [num2str(3) num2str(6) num2str(9)];
    case 'same'
        condition_num_str = [num2str(1) num2str(4) num2str(7)];
    case 'minc'
        condition_num_str = [num2str(2) num2str(5) num2str(8)];
    case 'con_same'
        condition_num_str = [num2str(1)];
    case 'con_minc'
        condition_num_str = [num2str(2)];
    case 'con_mcon'
        condition_num_str = [num2str(3)];
    case 'neu_same'
        condition_num_str = [num2str(4)];
    case 'neu_minc'
        condition_num_str = [num2str(5)];
    case 'neu_mcon'
        condition_num_str = [num2str(6)];
    case 'inc_same'
        condition_num_str = [num2str(7)];
    case 'inc_minc'
        condition_num_str = [num2str(8)];
    case 'inc_mcon'
        condition_num_str = [num2str(9)];
    otherwise
        error('Invalid condition_label');
end



end