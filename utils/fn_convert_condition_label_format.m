function alt_label = fn_convert_condition_label_format(label)
% Converts between CWH and KLA style condition labels
% Inputs:
%   label [str] - label in either KLA or CWH format
%       CWH format = ${trial_type}_${block_type}
%       KLA format = ${block_type}-${trial_type}
% Outputs:
%   alt_label [str] - label in the alternate format from label

if isempty(strfind(label,'_'))
    % KLA style, switch to CWH
    dash = strfind(label,'-');
    alt_label = [label(dash+1:end) '_' label(1:dash-1)];
else
    % CWH style, switch to KLA
    uscore = strfind(label,'_');
    alt_label = [label(uscore+1:end) '-' label(1:uscore-1)];
end
end
