function condition_idx = fn_condition_index(condition_label, condition_num)
% Returns binary truth array of trials matching a certain condition
% INPUTS:
%   condition_label [str] - one of [con, neu, inc, mcon, same, minc, ${trial-block}]
%       ${trial-block} can be any combination, e.g., "con-same"
%   condition_num [int array] - numeric condition labels, has length(n_trials)
%       con = 1-3, neu = 4-6, inc = 7-9
%       those are ordered: same, minc, mcon
%       e.g., 5 is minc-neu
% OUTPUTS:
%   condition_idx [binary array] - 0/1 for identity with condition_label, length(n_trials)


switch condition_label
    case 'con'
        condition_idx = [condition_num==1]+[condition_num==2]+[condition_num==3];
    case 'neu'
        condition_idx = [condition_num==4]+[condition_num==5]+[condition_num==6];
    case 'inc'
        condition_idx = [condition_num==7]+[condition_num==8]+[condition_num==9];
    case 'mcon'
        condition_idx = [condition_num==3]+[condition_num==6]+[condition_num==9];
    case 'same'
        condition_idx = [condition_num==1]+[condition_num==4]+[condition_num==7];
    case 'minc'
        condition_idx = [condition_num==2]+[condition_num==5]+[condition_num==8];
    case 'con_same'
        condition_idx = [condition_num==1];
    case 'con_minc'
        condition_idx = [condition_num==2];
    case 'con_mcon'
        condition_idx = [condition_num==3];
    case 'neu_same'
        condition_idx = [condition_num==4];
    case 'neu_minc'
        condition_idx = [condition_num==5];
    case 'neu_mcon'
        condition_idx = [condition_num==6];
    case 'inc_same'
        condition_idx = [condition_num==7];
    case 'inc_minc'
        condition_idx = [condition_num==8];
    case 'inc_mcon'
        condition_idx = [condition_num==9];
    otherwise
        error('Invalid condition_label');
end

end