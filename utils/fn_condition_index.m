function condition_idx = fn_condition_index(condition_label, condition_num)
% Returns binary truth array of trials matching a certain condition
% INPUTS:
%   condition_label [str] - one of [C, N, I, MC, EQ, MI, ${trial-block}]
%       ${trial-block} can be any combination, e.g., "C-EQ"
%   condition_num [int array] - numeric condition labels, has length(n_trials)
%       C = 1-3, N = 4-6, I = 7-9
%       those are ordered: EQ, MI, MC
%       e.g., 5 is MI-N
% OUTPUTS:
%   condition_idx [binary array] - 0/1 for identity with condition_label, length(n_trials)


switch condition_label
    case 'cI'
        C_idx = [condition_num==1]+[condition_num==2]+[condition_num==3];
        I_idx = [condition_num==7]+[condition_num==8]+[condition_num==9];
        condition_idx = zeros(size(condition_num));
        for t = 2:numel(condition_num)
            if C_idx(t-1)==1 && I_idx(t)==1
                condition_idx(t) = 1;
            end
        end
    case 'iI'
        C_idx = [condition_num==1]+[condition_num==2]+[condition_num==3];
        I_idx = [condition_num==7]+[condition_num==8]+[condition_num==9];
        condition_idx = zeros(size(condition_num));
        for t = 2:numel(condition_num)
            if I_idx(t-1)==1 && I_idx(t)==1
                condition_idx(t) = 1;
            end
        end
    case 'C'
        condition_idx = [condition_num==1]+[condition_num==2]+[condition_num==3];
    case 'N'
        condition_idx = [condition_num==4]+[condition_num==5]+[condition_num==6];
    case 'I'
        condition_idx = [condition_num==7]+[condition_num==8]+[condition_num==9];
    case 'MC'
        condition_idx = [condition_num==3]+[condition_num==6]+[condition_num==9];
    case 'EQ'
        condition_idx = [condition_num==1]+[condition_num==4]+[condition_num==7];
    case 'MI'
        condition_idx = [condition_num==2]+[condition_num==5]+[condition_num==8];
    case 'C_EQ'
        condition_idx = [condition_num==1];
    case 'C_MI'
        condition_idx = [condition_num==2];
    case 'C_MC'
        condition_idx = [condition_num==3];
    case 'N_EQ'
        condition_idx = [condition_num==4];
    case 'N_MI'
        condition_idx = [condition_num==5];
    case 'N_MC'
        condition_idx = [condition_num==6];
    case 'I_EQ'
        condition_idx = [condition_num==7];
    case 'I_MI'
        condition_idx = [condition_num==8];
    case 'I_MC'
        condition_idx = [condition_num==9];
    otherwise
        error('Invalid condition_label');
end

end
