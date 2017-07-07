function cond_num_idx = convert_trial_lab2num(labels)
% Converts string labels {[con, neu, inc]; [mcon, same, minc]} to numbers
% INPUTS:
%   labels [cell array of str] - string format ${trial-block}, e.g., "con-same"
% OUTPUTS:
%   cond_num_idx [int array] - numeric condition labels, has length(n_trials)
%       con = 1-3, neu = 4-6, inc = 7-9
%       those are ordered: same, minc, mcon
%       e.g., 5 is minc-neu

cond_num_idx = [];
for lab_ix = 1:length(labels)
    switch char(labels{lab_ix})
        case 'con-same'
            cond_num_idx = [cond_num_idx 1];
        case 'con-minc'
            cond_num_idx = [cond_num_idx 2];
        case 'con-mcon'
            cond_num_idx = [cond_num_idx 3];
        case 'neu-same'
            cond_num_idx = [cond_num_idx 4];
        case 'neu-minc'
            cond_num_idx = [cond_num_idx 5];
        case 'neu-mcon'
            cond_num_idx = [cond_num_idx 6];
        case 'inc-same'
            cond_num_idx = [cond_num_idx 7];
        case 'inc-minc'
            cond_num_idx = [cond_num_idx 8];
        case 'inc-mcon'
            cond_num_idx = [cond_num_idx 9];
        otherwise
            error('Invalid condition_label');
    end
end

%         case 'con'
%             condition_idx = [condition_num==1]+[condition_num==2]+[condition_num==3];
%         case 'neu'
%             condition_idx = [condition_num==4]+[condition_num==5]+[condition_num==6];
%         case 'inc'
%             condition_idx = [condition_num==7]+[condition_num==8]+[condition_num==9];
%         case 'mcon'
%             condition_idx = [condition_num==3]+[condition_num==6]+[condition_num==9];
%         case 'same'
%             condition_idx = [condition_num==2]+[condition_num==5]+[condition_num==8];
%         case 'minc'
%             condition_idx = [condition_num==1]+[condition_num==4]+[condition_num==7];

end