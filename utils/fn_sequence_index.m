function sequence_idx = fn_sequence_index(condition_lab1, condition_lab2, condition_num, block_num)
% Returns binary truth array of trials matching a certain condition
% INPUTS:
%   condition_lab1 [str] - 1st trial type in sequence
%       one of [C, N, I, MC, EQ, MI, ${trial-block}]
%       ${trial-block} can be any combination, e.g., "C-EQ"
%   condition_lab2 [str] - 2nd trial type in sequence; same as condition_lab2
%   condition_num [int array] - numeric condition labels, has length(n_trials)
%       C = 1-3, N = 4-6, I = 7-9
%       those are ordered: EQ, MI, MC
%       e.g., 5 is MI-N
%   block_num [int array] - numeric block number, length(n_trials)
%       used to exclude the first trial in each block
% OUTPUTS:
%   sequence_idx [binary array] - 0/1 for identity with condition_label, length(n_trials)

cond_idx1 = fn_condition_index(condition_lab1, condition_num);
cond_idx2 = fn_condition_index(condition_lab2, condition_num);
sequence_idx = zeros(size(condition_num));
block_start_ix = find(diff(block_num)>0)+1;
for ix = 2:length(condition_num)
    if isempty(find(block_start_ix==ix,1))
        if (cond_idx2(ix)==1) && (cond_idx1(ix-1)==1)
            sequence_idx(ix) = 1;
        end
    end
end

end
