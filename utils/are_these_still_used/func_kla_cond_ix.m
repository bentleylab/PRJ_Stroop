function cond_ix = func_kla_cond_ix(trial_info, LABELS)
% Returns struct of binary indices of trial condition labels from kla hdr
% INPUTS:   trial_info [struct] = kla assembled trial header
%                   .condition_types = numbering of all conditions
%                   .condition_n = [n_trials] vector of ints labeling by condition_n
%           LABELS [cell array, str] = condition supersets to return indices
%                   {con, inc, neu, mcon, same, minc}

% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
% trial_type = NaN(size(trial_info.condition_n));
for lab_ix = 1:length(LABELS)
    switch LABELS{lab_ix}
        % Trial Type
        case 'con'
            cond_ix.con = logical([trial_info.condition_n==1]+[trial_info.condition_n==2]+[trial_info.condition_n==3]);
        case 'neu'
            cond_ix.neu = logical([trial_info.condition_n==4]+[trial_info.condition_n==5]+[trial_info.condition_n==6]);
        case 'inc'
            cond_ix.inc = logical([trial_info.condition_n==7]+[trial_info.condition_n==8]+[trial_info.condition_n==9]);

        % Proportion Congruency
        case 'same'
            cond_ix.same = logical([trial_info.condition_n==1]+[trial_info.condition_n==4]+[trial_info.condition_n==7]);
        case 'minc'
            cond_ix.minc = logical([trial_info.condition_n==2]+[trial_info.condition_n==5]+[trial_info.condition_n==8]);
        case 'mcon'
            cond_ix.mcon = logical([trial_info.condition_n==3]+[trial_info.condition_n==6]+[trial_info.condition_n==9]);

        % Individual trial block combinations
        case 'same_con'
            cond_ix.same_con = logical(trial_info.condition_n==1);
        case 'minc_con'
            cond_ix.minc_con = logical(trial_info.condition_n==2);
        case 'mcon_con'
            cond_ix.mcon_con = logical(trial_info.condition_n==3);
        case 'same_neu'
            cond_ix.same_neu = logical(trial_info.condition_n==4);
        case 'minc_neu'
            cond_ix.minc_neu = logical(trial_info.condition_n==5);
        case 'mcon_neu'
            cond_ix.mcon_neu = logical(trial_info.condition_n==6);
        case 'same_inc'
            cond_ix.same_inc = logical(trial_info.condition_n==7);
        case 'minc_inc'
            cond_ix.minc_inc = logical(trial_info.condition_n==8);
        case 'mcon_inc'
            cond_ix.mcon_inc = logical(trial_info.condition_n==9);

        % Error
        otherwise
            disp(strcat('Invalid condition label: ',LABELS{lab_ix}));
            error(strcat('Invalid condition label: ',LABELS{lab_ix}));
    end
end

end



% for lab_ix = 1:length(LABELS)
%     if length(LABELS{lab_ix})==3   % Catch con, neu, inc
%         cond_ns = strfind(trial_info.condition_types, strcat('-',LABELS{lab_ix}));
%     else                            % Catch mcon, same, minc
%         cond_ns = strfind(trial_info.condition_types, LABELS{lab_ix});
%     end
%     cond_ns = 
%     % Loop over all conditions that aooly to super label
%     for cond_n_ix = 1:length(cond_ns)
%         cond_ix.LABELS{lab_ix} = strcmp(trial_info.condition_n, cond_ns(cond_n_ix));
%     end
% end

