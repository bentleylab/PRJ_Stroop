function baseDesign = fn_create_PERM_design(trial_info)
    
    % Initialize cell array for previous trial types
    prev_type = cell(1, numel(trial_info.block_n));
    
    % Determine the type of the previous trial for each trial
    for itrial = 1:numel(trial_info.block_n)
        if itrial == 1 || trial_info.block_n(itrial) ~= trial_info.block_n(itrial-1) || trial_info.trial_n(itrial) - trial_info.trial_n(itrial-1) ~= 1
            prev_type{itrial} = 'None';
        else
            prev_type{itrial} = trial_info.trialtype{itrial - 1};
        end
    end
    prev_type = prev_type'; % Convert to column vector to match the orientation of trial_info.trialtype
         
    % Exclude neutral trials and construct the base design matrix
    baseDesign = table();
    baseDesign.CurrentType = trial_info.trialtype;
    baseDesign.BlockType = trial_info.blocktype;
    baseDesign.PreviousType = prev_type;
    baseDesign.RT = trial_info.response_time;
    baseDesign.resp_onset = trial_info.resp_onset;
    baseDesign.error = trial_info.error;
    
    

end