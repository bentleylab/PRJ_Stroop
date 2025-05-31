function [designOutput1, designOutput2] = fn_create_LMM_design(trial_info, numElectrodesLPFC, numElectrodesMPFC, mode)
        
    % Initialize the PC column with zeros
    PC = zeros(size(trial_info.blocktype));
    
    % Loop through each blocktype and assign the corresponding PC value
    for i = 1:length(trial_info.blocktype)
        switch trial_info.blocktype{i}
            case 'minc'
                PC(i) = 0.167;
            case 'same'
                PC(i) = 0.333;
            case 'mcon'
                PC(i) = 0.5;
        end
    end

    % New trial_n per block
    block_end_idxs = find([diff(trial_info.block_n); 0]);
    block_end_idxs = [block_end_idxs; length(trial_info.block_n)];

    block_start_idxs = find([1; diff(trial_info.block_n)]);

    trial_info.norm_trial_n = NaN(length(trial_info.trial_n),1);
    for blocki = 1:length(block_start_idxs)
        
        start_n = double(trial_info.trial_n(block_start_idxs(blocki)));
        end_n = double(trial_info.trial_n(block_end_idxs(blocki)));

        current_block_trial_n = double(trial_info.trial_n(block_start_idxs(blocki):block_end_idxs(blocki)));
        
        norm_trial_n = (current_block_trial_n - start_n)./(end_n - start_n);

        trial_info.norm_trial_n(block_start_idxs(blocki):block_end_idxs(blocki)) = norm_trial_n;
    end
    % Add the PC column to the struct
    trial_info.PC = PC;
    
    % Initialize cell array for previous trial types
    prev_type = cell(numel(trial_info.block_n),1);
    prev_rt = NaN(numel(trial_info.block_n),1);

    % Determine the type of the previous trial for each trial
    for itrial = 1:numel(trial_info.block_n)
        if itrial == 1 || trial_info.block_n(itrial) ~= trial_info.block_n(itrial-1) || (trial_info.trial_n(itrial) - trial_info.trial_n(itrial-1)) ~= 1
            prev_type{itrial} = 'None';
        else
            prev_type{itrial} = trial_info.trialtype{itrial - 1};
            prev_rt(itrial) = trial_info.response_time(itrial - 1);
        end
    end
    
    % Initialize cell array for post-trial types
    post_type = cell(numel(trial_info.block_n),1);
    post_rt = NaN(numel(trial_info.block_n),1);

    % Determine the type of the post-trial for each trial
    for itrial = 1:numel(trial_info.block_n)
        if itrial == numel(trial_info.block_n) || trial_info.block_n(itrial) ~= trial_info.block_n(itrial+1) || (trial_info.trial_n(itrial+1) - trial_info.trial_n(itrial)) ~= 1
            post_type{itrial} = 'None';
        else
            post_type{itrial} = trial_info.trialtype{itrial + 1};
            post_rt(itrial) = trial_info.response_time(itrial + 1);
        end
    end
   
    % Make the base design matrix
    baseDesign = table();
    baseDesign.CurrentType = trial_info.trialtype;
    baseDesign.PC = trial_info.PC;
    baseDesign.PreviousType = prev_type;
    baseDesign.PostType = post_type;
    baseDesign.RT = trial_info.response_time;
    baseDesign.PostRT = post_rt;
    baseDesign.PrevRT = prev_rt;
    baseDesign.BlockType = trial_info.blocktype;
    baseDesign.nTrial = trial_info.trial_n;
    baseDesign.norm_nTrial = trial_info.norm_trial_n;
    
    
    if strcmp(mode, 'power')
        % The original functionality for power analysis
        lpfcDesign = [];
        for electrode = 1:numElectrodesLPFC
            tempDesign = baseDesign;
            tempDesign.Electrode = repmat(electrode, height(baseDesign), 1);
            lpfcDesign = [lpfcDesign; tempDesign];
        end
        
        mpfcDesign = [];
        for electrode = 1:numElectrodesMPFC
            tempDesign = baseDesign;
            tempDesign.Electrode = repmat(electrode, height(baseDesign), 1);
            mpfcDesign = [mpfcDesign; tempDesign];
        end
        
        designOutput1 = lpfcDesign;
        designOutput2 = mpfcDesign;
    else % Coherence mode
        coherenceDesign = [];
        electrodePairCounter = 1;
        for mpfcElectrode = 1:numElectrodesMPFC
            for lpfcElectrode = 1:numElectrodesLPFC
                % Instead of creating separate variables, encode the pair as a single value
                tempDesign = baseDesign;
                tempDesign.ElectrodePair = repmat(electrodePairCounter, height(baseDesign), 1);
                coherenceDesign = [coherenceDesign; tempDesign];
                electrodePairCounter = electrodePairCounter + 1;
            end
        end
        
        designOutput1 = coherenceDesign;
        designOutput2 = []; % Not used in coherence mode
    end
end