function DhamalaTFGranger(SBJ, proc_id, an_id,cond)

    % Set up path stuff for Cheaha
    restoredefaultpath;
    root_dir = '/data/user/anaskhan/';
    ft_dir=[root_dir 'fieldtrip/'];

    addpath([root_dir 'PRJ_Stroop/scripts']);
    addpath([root_dir 'PRJ_Stroop/scripts/utils']);
    addpath(ft_dir);
    ft_defaults;
   
    %% Load Data
    
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

   
    %% Select for Medial and Lateral PFC channels
   
    load(['/data/user/anaskhan/PRJ_Stroop/results/ThetaISPC/' SBJ '.mat'],'maxChan')

    cfgs.channel = maxChan';
    roi = ft_selectdata(cfgs,data);
    origsamp = roi.fsample;
    %% Bring all subjects' fs down to 500
    if origsamp ~= 500
        cfg_resamp = [];
        cfg_resamp.resamplefs = 500;
    
        roi = ft_resampledata(cfg_resamp,roi);
        trial_info.sample_rate = roi.fsample;
            if (origsamp/roi.fsample) == 2 
                trial_info.word_onset = round(trial_info.word_onset./(origsamp/roi.fsample));
                trial_info.resp_onset = round(trial_info.resp_onset./(origsamp/roi.fsample));
                trial_info.sample_rate = roi.fsample;
                if isfield(trial_info,'run_len')
                    trial_info.run_len = round(trial_info.run_len./(origsamp/roi.fsample));
                end
            else
                error('Discrepancies in sampling rate. Check again.')
            end
    end
    
    %% Cut into Trials (S)
    trial_lim_s_pad = [-1.5 2.75];

    stim_events = trial_info.word_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));
   
    %% Create Trials Table with Conflict Columns
    [trials,~] = fn_create_LMM_design(trial_info,1,1,'power');
    % Create CurrentConflict column
    trials.CurrentConflict = repmat({'NoConflict'}, height(trials), 1); % Default to 'NoConflict'
    trials.CurrentConflict(strcmp(trials.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc'

    % Create PreviousConflict column
    trials.PreviousConflict = repmat({'NoConflict'}, height(trials), 1); % Default to 'NoConflict'
    trials.PreviousConflict(strcmp(trials.PreviousType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where PreviousType is 'inc'

    freqrange = 2.^((10:50)/10);
    comparisons = {'conflict','GrattonC','BlockConflict'};
    for compi = 1:length(comparisons)
        if strcmpi(comparisons{compi},'conflict')

            % Determine whether Conflict or NoConflict is lower in quantity
            nConflict = sum(strcmpi(trials.CurrentConflict,'Conflict'));
            nNoConflict = sum(strcmpi(trials.CurrentConflict,'NoConflict'));

            ConflictTrials = find(strcmpi(trials.CurrentConflict,'Conflict'));
            NoConflictTrials = find(strcmpi(trials.CurrentConflict,'NoConflict'));

            ConflictIsSmaller = nConflict < nNoConflict;

            if ConflictIsSmaller

                % Get subsample of NoConflict trials
                subNoConflict = randsample(NoConflictTrials,nConflict);

                % Conflict
                cfg_select = [];
                cfg_select.trials = ConflictTrials;
                roi_trl_conflict = ft_selectdata(cfg_select,roi_trl);

                % NoConflict
                cfg_select = [];
                cfg_select.trials = subNoConflict;
                roi_trl_NoConflict = ft_selectdata(cfg_select,roi_trl);

            elseif ~ConflictIsSmaller
                
                % Get subsample of Conflict trials
                subConflict = randsample(ConflictTrials,nNoConflict);

                % NoConflict
                cfg_select = [];
                cfg_select.trials = NoConflictTrials;
                roi_trl_NoConflict = ft_selectdata(cfg_select,roi_trl);

                % Conflict
                cfg_select = [];
                cfg_select.trials = subConflict;
                roi_trl_conflict = ft_selectdata(cfg_select,roi_trl);

            else

                % Conflict
                cfg_select = [];
                cfg_select.trials = ConflictTrials;
                roi_trl_conflict = ft_selectdata(cfg_select,roi_trl);

                % NoConflict
                cfg_select = [];
                cfg_select.trials = NoConflictTrials;
                roi_trl_NoConflict = ft_selectdata(cfg_select,roi_trl);

            end

            % Compute NP GC 
            conflict_data = cat(3,roi_trl_conflict.trial{:});
            MPFC_conflictdata = conflict_data(strcmpi(roi_trl_conflict.label,maxChan(1)),:,:);
            LPFC_conflictdata = conflict_data(strcmpi(roi_trl_conflict.label,maxChan(2)),:,:);

            X_conflict = cat(1,MPFC_conflictdata,LPFC_conflictdata);
            X_conflict = permute(X_conflict,[2 3 1]);

            NoConflict_data = cat(3,roi_trl_NoConflict.trial{:});
            MPFC_NoConflictdata = NoConflict_data(strcmpi(roi_trl_NoConflict.label,maxChan(1)),:,:);
            LPFC_NoConflictdata = NoConflict_data(strcmpi(roi_trl_NoConflict.label,maxChan(2)),:,:);

            X_NoConflict = cat(1,MPFC_NoConflictdata,LPFC_NoConflictdata);
            X_NoConflict = permute(X_NoConflict,[2 3 1]);

            [frex, granger_conflict, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_conflict, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));
            [~, granger_NoConflict, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_NoConflict, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));


        elseif strcmpi(comparisons{compi},'GrattonC')

            % Determine whether iC or cC is lower in quantity

            niC = sum(strcmpi(trials.CurrentConflict,'NoConflict') & strcmpi(trials.PreviousConflict,'Conflict'));
            ncC = sum(strcmpi(trials.CurrentConflict,'NoConflict') & strcmpi(trials.PreviousConflict,'NoConflict'));


            iCTrials = find(strcmpi(trials.CurrentConflict,'NoConflict') & strcmpi(trials.PreviousConflict,'Conflict'));
            cCTrials = find(strcmpi(trials.CurrentConflict,'NoConflict') & strcmpi(trials.PreviousConflict,'NoConflict'));
            
            cCIsSmaller = ncC < niC;

            if cCIsSmaller

                % Get subsample of trials
                subiC = randsample(iCTrials,ncC);

                % iC
                cfg_select = [];
                cfg_select.trials = subiC;
                roi_trl_iC = ft_selectdata(cfg_select,roi_trl);

                % cC
                cfg_select = [];
                cfg_select.trials = cCTrials;
                roi_trl_cC = ft_selectdata(cfg_select,roi_trl);

            elseif ~cCIsSmaller
                
                % Get subsample of trials
                subcC = randsample(cCTrials,niC);

                % cC
                cfg_select = [];
                cfg_select.trials = subcC;
                roi_trl_cC = ft_selectdata(cfg_select,roi_trl);

                % iC
                cfg_select = [];
                cfg_select.trials = iCTrials;
                roi_trl_iC = ft_selectdata(cfg_select,roi_trl);

            else

                % cC
                cfg_select = [];
                cfg_select.trials = cCTrials;
                roi_trl_cC = ft_selectdata(cfg_select,roi_trl);

                % iC
                cfg_select = [];
                cfg_select.trials = iCTrials;
                roi_trl_iC = ft_selectdata(cfg_select,roi_trl);

            end

            % Compute NP GC     
            iC_data = cat(3,roi_trl_iC.trial{:});
            MPFC_iCdata = iC_data(strcmpi(roi_trl_iC.label,maxChan(1)),:,:);
            LPFC_iCdata = iC_data(strcmpi(roi_trl_iC.label,maxChan(2)),:,:);
            
            X_iC = cat(1,MPFC_iCdata,LPFC_iCdata);
            X_iC = permute(X_iC,[2 3 1]);

            cC_data = cat(3,roi_trl_cC.trial{:});
            MPFC_cCdata = cC_data(strcmpi(roi_trl_cC.label,maxChan(1)),:,:);
            LPFC_cCdata = cC_data(strcmpi(roi_trl_cC.label,maxChan(2)),:,:);

            X_cC = cat(1,MPFC_cCdata,LPFC_cCdata);
            X_cC = permute(X_cC,[2 3 1]);
            
            [frex, granger_iC, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_iC, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));
            [~, granger_cC, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_cC, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));

        elseif strcmpi(comparisons{compi},'BlockConflict')

            % Code block for conflict and noconflict trials separated by "same", "minc", and "mcon" blocks

            block_types = {'same', 'minc', 'mcon'};

            for blocki = 1:length(block_types)
                block_type = block_types{blocki};

                % Select trials in this block
                block_trials = strcmpi(trials.BlockType, block_type);

                % Now within this block, find Conflict and NoConflict trials
                nConflict = sum(block_trials & strcmpi(trials.CurrentConflict,'Conflict'));
                nNoConflict = sum(block_trials & strcmpi(trials.CurrentConflict,'NoConflict'));

                ConflictTrials = find(block_trials & strcmpi(trials.CurrentConflict,'Conflict'));
                NoConflictTrials = find(block_trials & strcmpi(trials.CurrentConflict,'NoConflict'));

                ConflictIsSmaller = nConflict < nNoConflict;

                if ConflictIsSmaller

                    % Get subsample of NoConflict trials
                    subNoConflict = randsample(NoConflictTrials,nConflict);

                    % Conflict
                    cfg_select = [];
                    cfg_select.trials = ConflictTrials;
                    roi_trl_conflict = ft_selectdata(cfg_select,roi_trl);

                    % NoConflict
                    cfg_select = [];
                    cfg_select.trials = subNoConflict;
                    roi_trl_NoConflict = ft_selectdata(cfg_select,roi_trl);

                elseif ~ConflictIsSmaller

                    % Get subsample of Conflict trials
                    subConflict = randsample(ConflictTrials,nNoConflict);

                    % NoConflict
                    cfg_select = [];
                    cfg_select.trials = NoConflictTrials;
                    roi_trl_NoConflict = ft_selectdata(cfg_select,roi_trl);

                    % Conflict
                    cfg_select = [];
                    cfg_select.trials = subConflict;
                    roi_trl_conflict = ft_selectdata(cfg_select,roi_trl);

                else

                    % Conflict
                    cfg_select = [];
                    cfg_select.trials = ConflictTrials;
                    roi_trl_conflict = ft_selectdata(cfg_select,roi_trl);

                    % NoConflict
                    cfg_select = [];
                    cfg_select.trials = NoConflictTrials;
                    roi_trl_NoConflict = ft_selectdata(cfg_select,roi_trl);

                end

                % Compute NP GC 
                conflict_data = cat(3,roi_trl_conflict.trial{:});
                MPFC_conflictdata = conflict_data(strcmpi(roi_trl_conflict.label,maxChan(1)),:,:);
                LPFC_conflictdata = conflict_data(strcmpi(roi_trl_conflict.label,maxChan(2)),:,:);

                X_conflict = cat(1,MPFC_conflictdata,LPFC_conflictdata);
                X_conflict = permute(X_conflict,[2 3 1]);

                NoConflict_data = cat(3,roi_trl_NoConflict.trial{:});
                MPFC_NoConflictdata = NoConflict_data(strcmpi(roi_trl_NoConflict.label,maxChan(1)),:,:);
                LPFC_NoConflictdata = NoConflict_data(strcmpi(roi_trl_NoConflict.label,maxChan(2)),:,:);

                X_NoConflict = cat(1,MPFC_NoConflictdata,LPFC_NoConflictdata);
                X_NoConflict = permute(X_NoConflict,[2 3 1]);

                [frex, granger_conflict_block, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_conflict, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));
                [~, granger_NoConflict_block, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_NoConflict, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));

                % Save or store the granger results for this block type
                % e.g., granger_conflict_same, granger_NoConflict_same, etc.
                eval(['granger_conflict_' block_type ' = granger_conflict_block;']);
                eval(['granger_NoConflict_' block_type ' = granger_NoConflict_block;']);

            end
        end
    end

    time = roi_trl.time{1};
    save_dir = strcat(root_dir, ['PRJ_Stroop/results/NPG/Stim/',SBJ,'/']);

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    NPG_save_filename = fullfile(save_dir, [SBJ '.mat']);
    save_vars = {'granger_iC','granger_cC',...
        'granger_conflict', 'granger_NoConflict',...
        'frex','time'};

    % Include block-specific granger variables
    block_types = {'same', 'minc', 'mcon'};
    for blocki = 1:length(block_types)
        block_type = block_types{blocki};
        save_vars{end+1} = ['granger_conflict_' block_type];
        save_vars{end+1} = ['granger_NoConflict_' block_type];
    end

    save(NPG_save_filename, save_vars{:}, '-v7.3');

end
