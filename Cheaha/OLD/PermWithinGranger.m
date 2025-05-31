function NonParGrangerTrialAnalysis(SBJ, proc_id, an_id, cond)

    %% Set up paths and load data
    restoredefaultpath;
    root_dir = '/data/user/anaskhan/';
    ft_dir = [root_dir 'fieldtrip/'];
    addpath([root_dir 'PRJ_Stroop/scripts']);
    addpath([root_dir 'PRJ_Stroop/scripts/utils']);
    addpath(ft_dir);
    ft_defaults;

    % Load subject-specific variables
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    
    % Load preprocessed data and trial information
    load(strcat(SBJ_vars.dirs.preproc, SBJ, '_preproc_', proc_id, '.mat'));
    load(strcat(SBJ_vars.dirs.events, SBJ, '_trial_info_final.mat'));

    % Load channels of interest (Medial and Lateral PFC)
    load(['/data/user/anaskhan/PRJ_Stroop/results/ThetaISPC/' SBJ '.mat'], 'maxChan');

    % Select channels of interest
    cfgs.channel = maxChan';
    roi = ft_selectdata(cfgs, data);
    origsamp = roi.fsample;

    %% Resample data if necessary
    if origsamp ~= 500
        cfg_resamp = [];
        cfg_resamp.resamplefs = 500;
        roi = ft_resampledata(cfg_resamp, roi);
        trial_info.sample_rate = roi.fsample;
        
        % Adjust event timings based on new sample rate
        if (origsamp / roi.fsample) == 2
            trial_info.word_onset = round(trial_info.word_onset / (origsamp / roi.fsample));
            trial_info.resp_onset = round(trial_info.resp_onset / (origsamp / roi.fsample));
            trial_info.sample_rate = roi.fsample;
            if isfield(trial_info, 'run_len')
                trial_info.run_len = round(trial_info.run_len / (origsamp / roi.fsample));
            end
        else
            error('Discrepancies in sampling rate. Check again.');
        end
    end

    %% Cut data into trials
    trial_lim_s_pad = [0 0.25];  % Define time window around events
    stim_events = trial_info.word_onset;  % Define events of interest
    roi_trl = fn_ft_cut_trials_equal_len(roi, stim_events, trial_info.condition_n', round([trial_lim_s_pad(1) trial_lim_s_pad(2)] * roi.fsample));

    %% Permutations
    n_permutations = 1000;  % Number of permutations
    perm_results = zeros(n_permutations, 1);  % Store permutation results

   
    comparisons = {'conflict', 'Gratton'};

    for perm = 1:n_permutations
        fprintf('Permutation %d/%d\n', perm, n_permutations);
        
        % Shuffle condition labels for conflict comparison
        shuffled_conflict_labels = trial_info.trialtype(randperm(length(trial_info.trialtype)));

        % Shuffle condition labels for Gratton comparison
        shuffled_gratton_labels = trial_info.PreviousType(randperm(length(trial_info.PreviousType)));

        % Initialize variables to store Granger causality results
        granger_con = []; granger_inc = []; granger_cC = []; granger_iC = [];
        
        for compi = 1:length(comparisons)
            if strcmpi(comparisons{compi}, 'conflict')
                % Conflict comparison (congruent vs. incongruent)
                % Determine the number of trials in each condition
                nCon = sum(strcmpi(shuffled_conflict_labels, 'con'));
                nInc = sum(strcmpi(shuffled_conflict_labels, 'inc'));
                conTrials = find(strcmpi(shuffled_conflict_labels, 'con'));
                incTrials = find(strcmpi(shuffled_conflict_labels, 'inc'));
                conIsSmaller = nCon < nInc;

                % Subsample trials to match the smaller condition
                if conIsSmaller
                    subInc = randsample(incTrials, nCon);
                    cfg_select = []; cfg_select.trials = subInc;
                    roi_trl_inc = ft_selectdata(cfg_select, roi_trl);
                    cfg_select.trials = conTrials;
                    roi_trl_con = ft_selectdata(cfg_select, roi_trl);
                elseif ~conIsSmaller
                    subCon = randsample(conTrials, nInc);
                    cfg_select = []; cfg_select.trials = subCon;
                    roi_trl_con = ft_selectdata(cfg_select, roi_trl);
                    cfg_select.trials = incTrials;
                    roi_trl_inc = ft_selectdata(cfg_select, roi_trl);
                else
                    cfg_select.trials = conTrials;
                    roi_trl_con = ft_selectdata(cfg_select, roi_trl);
                    cfg_select.trials = incTrials;
                    roi_trl_inc = ft_selectdata(cfg_select, roi_trl);
                end

                % Perform Fourier transform for frequency analysis
                cfg_spec = [];
                cfg_spec.method = 'mtmfft';
                cfg_spec.taper = 'hanning';
                cfg_spec.output = 'fourier';
                cfg_spec.pad = 1;
                freq_con = ft_freqanalysis(cfg_spec, roi_trl_con);
                freq_inc = ft_freqanalysis(cfg_spec, roi_trl_inc);

                % Compute Granger causality
                cfg = [];
                cfg.method = 'granger';
                cfg.granger.sfmethod = 'bivariate';
                cfg.granger.conditional = 'no';
                granger_con = ft_connectivityanalysis(cfg, freq_con);
                granger_inc = ft_connectivityanalysis(cfg, freq_inc);

            elseif strcmpi(comparisons{compi}, 'Gratton')
                % Gratton comparison (iC vs. cC)
                % Define trials for iC (incongruent followed by congruent) and cC (congruent followed by congruent)
                trials = fn_create_LMM_design(trial_info, 1, 1, 'power');
                niC = sum(strcmpi(trials.CurrentType, 'con') & strcmpi(shuffled_gratton_labels, 'inc'));
                ncC = sum(strcmpi(trials.CurrentType, 'con') & strcmpi(shuffled_gratton_labels, 'con'));
                iCTrials = find(strcmpi(trials.CurrentType, 'con') & strcmpi(shuffled_gratton_labels, 'inc'));
                cCTrials = find(strcmpi(trials.CurrentType, 'con') & strcmpi(shuffled_gratton_labels, 'con'));
                cCIsSmaller = ncC < niC;

                % Subsample trials to match the smaller condition
                if cCIsSmaller
                    subiC = randsample(iCTrials, ncC);
                    cfg_select = []; cfg_select.trials = subiC;
                    roi_trl_iC = ft_selectdata(cfg_select, roi_trl);
                    cfg_select.trials = cCTrials;
                    roi_trl_cC = ft_selectdata(cfg_select, roi_trl);
                elseif ~cCIsSmaller
                    subcC = randsample(cCTrials, niC);
                    cfg_select = []; cfg_select.trials = subcC;
                    roi_trl_cC = ft_selectdata(cfg_select, roi_trl);
                    cfg_select.trials = iCTrials;
                    roi_trl_iC = ft_selectdata(cfg_select, roi_trl);
                else
                    cfg_select.trials = cCTrials;
                    roi_trl_cC = ft_selectdata(cfg_select, roi_trl);
                    cfg_select.trials = iCTrials;
                    roi_trl_iC = ft_selectdata(cfg_select, roi_trl);
                end

                % Perform Fourier transform for frequency analysis
                cfg_spec = [];
                cfg_spec.method = 'mtmfft';
                cfg_spec.taper = 'hanning';
                cfg_spec.output = 'fourier';
                cfg_spec.pad = 1;
                freq_cC = ft_freqanalysis(cfg_spec, roi_trl_cC);
                freq_iC = ft_freqanalysis(cfg_spec, roi_trl_iC);

                % Compute Granger causality
                cfg = [];
                cfg.method = 'granger';
                cfg.granger.sfmethod = 'bivariate';
                cfg.granger.conditional = 'no';
                granger_cC = ft_connectivityanalysis(cfg, freq_cC);
                granger_iC = ft_connectivityanalysis(cfg, freq_iC);
            end
        end

        % Step 6: Compute net LPFC -> MPFC Granger Causality and Test Statistic
        % Calculate the net GC by subtracting MPFC -> LPFC from LPFC -> MPFC
        net_GC_con = calculate_net_GC(granger_con, LPFC, MPFC);
        net_GC_inc = calculate_net_GC(granger_inc, LPFC, MPFC);
        net_GC_cC = calculate_net_GC(granger_cC, LPFC, MPFC);
        net_GC_iC = calculate_net_GC(granger_iC, LPFC, MPFC);

        % Define frequency range of interest (4-8 Hz)
        freq_idx = find(frex >= 4 & frex <= 8);

        % Average net GC over the 4-8 Hz range
        avg_net_GC_conflict = mean(net_GC_inc(freq_idx, :) - net_GC_con(freq_idx, :), 1);
        avg_net_GC_gratton = mean(net_GC_iC(freq_idx, :) - net_GC_cC(freq_idx, :), 1);

        % Compute the test statistic: (iC - cC) - (inc - con)
        test_stat = avg_net_GC_gratton - avg_net_GC_conflict;

        % Step 7: Store the test statistic for this permutation
        perm_results(perm) = test_stat;
    end

    % Step 8: Save the permutation results
    save_dir = strcat(root_dir, ['PRJ_Stroop/results/NPG/Stim/', SBJ, '/']);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    NPG_save_filename = fullfile(save_dir, [SBJ '_perm_results.mat']);
    save(NPG_save_filename, 'perm_results', '-v7.3');

end

function net_GC = calculate_net_GC(granger, LPFC, MPFC)
    % Helper function to calculate net LPFC -> MPFC Granger Causality
    net_GC = zeros(size(granger.freq));

    for idx = 1:size(granger.labelcmb, 1)
        from_label = strtok(granger.labelcmb{idx, 1}, '[');
        to_label = strtok(granger.labelcmb{idx, 2}, '[');

        if any(strcmp(from_label, LPFC)) && any(strcmp(to_label, MPFC))
            lpfc_to_mpfc_gc = granger.grangerspctrm(idx, :);
        elseif any(strcmp(from_label, MPFC)) && any(strcmp(to_label, LPFC))
            mpfc_to_lpfc_gc = granger.grangerspctrm(idx, :);
        end
    end

    net_GC = lpfc_to_mpfc_gc - mpfc_to_lpfc_gc;
end
