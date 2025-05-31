function SequencePACAnalysis(SBJ, proc_id, an_id,cond)

    % Set up path stuff for Cheaha
    restoredefaultpath;
    root_dir = '/data/user/anaskhan/';
    ft_dir=[root_dir 'fieldtrip/'];

    addpath([root_dir 'PRJ_Stroop/scripts']);
    addpath([root_dir 'PRJ_Stroop/scripts/utils']);
    addpath(ft_dir);
    ft_defaults;
    startup;
    %% Load Data
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

    %% Select Conditions of Interest
    [cond_lab, ~, ~] = fn_condition_label_styles(cond);
    cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
    for cond_ix = 1:length(cond_lab)
        % Get binary condition index
        cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
            trial_info.condition_n));
    end

    %% Select for Medial and Lateral PFC channels
    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);
    
    MPFC = fn_select_elec_lab_match(elec, 'b', 'Dx', 'MPFC');
    LPFC = fn_select_elec_lab_match(elec, 'b', 'Dx', 'LPFC');
    
    % MPFC = elec.label(ismember(elec.ROI,{'ACC','aMCC'}));
    % LPFC = elec.label(ismember(elec.ROI,'dlPFC'));

    if isempty(MPFC) || isempty(LPFC)
    error('No MPFC or LPFC channels found.');
    end

    cfgs.channel = [MPFC;LPFC];
    roi = ft_selectdata(cfgs,data);
    origsamp = roi.fsample;
               
    %% Get task-active electrodes
    hfa = extractHFA(roi,-0.5,-0.2,1000,trial_info); % extract 70-150 Hz HFA and z-score to -0.5:-0.2s using boostrapping
    
    % Get only congruent and incongruent trials (ignore neutrals)
    CI_trials = logical(sum(cond_idx,1));

    % Start counting at stimulus onset
    start_trial_idx = dsearchn(hfa.time{1}',0);

    % Get mean over all trials (con + inc) for each channel
    hfa.meanHFA = squeeze(mean(hfa.ztrial(CI_trials,:,start_trial_idx:end),1));

    % Two-way test against baseline
    meanHFA_thresh = abs(hfa.meanHFA) >= 1.96;

    % Initialize channels matrix
    actv_chans = false(size(meanHFA_thresh,1),1);

    for chani = 1:size(meanHFA_thresh,1)
        islands = bwconncomp(meanHFA_thresh(chani,:));
        if islands.NumObjects == 0
            continue
        else
            if any(cellfun(@numel,islands.PixelIdxList) >= 0.1*roi.fsample) % > 100 ms
                actv_chans(chani,1) = true;
            end
        end        
    end

     % Logical indices for MPFC and LPFC channels within roi.label
    mpfc_indices = ismember(roi.label, MPFC);
    lpfc_indices = ismember(roi.label, LPFC);
    
    % Check if there's at least one active channel in MPFC and LPFC
    mpfc_active = any(actv_chans(mpfc_indices));
    lpfc_active = any(actv_chans(lpfc_indices));
    
    if ~(mpfc_active && lpfc_active)
        error(['Subject ' SBJ ' does not contain active channels in both LPFC and MPFC']);
    end

    % Select for active MPFC and LPFC channels
    aMPFC = roi.label(mpfc_indices & actv_chans);
    aLPFC = roi.label(lpfc_indices & actv_chans);

    cfgs.channel = [aMPFC;aLPFC];
    roi = ft_selectdata(cfgs,roi);
     %% Cut into Trials
    trial_lim_s_pad = [-0.75 2.25];

    stim_events = trial_info.word_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));
   
    %% Bring all subjects' fs down to 500 Hz
    
    if origsamp ~= 500
        cfg_resamp = [];
        cfg_resamp.resamplefs = 500;
        cfg_resamp.method = 'resample';
    
    
        roi_trl = ft_resampledata(cfg_resamp,roi_trl);
        trial_info.sample_rate = roi_trl.fsample;
    
        trial_info.word_onset = round(trial_info.word_onset./(origsamp/roi.fsample));
        trial_info.resp_onset = round(trial_info.resp_onset./(origsamp/roi.fsample));
    
        if isfield(trial_info,'run_len')
            trial_info.run_len = round(trial_info.run_len./(origsamp/roi.fsample));
        end
    end
    %% Create Deisgn Matrix and filter trials
    
    [coherenceDesign,~,isNeuCurrentOrPrevious] = fn_create_LMM_design(trial_info, length(aLPFC), length(aMPFC),'power');
    coherenceDesign(coherenceDesign.Electrode ~= 1,:) = [];

    allData = cat(3, roi_trl.trial{:});
   
    errors = logical(trial_info.error(~isNeuCurrentOrPrevious));

    allData = allData(:,:,~isNeuCurrentOrPrevious);
    
    allData = allData(:,:,~errors);

    % Reorder channels in X so that all MPFC come first followed by LPFC

    new_order_channels = [aMPFC; aLPFC];

    new_order_indices = zeros(length(new_order_channels), 1);

    % Loop over the new order of channel names to find their current indices
    for i = 1:length(new_order_channels)
        % Find the index of each channel name in the current ordering
        idx = find(strcmp(roi_trl.label, new_order_channels{i}));
        new_order_indices(i) = idx;
    end

    allData = allData(new_order_indices,:,:);

    X_cI = allData(:,:,strcmpi('inc',coherenceDesign.CurrentType) & strcmpi('con',coherenceDesign.PreviousType));
    X_iI = allData(:,:,strcmpi('inc',coherenceDesign.CurrentType) & strcmpi('inc',coherenceDesign.PreviousType));


    %% Compute PAC      
    fprintf('------------- Granger Calculations ----------\n');
    
    for mchani = 1:n_MPFC_channels
        for lchani = 1:n_LPFC_channels
            mpfc_index = mchani; % MPFC channel index
            lpfc_index = n_MPFC_channels + lchani;

            

            


    %% Baseline Correction
    fprintf('===================================================\n');
    fprintf('------------- Baseline Correction ----------\n');
    fprintf('===================================================\n');
  
    % Baseline indices calculation
    % 
    % baseline_start = dsearchn(grangerTime',-0.5); % -0.5s in index
    % baseline_end = dsearchn(grangerTime',-0.2); % -0.2s in index
    % 
    % Given: baseline_start, baseline_end, GC_cI, GC_iI
    % 
    % Extract Baseline Periods for Both Conditions
    % baseline_GC_cI = m_to_l_GC_cI(:,:,:, baseline_start:baseline_end);
    % baseline_GC_iI = m_to_l_GC_iI(:,:,:, baseline_start:baseline_end);
    % 
    % Concatenate along the time dimension (2nd dimension)
    % pooled_baseline = cat(4, baseline_GC_cI, baseline_GC_iI);
    % 
    % Step 3: Calculate Combined Mean Baseline
    % Since we concatenated along the 2nd dimension, the new size of the 2nd dimension is doubled
    % Compute mean across this combined baseline period (now the 2nd dimension)
    % mean_baseline_GC = mean(pooled_baseline, 4); % Mean over the combined time dimension
    % 
    % Step 4: Perform Percent Change Normalization on GC_cI and GC_iI using the combined baseline
    % Initialize matrices for normalized GC values
    % GC_cI_norm = 100*((m_to_l_GC_cI - mean_baseline_GC)./mean_baseline_GC);
    % GC_iI_norm = 100*((m_to_l_GC_iI - mean_baseline_GC)./mean_baseline_GC);
    % 
    % GC_cI_norm = GC_cI_norm(:,:,:,baseline_end:end);
    % GC_iI_norm = GC_iI_norm(:,:,:,baseline_end:end);
    % grangerTime = grangerTime(baseline_end:end);

    %% Save the results for each subject
    save_dir = strcat(root_dir, ['PRJ_Stroop/results/PGranger/',SBJ,'/']);

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    NPG_save_filename = fullfile(save_dir, [SBJ '.mat']);

    save(NPG_save_filename,'frex','grangerTime','l_to_m_GC_cI','l_to_m_GC_iI','MPFC','LPFC','aMPFC','aLPFC','-v7.3');

end