function BlockPACAnalysis(SBJ, proc_id, placeholder1, placeholder2)

    % Set up path stuff for Cheaha
    restoredefaultpath;
    root_dir = '/data/user/anaskhan/';
    ft_dir = [root_dir 'fieldtrip/'];
    
    addpath([root_dir 'PRJ_Stroop/scripts']);
    addpath([root_dir 'PRJ_Stroop/scripts/utils']);
    addpath(ft_dir);
    ft_defaults;
   
    %% Load Data
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

    load(strcat(SBJ_vars.dirs.preproc, SBJ, '_preproc_', proc_id, '.mat'));
    load(strcat(SBJ_vars.dirs.events, SBJ, '_trial_info_final.mat'));

    %% Select for Medial and Lateral PFC channels
    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon, SBJ, '_elec_', proc_id, '_pat_', atlas_id, '_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    LPFC = intersect(elec.label(ismember(elec.ROI, 'dlPFC')), w2{1,1}.label(w2{1,1}.sig_chans));
    MPFC = intersect(elec.label(ismember(elec.ROI, 'aMCC') | ismember(elec.ROI, 'SMC')), w2{1,1}.label(w2{1,1}.sig_chans));

    %% Bring all subjects' fs down to 500
    origsamp = data.fsample;

    if origsamp ~= 500
        cfg_resamp = [];
        cfg_resamp.resamplefs = 500;
    
        data = ft_resampledata(cfg_resamp, data);
        trial_info.sample_rate = data.fsample;
        if (origsamp / data.fsample) == 2 
            trial_info.word_onset = round(trial_info.word_onset ./ (origsamp / data.fsample));
            trial_info.resp_onset = round(trial_info.resp_onset ./ (origsamp / data.fsample));
            trial_info.sample_rate = data.fsample;
            if isfield(trial_info, 'run_len')
                trial_info.run_len = round(trial_info.run_len ./ (origsamp / data.fsample));
            end
        else
            error('Discrepancies in sampling rate. Check again.');
        end
    end

    %% Select channels
    cfgs = [];
    cfgs.channel = [MPFC; LPFC];

    roi = ft_selectdata(cfgs, data);

    %% Cut into Trials
    trial_lim_s_pad = [-2.5 1.3]; % 1.5s buffers on each side

    events = trial_info.word_onset;
   
    roi_trl = fn_ft_cut_trials_equal_len(roi, events, trial_info.condition_n', ...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)] * roi.fsample));

    %% Get data for dlPFC and dmPFC
    data_all = cat(3, roi_trl.trial{:});
    mpfc_idx = ismember(roi_trl.label, MPFC);
    mpfc_data = data_all(mpfc_idx, :, :);
    lpfc_idx = ismember(roi_trl.label, LPFC);
    lpfc_data = data_all(lpfc_idx, :, :);
    clear data_all

    %% Create design matrix for LMM
    trials = fn_create_LMM_design(trial_info, 1, 1, 'power');
    
    %% Create mpfc and lpfc variables for the 3 block types: 'same', 'minc', 'mcon'
    block_types = {'minc', 'mcon'};
    
    % Initialize cell arrays to store data for each block type
    mpfc_block_data = cell(length(block_types), 1);
    lpfc_block_data = cell(length(block_types), 1);
    
    for b = 1:length(block_types)
        block_type = block_types{b};
        
        % Find trials corresponding to this block type
        block_trials = find(strcmpi(trials.BlockType, block_type) & strcmpi(trials.PreviousType, 'neu'));
        
        % Extract data for these trials
        mpfc_block_data{b} = mpfc_data(:, :, block_trials);
        lpfc_block_data{b} = lpfc_data(:, :, block_trials);
        
    end

    % Balancing trials across block types
       
    min_trials = min(cellfun(@(x) size(x, 3), mpfc_block_data));
    for b = 1:length(block_types)
        if size(mpfc_block_data{b}, 3) > min_trials
            rand_idx = randsample(size(mpfc_block_data{b}, 3), min_trials);
            mpfc_block_data{b} = mpfc_block_data{b}(:, :, rand_idx);
            lpfc_block_data{b} = lpfc_block_data{b}(:, :, rand_idx);
        end
    end

    % Assign variables for each block type
   
    % 'minc' block type
    mpfc_minc = mpfc_block_data{strcmp(block_types, 'minc')};
    lpfc_minc = lpfc_block_data{strcmp(block_types, 'minc')};

    % 'mcon' block type
    mpfc_mcon = mpfc_block_data{strcmp(block_types, 'mcon')};
    lpfc_mcon = lpfc_block_data{strcmp(block_types, 'mcon')};

    %% Start PAC computations
    n_surrogates = 200; % 200 by default
    n_bins = 18;
    LF_steps = 2:2:12;
    LF_bw = 2;
    HF_steps = 30:5:150;
    tcutEdge = 1.5;

    %% Initialize variables to store results
    PACLtM = [];  % Cross-region LPFC phase -> MPFC amplitude
    PACMtL = [];  % Cross-region MPFC phase -> LPFC amplitude

    %% Cross-region PAC
    if ~isempty(LPFC) && ~isempty(MPFC)
        % For each block type, compute cross-region PAC between all pairs of LPFC and MPFC channels
        
        % Initialize cell arrays to store results
        PACLtM = cell(length(block_types), 1);  % LPFC phase -> MPFC amplitude
        PACMtL = cell(length(block_types), 1);  % MPFC phase -> LPFC amplitude
        
        for idx = 1:length(block_types)
            block = block_types{idx};
            
            % Access data for the current block type
            lpfc_block = eval(['lpfc_' block]);  % Dimensions: [num_LPFC_channels x time x trials]
            mpfc_block = eval(['mpfc_' block]);  % Dimensions: [num_MPFC_channels x time x trials]
            
            num_lpfc_ch = size(lpfc_block, 1);
            num_mpfc_ch = size(mpfc_block, 1);
            
            % Initialize storage for results
            PACLtM_block = cell(num_lpfc_ch*num_mpfc_ch,1);  % For LPFC phase -> MPFC amplitude
            PACMtL_block = cell(num_lpfc_ch*num_mpfc_ch,1);  % For MPFC phase -> LPFC amplitude
            
            cnt = 1;
            % Loop over all pairs of LPFC and MPFC channels
            for lpfc_ch = 1:num_lpfc_ch
                for mpfc_ch = 1:num_mpfc_ch

                    % LPFC channel provides phase, MPFC channel provides amplitude
                    [comdlgrm, comdlgrm_z, phase2power] = crossROI_cfc_tort_comodulogram(...
                        squeeze(lpfc_block(lpfc_ch, :, :)), squeeze(mpfc_block(mpfc_ch, :, :)), ...
                        roi_trl.fsample, n_surrogates, n_bins, LF_steps, LF_bw, HF_steps, tcutEdge);
                    
                    PACLtM_block{cnt,1}.label = [LPFC{lpfc_ch} '_' MPFC{mpfc_ch}];
                    PACLtM_block{cnt,1}.zMI = comdlgrm_z;
                    PACLtM_block{cnt,1}.MI = comdlgrm;
                    PACLtM_block{cnt,1}.PPHist = phase2power;
                    PACLtM_block{cnt,1}.labelinfo = 'dlPFC_dmPFC';
                    PACLtM_block{cnt,1}.modDirection = "dlPFCphase->dmPFCamplitude";
                    PACLtM_block{cnt,1}.Condition = block;
                    
                    % MPFC channel provides phase, LPFC channel provides amplitude
                    [comdlgrm, comdlgrm_z, phase2power] = crossROI_cfc_tort_comodulogram(...
                        squeeze(mpfc_block(mpfc_ch, :, :)), squeeze(lpfc_block(lpfc_ch, :, :)), ...
                        roi_trl.fsample, n_surrogates, n_bins, LF_steps, LF_bw, HF_steps, tcutEdge);
                    
                    PACMtL_block{cnt,1}.label = [MPFC{mpfc_ch} '_' LPFC{lpfc_ch}];
                    PACMtL_block{cnt,1}.zMI = comdlgrm_z;
                    PACMtL_block{cnt,1}.MI = comdlgrm;
                    PACMtL_block{cnt,1}.PPHist = phase2power;
                    PACMtL_block{cnt,1}.labelinfo = 'dmPFC_dlPFC';
                    PACMtL_block{cnt,1}.modDirection = "dmPFCphase->dlPFCamplitude";
                    PACMtL_block{cnt,1}.Condition = block;

                    cnt = cnt + 1;
                end
            end
            
            % Store the results for the current block type
            PACLtM{idx} = PACLtM_block;
            PACMtL{idx} = PACMtL_block;
        end
    end

    %% Save the results for each subject
    save_dir = strcat(root_dir, 'PRJ_Stroop/results/PAC/');

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    psave_filename = fullfile(save_dir, [SBJ '.mat']);

    % Prepare variables to save
    save_vars = {};

    if exist('PACM', 'var') && ~isempty(PACM)
        save_vars{end+1} = 'PACM';
    end

    if exist('PACL', 'var') && ~isempty(PACL)
        save_vars{end+1} = 'PACL';
    end

    if exist('PACLtM', 'var') && ~isempty(PACLtM)
        save_vars{end+1} = 'PACLtM';
    end

    if exist('PACMtL', 'var') && ~isempty(PACMtL)
        save_vars{end+1} = 'PACMtL';
    end

    % Save only if there are variables to save
    if ~isempty(save_vars)
        save(psave_filename, save_vars{:}, '-v7.3');
    else
        warning('No PAC results to save for subject %s.', SBJ);
    end

end