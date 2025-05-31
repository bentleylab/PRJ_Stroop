function BlockPACAnalysisMPFC(SBJ, proc_id, placeholder1, placeholder2)

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
    MPFC = intersect(elec.label(ismember(elec.ROI, 'aMCC') | ismember(elec.ROI, 'SMC')), w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(MPFC)
        error('No MPFC channels found.')
    end
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
    cfgs.channel = MPFC;

    roi = ft_selectdata(cfgs, data);

    %% Cut into Trials
    trial_lim_s_pad = [-2.5 1.3]; % 1.5s buffers on each side

    events = trial_info.word_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi, events, trial_info.condition_n', ...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)] * roi.fsample));

    %% Extract MPFC data
    data_all = cat(3, roi_trl.trial{:});
    mpfc_idx = ismember(roi_trl.label, MPFC);
    mpfc_data = data_all(mpfc_idx, :, :);
    clear data_all

    %% Start PAC computations for MPFC
    n_surrogates = 200; % 200 by default
    n_bins = 18;
    LF_steps = 2:2:12;
    LF_bw = 2;
    HF_steps = 30:5:150;
    tcutEdge = 1.5;

    % Initialize cell array to store MPFC PAC results
    PACM = cell(length(MPFC), 1);

    % Loop through each MPFC channel
    for mpfc_ch = 1:length(MPFC)
        % Compute within-region PAC for MPFC
        [comdlgrm, comdlgrm_z, phase2power] = cfc_tort_comodulogram(...
            squeeze(mpfc_data(mpfc_ch, :, :)), roi_trl.fsample, n_surrogates, n_bins, LF_steps, LF_bw, HF_steps, tcutEdge);

        % Store results
        PACM{mpfc_ch}.label = MPFC{mpfc_ch};
        PACM{mpfc_ch}.zMI = comdlgrm_z;
        PACM{mpfc_ch}.MI = comdlgrm;
        PACM{mpfc_ch}.PPHist = phase2power;
        PACM{mpfc_ch}.labelinfo = 'dmPFC';
        PACM{mpfc_ch}.modDirection = "dmPFCphase->dmPFCamplitude";
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

    % Save only if there are variables to save
    if ~isempty(save_vars)
        save(psave_filename, save_vars{:}, '-v7.3');
    else
        warning('No PAC results to save for subject %s.', SBJ);
    end

end