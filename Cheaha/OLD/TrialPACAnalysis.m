function TrialPACAnalysis(SBJ, proc_id, placeholder1, placeholder2)

    % Set up path stuff for Cheaha
    restoredefaultpath;
    root_dir = '/data/user/anaskhan/';
    ft_dir=[root_dir 'fieldtrip/'];
    
    addpath([root_dir 'PRJ_Stroop/scripts']);
    addpath([root_dir 'PRJ_Stroop/scripts/utils']);
    addpath(ft_dir)
    ft_defaults;
   
    % Load Data
    
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

    
    % Select for Medial and Lateral PFC channels
    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    LPFC = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));
    MPFC = intersect(elec.label(ismember(elec.ROI,'aMCC')),w2{1,1}.label(w2{1,1}.sig_chans));


    %% Bring all subjects' fs down to 500
    origsamp = data.fsample;

    if origsamp ~= 500
        cfg_resamp = [];
        cfg_resamp.resamplefs = 500;
    
        data = ft_resampledata(cfg_resamp,data);
        trial_info.sample_rate = data.fsample;
            if (origsamp/data.fsample) == 2 
                trial_info.word_onset = round(trial_info.word_onset./(origsamp/data.fsample));
                trial_info.resp_onset = round(trial_info.resp_onset./(origsamp/data.fsample));
                trial_info.sample_rate = data.fsample;
                if isfield(trial_info,'run_len')
                    trial_info.run_len = round(trial_info.run_len./(origsamp/data.fsample));
                end
            else
                error('Discrepancies in sampling rate. Check again.')
            end
    end


    % Select for channels
    cfgs = [];
    cfgs.channel = [MPFC;LPFC];

    roi = ft_selectdata(cfgs,data);


    % Cut into Trials
    trial_lim_s_pad = [-1.5 2]; % 1.5s buffers on each side

    resp_events = trial_info.resp_onset;
   
    roi_trl = fn_ft_cut_trials_equal_len(roi,resp_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));


    %% Start PAC computations
   
    n_surrogates = 200; % 200 by default
    n_bins = 18;
    LF_steps = 2:2:12;
    LF_bw = 2;
    HF_steps = 30:5:150;
    tcutEdge = 1.5;

    %% Get data for dlPFC and dACC
    data = cat(3,roi_trl.trial{:});
    mpfc_idx = ismember(roi_trl.label,MPFC);
    mpfc_data = data(mpfc_idx,:,:);
    lpfc_idx = ismember(roi_trl.label,LPFC);
    lpfc_data = data(lpfc_idx,:,:);
    clear data

    
    % Create design matrix for LMM (assuming the function and data are similar)
    trials = fn_create_LMM_design(trial_info,1,1,'power');
    
    % Create CurrentConflict column
    trials.CurrentConflict = repmat({'NoConflict'}, height(trials), 1); % Default to 'NoConflict'
    trials.CurrentConflict(strcmp(trials.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc' 
    
    % Determine which condition has fewer trials
    nConflict = sum(strcmpi(trials.CurrentConflict, 'Conflict'));
    nNoConflict = sum(strcmpi(trials.CurrentConflict, 'NoConflict'));
    
    % Get trial indices for each condition
    conflict_trials = find(strcmpi(trials.CurrentConflict, 'Conflict'));
    noconflict_trials = find(strcmpi(trials.CurrentConflict, 'NoConflict'));
    
    % Extract trials for each condition
    lpfc_Conflict = lpfc_data(:, :, conflict_trials);
    mpfc_Conflict = mpfc_data(:, :, conflict_trials);
    
    lpfc_NoConflict = lpfc_data(:, :, noconflict_trials);
    mpfc_NoConflict = mpfc_data(:, :, noconflict_trials);
    
    % Balance the number of trials by subsampling the larger set
    if nConflict > nNoConflict
    
        % Get sub-sample from "Conflict" trials equal to the number of "NoConflict" trials
        subConflictIdx = randsample(nConflict, nNoConflict);
        subConflict = conflict_trials(subConflictIdx);
    
        % Subsample data
        lpfc_Conflict = lpfc_data(:, :, subConflict);
        mpfc_Conflict = mpfc_data(:, :, subConflict);
    
    elseif nConflict < nNoConflict
    
        % Get sub-sample from "NoConflict" trials equal to the number of "Conflict" trials
        subNoConflictIdx = randsample(nNoConflict, nConflict);
        subNoConflict = noconflict_trials(subNoConflictIdx);
    
        % Subsample data
        lpfc_NoConflict = lpfc_data(:, :, subNoConflict);
        mpfc_NoConflict = mpfc_data(:, :, subNoConflict);
       
    end


    %% Cross-region
    % LPFC phase & MPFC amp

    % % cC
    % [comdlgrm, comdlgrm_z, phase2power] = crossROI_cfc_tort_comodulogram(lpfc_cC,mpfc_cC,roi_trl.fsample,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge);
    % PACLtM{1}.label = [LPFC_chan{1} '_' MPFC_chan{1}];
    % PACLtM{1}.zMI = comdlgrm_z;
    % PACLtM{1}.MI = comdlgrm;
    % PACLtM{1}.PPHist = phase2power;
    % PACLtM{1}.labelinfo = 'dlPFC_dACC';
    % PACLtM{1}.modDirection = "dlPFCphase->dACCamplitude";
    % PACLtM{1}.Condition = 'cC';
    % 
    % % iC
    % [comdlgrm, comdlgrm_z, phase2power] = crossROI_cfc_tort_comodulogram(lpfc_iC,mpfc_iC,roi_trl.fsample,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge);
    % PACLtM{2}.label = [LPFC_chan{1} '_' MPFC_chan{1}];
    % PACLtM{2}.zMI = comdlgrm_z;
    % PACLtM{2}.MI = comdlgrm;
    % PACLtM{2}.PPHist = phase2power;
    % PACLtM{2}.labelinfo = 'dlPFC_dACC';
    % PACLtM{2}.modDirection = "dlPFCphase->dACCamplitude";
    % PACLtM{2}.Condition = 'iC';
    % 
    % % MPFC phase & LPFC amp
    % 
    % % cC
    % [comdlgrm, comdlgrm_z, phase2power] = crossROI_cfc_tort_comodulogram(mpfc_cC,lpfc_cC,roi_trl.fsample,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge);
    % PACMtL{1}.label = [MPFC_chan{1} '_' LPFC_chan{1}];
    % PACMtL{1}.zMI = comdlgrm_z;
    % PACMtL{1}.MI = comdlgrm;
    % PACMtL{1}.PPHist = phase2power;
    % PACMtL{1}.labelinfo = 'dACC_dlPFC';
    % PACMtL{1}.modDirection = "dACCphase->dlPFCamplitude";
    % PACMtL{1}.Condition = 'cC';
    % 
    % % iC
    % [comdlgrm, comdlgrm_z, phase2power] = crossROI_cfc_tort_comodulogram(mpfc_iC,lpfc_iC,roi_trl.fsample,n_surrogates,n_bins,LF_steps,LF_bw,HF_steps,tcutEdge);
    % PACMtL{2}.label = [MPFC_chan{1} '_' LPFC_chan{1}];
    % PACMtL{2}.zMI = comdlgrm_z;
    % PACMtL{2}.MI = comdlgrm;
    % PACMtL{2}.PPHist = phase2power;
    % PACMtL{2}.labelinfo = 'dACC_dlPFC';
    % PACMtL{2}.modDirection = "dACCphase->dlPFCamplitude";
    % PACMtL{2}.Condition = 'iC';

    %% Within-region

   
    % MPFC
    for ch = 1:length(MPFC)
        % Condition: "con"
        [comdlgrm, comdlgrm_z, phase2power] = cfc_tort_comodulogram(squeeze(mpfc_con(ch, :, :)), roi_trl.fsample, n_surrogates, n_bins, LF_steps, LF_bw, HF_steps, tcutEdge);
        PACM{1}.label{ch,1} = MPFC{ch};
        PACM{1}.zMI{ch,1} = comdlgrm_z;
        PACM{1}.MI{ch,1} = comdlgrm;
        PACM{1}.PPHist{ch,1} = phase2power;
        PACM{1}.labelinfo = 'dACC';
        PACM{1}.Condition = 'con';
    
        % Condition: "inc"
        [comdlgrm, comdlgrm_z, phase2power] = cfc_tort_comodulogram(squeeze(mpfc_inc(ch, :, :)), roi_trl.fsample, n_surrogates, n_bins, LF_steps, LF_bw, HF_steps, tcutEdge);
        PACM{2}.label{ch,1} = MPFC{ch};
        PACM{2}.zMI{ch,1} = comdlgrm_z;
        PACM{2}.MI{ch,1} = comdlgrm;
        PACM{2}.PPHist{ch,1} = phase2power;
        PACM{2}.labelinfo = 'dACC';
        PACM{2}.Condition = 'inc';
    end
    
    % LPFC
    for ch = 1:length(LPFC)
        % Condition: "con"
        [comdlgrm, comdlgrm_z, phase2power] = cfc_tort_comodulogram(squeeze(lpfc_con(ch, :, :)), roi_trl.fsample, n_surrogates, n_bins, LF_steps, LF_bw, HF_steps, tcutEdge);
        PACL{1}.label{ch,1} = LPFC{ch};
        PACL{1}.zMI{ch,1} = comdlgrm_z;
        PACL{1}.MI{ch,1} = comdlgrm;
        PACL{1}.PPHist{ch,1} = phase2power;
        PACL{1}.labelinfo = 'dlPFC';
        PACL{1}.Condition = 'con';
    
        % Condition: "inc"
        [comdlgrm, comdlgrm_z, phase2power] = cfc_tort_comodulogram(squeeze(lpfc_inc(ch, :, :)), roi_trl.fsample, n_surrogates, n_bins, LF_steps, LF_bw, HF_steps, tcutEdge);
        PACL{2}.label{ch,1} = LPFC{ch};
        PACL{2}.zMI{ch,1} = comdlgrm_z;
        PACL{2}.MI{ch,1} = comdlgrm;
        PACL{2}.PPHist{ch,1} = phase2power;
        PACL{2}.labelinfo = 'dlPFC';
        PACL{2}.Condition = 'inc';
    end

    % Save the results for each subject
    save_dir = strcat(root_dir, 'PRJ_Stroop/results/PAC/');
    
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    psave_filename = fullfile(save_dir,[SBJ '.mat']);
    
    save(psave_filename,'PACM','PACL','-v7.3');
end