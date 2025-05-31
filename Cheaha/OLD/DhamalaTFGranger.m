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
   
    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    MPFC = intersect(elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC')),w2{1,1}.label(w2{1,1}.sig_chans));
    LPFC = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(LPFC) || isempty(MPFC)
        error('Subject does not have both LPFC and MPFC channels.')
    end

    cfgs.channel = [MPFC;LPFC];
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
    % trial_lim_s_pad = [-1.5 2.75];
    % 
    % stim_events = trial_info.word_onset;
    % 
    % roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
    %     round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));
    % 
    %% Cut into Trials (R)
    trial_lim_s_pad = [-1 1.5];

    stim_events = trial_info.resp_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));


    [trials,~] = fn_create_LMM_design(trial_info,1,1,'power');
    % Create CurrentConflict column
    trials.CurrentConflict = repmat({'NoConflict'}, height(trials), 1); % Default to 'NoConflict'
    trials.CurrentConflict(strcmp(trials.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc'

    % Create PreviousConflict column
    trials.PreviousConflict = repmat({'NoConflict'}, height(trials), 1); % Default to 'NoConflict'
    trials.PreviousConflict(strcmp(trials.PreviousType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where PreviousType is 'inc'

    freqrange = 20:150;
    
    % Determine whether congruent or incongruent is lower in
    % quantity
    nCon = sum(strcmpi(trials.CurrentConflict,'NoConflict'));
    nInc = sum(strcmpi(trials.CurrentConflict,'Conflict'));

    conTrials = find(strcmpi(trials.CurrentConflict,'NoConflict'));
    incTrials = find(strcmpi(trials.CurrentConflict,'Conflict'));

    conIsSmaller = nCon < nInc;

    if conIsSmaller

        % Get subsample of trials
        subInc = randsample(incTrials,nCon);

        % Incongruent
        cfg_select = [];
        cfg_select.trials = subInc;
        roi_trl_inc = ft_selectdata(cfg_select,roi_trl);

        % Congruent
        cfg_select = [];
        cfg_select.trials = conTrials;
        roi_trl_con = ft_selectdata(cfg_select,roi_trl);

    elseif ~conIsSmaller
        
        % Get subsample of trials
        subCon = randsample(conTrials,nInc);

        % Congruent
        cfg_select = [];
        cfg_select.trials = subCon;
        roi_trl_con = ft_selectdata(cfg_select,roi_trl);

        % Incongruent
        cfg_select = [];
        cfg_select.trials = incTrials;
        roi_trl_inc = ft_selectdata(cfg_select,roi_trl);

    else

        % Congruent
        cfg_select = [];
        cfg_select.trials = conTrials;
        roi_trl_con = ft_selectdata(cfg_select,roi_trl);

        % Incongruent
        cfg_select = [];
        cfg_select.trials = incTrials;
        roi_trl_inc = ft_selectdata(cfg_select,roi_trl);

    end

    % Compute NP GC 
    inc_data = cat(3,roi_trl_inc.trial{:});
    MPFC_incdata = inc_data(ismember(roi_trl_inc.label,MPFC),:,:);
    LPFC_incdata = inc_data(ismember(roi_trl_inc.label,LPFC),:,:);

    X_inc = cat(1,MPFC_incdata,LPFC_incdata);
    X_inc = permute(X_inc,[2 3 1]);

    con_data = cat(3,roi_trl_con.trial{:});
    MPFC_condata = con_data(ismember(roi_trl_inc.label,MPFC),:,:);
    LPFC_condata = con_data(ismember(roi_trl_inc.label,LPFC),:,:);

    X_con = cat(1,MPFC_condata,LPFC_condata);
    X_con = permute(X_con,[2 3 1]);

    [frex, granger_inc, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_inc, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));
    [~, granger_con, ~, ~, ~, ~] = compute_nonparGGC_wavelet_Morlet(X_con, 500, freqrange, 0, 0, (8*pi^2-1)/(4*pi));

    time = roi_trl_con.time{1};
    save_dir = strcat(root_dir, ['PRJ_Stroop/results/NPG/Resp/',SBJ,'/']);

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    NPG_save_filename = fullfile(save_dir, [SBJ '.mat']);
    save(NPG_save_filename,'granger_con', 'granger_inc','frex','time','MPFC','LPFC','-v7.3');

end