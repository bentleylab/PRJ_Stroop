function PLVPACAnalysis(SBJ, proc_id, placeholder1,placeholder2)

    % Set up path stuff for Cheaha
    restoredefaultpath;
    root_dir = '/data/user/anaskhan/';
    ft_dir=[root_dir 'fieldtrip/'];
    
    addpath([root_dir 'PRJ_Stroop/scripts']);
    addpath([root_dir 'PRJ_Stroop/scripts/utils']);
    addpath(ft_dir);
    ft_defaults;
    
    %% Load Data and vars
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])
    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    
    %% Remove bad channels from data
    cfg = [];
    cfg.channel = SBJ_vars.ch_lab.ROI;
    data = ft_selectdata(cfg,data);
    
    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    
    LPFC = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));
    MPFC = intersect(elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC')),w2{1,1}.label(w2{1,1}.sig_chans));
    
    origsamp = data.fsample;
    %% Bring all subjects' fs down to 500
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
    %% Get task-active electrodes
    trial_lim_s_pad = [-0.5 - 1.5 1.25 + 1.5];
    
    stim_events = trial_info.word_onset;
    
    roi_trl = fn_ft_cut_trials_equal_len(data,stim_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*data.fsample));
    
    hg_vec  = 75:10:145;      
    win     = 5;              
    bl_start = -0.5;          
    bl_end   = -0.2;          
    niter    = 1000;   
    
    for f = 1:length(hg_vec)
        
        cfg            = [];
        cfg.bpfilter   = 'yes';
        cfg.bpfreq     = [hg_vec(f)-win, hg_vec(f)+win];  % e.g., [70 80], [80 90], ...
        cfg.hilbert    = 'abs';  % amplitude envelope
        tmpHFA         = ft_preprocessing(cfg, roi_trl);
        
        tmpHFA = zbaseline_HFA(tmpHFA, bl_start, bl_end, niter);
        
        hgtrace(f,:,:,:) = permute(cat(3,tmpHFA.trial{:}),[3 1 2]);
    end
    
    hgtrace = squeeze(mean(hgtrace));
    
    for itrial = 1:size(hgtrace,1)
        hfa.trial{1,itrial} = squeeze(hgtrace(itrial,:,:));
    end
    
    hfa.label = tmpHFA.label;
    hfa.time = tmpHFA.time;
    
    lowFreqs = 3:12;
    
    for f = 1:length(lowFreqs)
        cfg              = [];
        cfg.bpfilter     = 'yes';
        cfg.bpfreq       = [lowFreqs(f)-1 lowFreqs(f)+1];
        cfg.hilbert      = 'angle';
        theta_raw        = ft_preprocessing(cfg, roi_trl);
        
        cfg_trim         = [];
        cfg_trim.latency = [-0.5 1.25];
        theta_raw        = ft_selectdata(cfg_trim, theta_raw);
    
        theta_hfa        = ft_preprocessing(cfg, hfa);
        theta_hfa        = ft_selectdata(cfg_trim, theta_hfa);
    
        % Convert to [trials x channels x time]
        theta_raw_data = permute(cat(3, theta_raw.trial{:}), [3, 1, 2]);
        theta_hfa_data = permute(cat(3, theta_hfa.trial{:}), [3, 1, 2]);
        
        % Get channel labels
        labels = theta_raw.label;
        time = theta_raw.time{1,1};
        clear theta_raw theta_hfa
    
        %% ---- Balance conflict vs. non-conflict trials ----
    
        % Identify conflict vs. non-conflict trials
        conflictTrials   = strcmpi(trial_info.trialtype,'inc');
    
        nConflict    = sum(conflictTrials);
        nNoConflict  = sum(~conflictTrials);
    
        if nConflict < nNoConflict
            % Keep all conflict, randomly select the same number of non-conflict
            idxConflict       = find(conflictTrials);
            idxNonConflict    = find(~conflictTrials);
            
            % Randomly sample from the non-conflict trials
            idxNonConflictSub = randsample(idxNonConflict, nConflict);
            
            finalTrialIdx     = [idxConflict; idxNonConflictSub];
        else
            % Keep all non-conflict, randomly select the same number of conflict
            idxConflict       = find(conflictTrials);
            idxNonConflict    = find(~conflictTrials);
            
            % Randomly sample from the conflict trials
            idxConflictSub    = randsample(idxConflict, nNoConflict);
            
            finalTrialIdx     = [idxNonConflict; idxConflictSub];
        end
        
        % Subset the data so conflict and non-conflict are balanced
        theta_raw_data = theta_raw_data(finalTrialIdx, :, :);
        theta_hfa_data = theta_hfa_data(finalTrialIdx, :, :);
    
        balanceConflictMask = conflictTrials(finalTrialIdx);
    
        for mpfci = 1:length(MPFC)
            for lpfci = 1:length(LPFC)
    
                % Identify channel indices in the data's label list
                mpfc_chan_idx = strcmpi(labels, MPFC{mpfci});
                lpfc_chan_idx = strcmpi(labels, LPFC{lpfci});
        
                % Conflict trials
                MPFC_phase = theta_raw_data(balanceConflictMask, mpfc_chan_idx, :);
                LPFC_phase = theta_hfa_data(balanceConflictMask, lpfc_chan_idx, :);
        
                
                PLV_con = squeeze(abs(mean(exp(1i*(MPFC_phase - LPFC_phase)), 1)));
        
                % Non-conflict trials
                MPFC_phase = theta_raw_data(~balanceConflictMask, mpfc_chan_idx, :);
                LPFC_phase = theta_hfa_data(~balanceConflictMask, lpfc_chan_idx, :);
        
                PLV_ncon = squeeze(abs(mean(exp(1i*(MPFC_phase - LPFC_phase)), 1)));
        
                % Store
               
                conflict(f,mpfci,lpfci,:)  = PLV_con;
                noconflict(f,mpfci,lpfci,:) = PLV_ncon;
    
            end % lpfci
        end % mpfci
    end % f  
    
    % Baseline correction
    baseline_idx = time >= -0.5 & time <= -0.2;
    
    BL_PAC = mean(cat(4,conflict(:,:,:,baseline_idx),noconflict(:,:,:,baseline_idx)),4);
    
    conflict = ( (conflict - BL_PAC) ./ BL_PAC ) * 100;
    noconflict = ( (noconflict - BL_PAC) ./ BL_PAC ) * 100;
    
    conflict = squeeze(mean(mean(conflict,3),2));
    noconflict = squeeze(mean(mean(noconflict,3),2));
    
    time2save = time >= -0.1;

    time = time(time2save);

    conflict = conflict(:,time2save);
    noconflict = noconflict(:,time2save);

    save_dir = strcat(root_dir, 'PRJ_Stroop/results/PLVPAC/Stim/');
        
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    save_filename = fullfile(save_dir, [SBJ '.mat']);

    save(save_filename,'lowFreqs','time','MPFC','LPFC','conflict','noconflict','-v7.3');

end