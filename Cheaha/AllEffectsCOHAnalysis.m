function AllEffectsCOHAnalysis(SBJ, proc_id, an_id,cond)

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
    eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));


    %% Select for Medial and Lateral PFC channels
    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    MPFC = intersect(elec.label(ismember(elec.ROI,'aMCC')),w2{1,1}.label(w2{1,1}.sig_chans));
    LPFC = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(MPFC) || isempty(LPFC)
        error('No MPFC or LPFC channels found.');
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
    
     %% Cut into Trials
    trial_lim_s_pad = [an.trial_lim_s(1)-(1/min(cfg_tfr.foi))*3 ...
        an.trial_lim_s(2)+(1/min(cfg_tfr.foi))*3];

    stim_events = trial_info.word_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));
   
    %% Compute TFRs        
    fprintf('------------- TFR Calculations ----------\n');
    tfr  = ft_freqanalysis(cfg_tfr, roi_trl);

    % Trim back down to original trial_lim_s to exclude NaNs
    cfg_trim = [];
    cfg_trim.latency = an.trial_lim_s;
    tfr = ft_selectdata(cfg_trim,tfr);
    
    %% Compute pairwise COH
    windowSize = 250; % in ms
    stepSize = 2; % in ms

    windowSamples = windowSize * (roi.fsample/1000);
    stepSamples = stepSize * (roi.fsample/1000);
    halfSamples = floor(windowSamples/2);
    
    t = round(tfr.time,3);

    newt = t(1+halfSamples:stepSamples:length(tfr.time)-halfSamples);
    frex = tfr.freq;

    COH = NaN(size(tfr.fourierspctrm,1),length(MPFC),length(LPFC),length(tfr.freq),length(newt));

    for chaniMPFC = 1:length(MPFC)
        for chaniLPFC = 1:length(LPFC)
            cnt = 1;
            for k = 1+halfSamples:stepSamples:length(tfr.time)-halfSamples
                
                % Cross-spectrum for current window
                Sml = tfr.fourierspctrm(:,strcmpi(tfr.label,MPFC{chaniMPFC}),:,k-halfSamples:k+halfSamples).*...
                    conj(tfr.fourierspctrm(:,strcmpi(tfr.label,LPFC{chaniLPFC}),:,k-halfSamples:k+halfSamples));
                
                % MPFC auto-sprectum
                Smm = tfr.fourierspctrm(:,strcmpi(tfr.label,MPFC{chaniMPFC}),:,k-halfSamples:k+halfSamples).*...
                    conj(tfr.fourierspctrm(:,strcmpi(tfr.label,MPFC{chaniMPFC}),:,k-halfSamples:k+halfSamples));

                % LPFC auto-spectrum
                Sll = tfr.fourierspctrm(:,strcmpi(tfr.label,LPFC{chaniLPFC}),:,k-halfSamples:k+halfSamples).*...
                    conj(tfr.fourierspctrm(:,strcmpi(tfr.label,LPFC{chaniLPFC}),:,k-halfSamples:k+halfSamples));

                COH(:,chaniMPFC,chaniLPFC,:,cnt) = abs(mean(Sml,4))./sqrt(mean(Smm,4).*mean(Sll,4));
                    
                cnt = cnt + 1;
            end
        end
    end

    
    %% Baseline Correction
    fprintf('===================================================\n');
    fprintf('------------- Baseline Correction ----------\n');
    fprintf('===================================================\n');

    zCOH = zbaseline_ISPC(COH,newt, -0.5, -0.2, 1000);

    clear COH
    
    %% Create Deisgn Matrix and filter trials
    sbj_list = {'CP24',
                'IR21',
                'IR26',
                'IR27',
                'IR31',
                'IR32',
                'IR35',
                'IR37',
                'IR39',
                'IR41',
                'IR48',
                'IR52',
                'IR54',
                'IR57',
                'IR61',
                'IR65',
                'IR67',
                'IR68',
                'IR72',
                'IR74',
                'IR79',
                'IR82',
                'IR84',
                'IR85'};

    [coherenceDesign,~] = fn_create_LMM_design(trial_info, length(LPFC), length(MPFC),'coherence');
    
    coherenceDesign.Subject = find(strcmpi(sbj_list,SBJ)).*ones(height(coherenceDesign), 1);
    

    %% Realign to response if desired
    if strcmpi(an.evnt_lab,'R')
        
        zCOH_r = NaN(size(zCOH,1),size(zCOH,2),size(zCOH,3),size(zCOH,4),length(-1:1/roi_trl.fsample:1));
        for triali = 1:size(zCOH,1)
            rt_idx = nearest(newt,trial_info.response_time(triali));
            zCOH_r(triali,:,:,:,:) = zCOH(triali,:,:,:,rt_idx-500:rt_idx+500);
        end
        newt = -1:1/roi_trl.fsample:1;
        zCOH = zCOH_r;
        clear zCOH_r;

        times2save = -1:.025:1; % in seconds

        tidx = dsearchn(newt',times2save');
    
        newt = newt(tidx);
        zCOH = zCOH(:,:,:,:,tidx);
        % Save the results for each subject
        save_dir = strcat(root_dir, ['PRJ_Stroop/results/COH/Resp/',SBJ,'/']);
    else
        times2save = -.2:.025:2; % in seconds
    
        tidx = dsearchn(newt',times2save');
    
        zCOH = zCOH(:,:,:,:,tidx);
        
        newt = newt(tidx);
    
        % Save the results for each subject
        save_dir = strcat(root_dir, ['PRJ_Stroop/results/COH/Stim/',SBJ,'/']);
    end
    
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    COH_save_filename = fullfile(save_dir, [SBJ '.mat']);

    save(COH_save_filename,'coherenceDesign','frex','newt', 'zCOH','MPFC','LPFC','elec','-v7.3');

end