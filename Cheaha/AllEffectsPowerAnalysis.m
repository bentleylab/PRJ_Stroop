function AllEffectsPowerAnalysis(SBJ, proc_id, an_id,cond)

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
    
    MPFC = elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC'));
    LPFC = elec.label(ismember(elec.ROI,'dlPFC'));

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
    tfr.powspctrm = log10(tfr.powspctrm);
    

    %% Baseline Correction
    fprintf('===================================================\n');
    fprintf('------------- Baseline Correction ----------\n');
    fprintf('===================================================\n');

    switch an.bsln_type
        case {'zscore', 'demean', 'my_relchange','zboot'}
            tfr = fn_bsln_ft_tfr(tfr,an.bsln_lim,an.bsln_type,an.bsln_boots);
        case 'relchange'
            cfgbsln = [];
            cfgbsln.baseline     = an.bsln_lim;
            cfgbsln.baselinetype = an.bsln_type;
            cfgbsln.parameter    = 'powspctrm';
            tfr = ft_freqbaseline(cfgbsln,tfr);
        otherwise
            error(['No baseline implemented for an.bsln_type: ' an.bsln_type]);
    end

    %% Create Deisgn Matrix and filter trials
    sbj_list = {'CP24',
                'IR21',
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

    [lpfcDesign,mpfcDesign] = fn_create_LMM_design(trial_info, length(LPFC), length(MPFC),'power');
    
    lpfcDesign.Subject = find(strcmpi(sbj_list,SBJ)).*ones(height(lpfcDesign), 1);
    mpfcDesign.Subject = find(strcmpi(sbj_list,SBJ)).*ones(height(mpfcDesign), 1);

    frex = tfr.freq;
    %% Realign to response if desired
    if strcmpi(an.evnt_lab,'R')
        tfr_r = tfr;
        tfr_r.powspctrm = NaN(size(tfr.powspctrm,1),size(tfr.powspctrm,2),size(tfr.powspctrm,3),length(-1:1/roi_trl.fsample:1));
        for triali = 1:size(tfr.powspctrm,1)
            rt_idx = nearest(tfr.time,trial_info.response_time(triali));
            tfr_r.powspctrm(triali,:,:,:) = tfr.powspctrm(triali,:,:,rt_idx-500:rt_idx+500);
        end
        tfr_r.time = -1:1/roi_trl.fsample:1;
        tfr = tfr_r;
        clear tfr_r;

        t = round(tfr.time,3);

        times2save = -0.75:.025:0.75; % in seconds
        times2plot = -0.75:.01:0.75; % in seconds

        tidx = dsearchn(t',times2save');
        plotidx = dsearchn(t',times2plot');
    
        newt = t(tidx);
        tplot = t(plotidx);
    
        tfr2plot = tfr; 
        tfr2plot.powspctrm = tfr.powspctrm(:,:,:,plotidx);

        tfr.powspctrm = tfr.powspctrm(:,:,:,tidx);
    else % if stim-locked

        t = round(tfr.time,3);

        times2save = -.1:.025:1.25; % in seconds
        times2plot = -.1:.01:1.25; % in seconds

        tidx = dsearchn(t',times2save');
        plotidx = dsearchn(t',times2plot');
    
        newt = t(tidx);
        tplot = t(plotidx);
    
        tfr2plot = tfr; 
        tfr2plot.powspctrm = tfr.powspctrm(:,:,:,plotidx);
        tfr2plot.time = tplot;
        
        tfr.powspctrm = tfr.powspctrm(:,:,:,tidx);
    end
   
    if strcmpi(an.evnt_lab,'S')
        % Save the results for each subject
        save_dir = strcat(root_dir, 'PRJ_Stroop/results/Power/Stim/');
    
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        power_save_filename = fullfile(save_dir, [SBJ '.mat']);
    else
         % Save the results for each subject
        save_dir = strcat(root_dir, 'PRJ_Stroop/results/Power/Resp/');
    
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        power_save_filename = fullfile(save_dir, [SBJ '.mat']);
    end

    save(power_save_filename,'frex','newt','tplot', 'tfr','tfr2plot','MPFC','LPFC','elec','lpfcDesign','mpfcDesign','-v7.3');

end