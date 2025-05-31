function AllEffectsITPCAnalysis(SBJ, proc_id, an_id,cond)

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

    MPFC = elec.label(ismember(elec.ROI,'aMCC'));
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

    % Trim back down to original trial_lim_s
    cfg_trim = [];
    cfg_trim.latency = an.trial_lim_s;
    tfr = ft_selectdata(cfg_trim,tfr);
    tfr.angles = angle(tfr.fourierspctrm);
    tfr = rmfield(tfr,'fourierspctrm');
    
    %% Compute ITPC and jackknife procedure
    
    full_itpc = abs(mean(exp(1i*tfr.angles),1));

    % Number of trials
    N_trials = size(tfr.angles, 1);
    
    % Pre-allocate storage for LOO ITPC
    LOO_itpc = zeros(N_trials, size(tfr.angles, 2), size(tfr.angles, 3), size(tfr.angles, 4));
    
    for itrial = 1:N_trials
        % Create a logical index for trials to include
        trial_idx = true(N_trials, 1);
        trial_idx(itrial) = false;  % Exclude the current trial
        
        % Extract angles excluding the current trial
        LOO_angles = tfr.angles(trial_idx, :, :, :);
        
        % Compute LOO ITPC
        LOO_itpc_curr = abs(mean(exp(1i * LOO_angles), 1));
        
        % Store the LOO ITPC for the current trial
        LOO_itpc(itrial, :, :, :) = full_itpc - LOO_itpc_curr;
    end

    t = round(tfr.time,3);

   
    frex = tfr.freq;

    
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

    [lpfcDesign,mpfcDesign] = fn_create_LMM_design(trial_info, length(LPFC), length(MPFC),'power');
    
    lpfcDesign.Subject = find(strcmpi(sbj_list,SBJ)).*ones(height(lpfcDesign), 1);
    mpfcDesign.Subject = find(strcmpi(sbj_list,SBJ)).*ones(height(mpfcDesign), 1);    
%%  Downsample

    t = round(tfr.time,3);

    times2save = -.2:.025:2; % in seconds
    times2plot = -.2:.01:2; % in seconds

    tidx = dsearchn(t',times2save');
    plotidx = dsearchn(t',times2plot');

    newt = t(tidx);
    tplot = t(plotidx);

    tfr = rmfield(tfr,'angles');
    tfr.itpc = LOO_itpc;
    tfr2plot = tfr; 
    tfr2plot.itpc = tfr.itpc(:,:,:,plotidx);
    tfr2plot.time = tplot;
    
    tfr.itpc = tfr.itpc(:,:,:,tidx);
    tfr.time = newt;

    % Save the results for each subject
    save_dir = strcat(root_dir, ['PRJ_Stroop/results/ITPC/Stim/',SBJ,'/']);

    
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    ITPC_save_filename = fullfile(save_dir, [SBJ '.mat']);

    save(ITPC_save_filename,'lpfcDesign','mpfcDesign','frex','newt', 'MPFC','LPFC','elec','tfr','tfr2plot','-v7.3');

end