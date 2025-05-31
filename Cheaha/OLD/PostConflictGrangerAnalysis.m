function PostConflictGrangerAnalysis(SBJ, proc_id, an_id,cond)

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

    % MPFC = elec.label(ismember(elec.gROI,'MPFC'));
    % LPFC = elec.label(ismember(elec.gROI,'LPFC'));

    MPFC = elec.label(ismember(elec.ROI,'aMCC'));
    LPFC = elec.label(ismember(elec.ROI,'dlPFC'));

    if isempty(MPFC) || isempty(LPFC)
    error('No MPFC or LPFC channels found.');
    end

    % load(strcat("/data/user/anaskhan/PRJ_Stroop/results/Power/",SBJ,"/",SBJ,".mat"),'m_theta_elec',"l_theta_elec");

    % cfgs.channel = [m_theta_elec;l_theta_elec];

    cfgs.channel = [MPFC;LPFC];
    roi = ft_selectdata(cfgs,data);
    origsamp = roi.fsample;
    
           
    %% Get task-active electrodes
    % hfa = extractHFA(roi,-0.5,-0.2,1000,trial_info); % extract 70-150 Hz HFA and z-score to -0.5:-0.2s using boostrapping
    % 
    % % Get only congruent and incongruent trials (ignore neutrals)
    % CI_trials = logical(sum(cond_idx,1));
    % 
    % % Start counting at stimulus onset
    % start_trial_idx = dsearchn(hfa.time{1}',0);
    % 
    % % Get mean over all trials (con + inc) for each channel
    % hfa.meanHFA = squeeze(mean(hfa.ztrial(CI_trials,:,start_trial_idx:end),1));
    % 
    % % Two-way test against baseline
    % meanHFA_thresh = abs(hfa.meanHFA) >= 1.96;
    % 
    % % Initialize channels matrix
    % actv_chans = false(size(meanHFA_thresh,1),1);
    % 
    % for chani = 1:size(meanHFA_thresh,1)
    %     islands = bwconncomp(meanHFA_thresh(chani,:));
    %     if islands.NumObjects == 0
    %         continue
    %     else
    %         if any(cellfun(@numel,islands.PixelIdxList) >= 0.1*roi.fsample) % > 100 ms
    %             actv_chans(chani,1) = true;
    %         end
    %     end        
    % end
    % 
    %  % Logical indices for MPFC and LPFC channels within roi.label
    % mpfc_indices = ismember(roi.label, MPFC);
    % lpfc_indices = ismember(roi.label, LPFC);
    % 
    % % Check if there's at least one active channel in MPFC and LPFC
    % mpfc_active = any(actv_chans(mpfc_indices));
    % lpfc_active = any(actv_chans(lpfc_indices));
    % 
    % if ~(mpfc_active && lpfc_active)
    %     error(['Subject ' SBJ ' does not contain active channels in both LPFC and MPFC']);
    % end
    % 
    % % Select for active MPFC and LPFC channels
    % aMPFC = roi.label(mpfc_indices & actv_chans);
    % aLPFC = roi.label(lpfc_indices & actv_chans);

    % cfgs.channel = [aMPFC;aLPFC];
    % roi = ft_selectdata(cfgs,roi);

     %% Cut into Trials
    % trial_lim_s_pad = [-0.75 2.25];
    trial_lim_s_pad = [0 1];
    stim_events = trial_info.resp_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));
   
    %% Bring all subjects' fs down to 250 Hz

    cfg_resamp = [];
    cfg_resamp.resamplefs = 250;
    cfg_resamp.method = 'downsample';
    cfg_resamp.lpfilter = 'no';

    roi_trl = ft_resampledata(cfg_resamp,roi_trl);
    trial_info.sample_rate = roi_trl.fsample;

    trial_info.word_onset = round(trial_info.word_onset./(origsamp/roi.fsample));
    trial_info.resp_onset = round(trial_info.resp_onset./(origsamp/roi.fsample));

    if isfield(trial_info,'run_len')
        trial_info.run_len = round(trial_info.run_len./(origsamp/roi.fsample));
    end
    %% Create Deisgn Matrix and filter trials
  
    % aMPFC = MPFC;
    % aLPFC = LPFC;

    % aMPFC = m_theta_elec;
    % aLPFC = l_theta_elec;

    [coherenceDesign,~] = fn_create_LMM_design(trial_info, 1, 1,'power');
    coherenceDesign(coherenceDesign.Electrode ~= 1,:) = [];

    allData = cat(3, roi_trl.trial{:});
   
    errors = logical(trial_info.error);
    allData = allData(:,:,~errors);

    % Reorder channels in X so that all MPFC come first followed by LPFC
    
    new_order_channels = [MPFC; LPFC];

    new_order_indices = zeros(length(new_order_channels), 1);

    % Loop over the new order of channel names to find their current indices
    for i = 1:length(new_order_channels)
        % Find the index of each channel name in the current ordering
        idx = find(strcmp(roi_trl.label, new_order_channels{i}));
        new_order_indices(i) = idx;
    end

    allData = allData(new_order_indices,:,:);

    % X_cI = allData(:,:,strcmpi('inc',coherenceDesign.CurrentType) & strcmpi('con',coherenceDesign.PreviousType));
    % X_iI = allData(:,:,strcmpi('inc',coherenceDesign.CurrentType) & strcmpi('inc',coherenceDesign.PreviousType));
    allData = (allData - mean(allData,2))./std(allData,[],2);
    allData = (allData - mean(allData,3))./std(allData,[],3);

    %% Compute Granger Causality       
    fprintf('------------- Granger Calculations ----------\n');
    
    n_MPFC_channels = length(MPFC);
    n_LPFC_channels = length(LPFC);

    % window_size = 0.5; % in sec
    % halfWindow_size = window_size/2;
    % half_samples = round(halfWindow_size*roi_trl.fsample);
    % 
    % step_size = 0.050; % in sec
    % step_samples = round(step_size*roi_trl.fsample);
    % 
    % grangerTime = roi_trl.time{1,1}(1+half_samples:step_samples:end-half_samples);

    % for mchani = 1:n_MPFC_channels
    %     for lchani = 1:n_LPFC_channels
    %         mpfc_index = mchani; % MPFC channel index
    %         lpfc_index = n_MPFC_channels + lchani; % LPFC channel index, adjusted for position in the combined array
    %         cnt = 1;
    %         for tidx = 1+half_samples:step_samples:length(roi_trl.time{1,1})-half_samples
    %             [~,~,~,moBIC(mchani,lchani,cnt)] = tsdata_to_infocrit(X_cI([mpfc_index lpfc_index],tidx-half_samples:tidx+half_samples,:),20,'LWR');
    %             cnt = cnt + 1;
    %         end
    %     end
    % end

    % for mchani = 1:n_MPFC_channels
    %     for lchani = 1:n_LPFC_channels
    %         skipIteration = false; % Initialize the flag variable
    %         mpfc_index = mchani; % MPFC channel index
    %         lpfc_index = n_MPFC_channels + lchani;
    %         cnt = 1;
    %         for tidx = 1+half_samples:step_samples:length(roi_trl.time{1,1})-half_samples
    %             [A_cI, SIG_cI] = tsdata_to_var(X_cI([mpfc_index lpfc_index], tidx-half_samples:tidx+half_samples,:), 12, 'LWR');
    %             [A_iI, SIG_iI] = tsdata_to_var(X_iI([mpfc_index lpfc_index], tidx-half_samples:tidx+half_samples,:), 12, 'LWR');
    %             if ~isbad(A_cI) && ~isbad(A_iI)
    %                 infocI = var_info(A_cI, SIG_cI);
    %                 assert(~infocI.error, 'VAR error(s) found - bailing out');
    % 
    %                 infoiI = var_info(A_iI, SIG_iI);
    %                 assert(~infoiI.error, 'VAR error(s) found - bailing out');
    % 
    %                 fcI = var_to_spwcgc(A_cI, SIG_cI, 2^nextpow2(length(X_cI)));
    %                 m_to_l_GC_cI(mchani, lchani, :, cnt) = fcI(2, 1, :) - fcI(1, 2, :);
    % 
    %                 fiI = var_to_spwcgc(A_iI, SIG_iI, 2^nextpow2(length(X_iI)));
    %                 m_to_l_GC_iI(mchani, lchani, :, cnt) = fiI(2, 1, :) - fiI(1, 2, :);
    %                 cnt = cnt + 1;
    %             else
    %                 skipIteration = true; % Set the flag to skip the rest of this lchani iteration
    %                 break; % Exit the tidx loop
    %             end
    %         end
    %         if skipIteration
    %             m_to_l_GC_cI(mchani, lchani, :, :) = NaN;
    %             m_to_l_GC_iI(mchani, lchani, :, :) = NaN;
    %             continue; % Skip to the next iteration of lchani
    %         end
    %     end
    % end

    % for mchani = 1:n_MPFC_channels
    %     for lchani = 1:n_LPFC_channels
    %         mpfc_index = mchani; % MPFC channel index
    %         lpfc_index = n_MPFC_channels + lchani;
    %         cnt = 1;
    %         for tidx = 1+half_samples:step_samples:length(roi_trl.time{1,1})-half_samples
    %             temp_cI = reshape(X_cI([mpfc_index lpfc_index],tidx-half_samples:tidx+half_samples,:),2,[]);
    %             temp_iI = reshape(X_iI([mpfc_index lpfc_index],tidx-half_samples:tidx+half_samples,:),2,[]);
    % 
    %             temp_cI = zscore(temp_cI,[],2);
    %             temp_iI = zscore(temp_iI,[],2);
    % 
    %             [A_cI, SIG_cI] = tsdata_to_var(temp_cI, 16);
    %             [A_iI, SIG_iI] = tsdata_to_var(temp_iI, 16);
    % 
    %             fcI = var_to_spwcgc(A_cI, SIG_cI, 125);
    %             l_to_m_GC_cI(mchani, lchani, :, cnt) = squeeze(fcI(1, 2, :));
    %             m_to_l_GC_cI(mchani, lchani, :, cnt) = squeeze(fcI(2, 1, :));
    % 
    %             fiI = var_to_spwcgc(A_iI, SIG_iI, 125);
    %             l_to_m_GC_iI(mchani, lchani, :, cnt) = squeeze(fiI(1, 2, :));
    %             m_to_l_GC_iI(mchani, lchani, :, cnt) = squeeze(fiI(2, 1, :));
    %             cnt = cnt + 1; 
    %         end
    % 
    %     end
    % end

    % for mchani = 1:n_MPFC_channels
    %     for lchani = 1:n_LPFC_channels
    %         mpfc_index = mchani; % MPFC channel index
    %         lpfc_index = n_MPFC_channels + lchani;
    %         cnt = 1;
    %         for triali = 1:size(allData,3)
    %             for tidx = 1+half_samples:step_samples:length(roi_trl.time{1,1})-half_samples
    %                 [A, SIG] = tsdata_to_var(allData([mpfc_index lpfc_index],tidx-half_samples:tidx+half_samples,triali), 16);
    % 
    %                 f = var_to_spwcgc(A, SIG, 125);
    %                 l_to_m_GC(mchani, lchani, :, cnt,triali) = squeeze(f(1, 2, :));
    %                 m_to_l_GC(mchani, lchani, :, cnt,triali) = squeeze(f(2, 1, :));
    %                 cnt = cnt + 1; 
    %             end
    % 
    %         end
    %     end
    % end
    
    for mchani = 1:n_MPFC_channels
        for lchani = 1:n_LPFC_channels
            mpfc_index = mchani; % MPFC channel index
            lpfc_index = n_MPFC_channels + lchani;
            for triali = 1:size(allData,3)
                [A, SIG] = tsdata_to_var(allData([mpfc_index lpfc_index],:,triali), 12);
                    
                f = var_to_spwcgc(A, SIG, 125);
                l_to_m_GC(mchani, lchani, :,triali) = squeeze(f(1, 2, :));
                m_to_l_GC(mchani, lchani, :,triali) = squeeze(f(2, 1, :));
            end
        end
    end
         
    frex = sfreqs(125,roi_trl.fsample);

    %% Save the results for each subject
    save_dir = strcat(root_dir, ['PRJ_Stroop/results/PGranger/',SBJ,'/']);

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    NPG_save_filename = fullfile(save_dir, [SBJ '.mat']);

    save(NPG_save_filename,'frex','l_to_m_GC','m_to_l_GC',"coherenceDesign",'-v7.3');

end