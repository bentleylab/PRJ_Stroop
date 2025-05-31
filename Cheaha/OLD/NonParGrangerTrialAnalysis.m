function NonParGrangerTrialAnalysis(SBJ, proc_id, an_id,cond)

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
   
    load(['/data/user/anaskhan/PRJ_Stroop/results/ThetaISPC/' SBJ '.mat'],'maxChan')

    cfgs.channel = maxChan';
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
    trial_lim_s_pad = [0 0.25];

    stim_events = trial_info.word_onset;

    roi_trl = fn_ft_cut_trials_equal_len(roi,stim_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) trial_lim_s_pad(2)]*roi.fsample));
   
    comparisons = {'conflict','Gratton'};
    for compi = 1:length(comparisons)
        if strcmpi(comparisons{compi},'conflict')

            % Determine whether congruent or incongruent is lower in
            % quantity
            nCon = sum(strcmpi(trial_info.trialtype,'con'));
            nInc = sum(strcmpi(trial_info.trialtype,'inc'));

            conTrials = find(strcmpi(trial_info.trialtype,'con'));
            incTrials = find(strcmpi(trial_info.trialtype,'inc'));

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

            % Do FFT      
            cfg_spec            = [];
            cfg_spec.method     = 'mtmfft';
            cfg_spec.taper      = 'hanning';
            cfg_spec.output     = 'fourier';
            cfg_spec.pad        = 1;
            freq_con            = ft_freqanalysis(cfg_spec, roi_trl_con);
            freq_inc            = ft_freqanalysis(cfg_spec, roi_trl_inc);

            % Get Granger Causality Spectrum
            cfg                     = [];
            cfg.method              = 'granger';
            cfg.granger.sfmethod    = 'bivariate';
            cfg.granger.conditional = 'no';
            granger_con             = ft_connectivityanalysis(cfg, freq_con);
            granger_inc             = ft_connectivityanalysis(cfg, freq_inc);

        elseif strcmpi(comparisons{compi},'Gratton')

            % Determine whether iC or cC is lower in quantity

            trials = fn_create_LMM_design(trial_info,1,1,'power');

            niC = sum(strcmpi(trials.CurrentType,'con') & strcmpi(trials.PreviousType,'inc'));
            ncC = sum(strcmpi(trials.CurrentType,'con') & strcmpi(trials.PreviousType,'con'));


            iCTrials = find(strcmpi(trials.CurrentType,'con') & strcmpi(trials.PreviousType,'inc'));
            cCTrials = find(strcmpi(trials.CurrentType,'con') & strcmpi(trials.PreviousType,'con'));
            
            cCIsSmaller = ncC < niC;

            if cCIsSmaller

                % Get subsample of trials
                subiC = randsample(iCTrials,ncC);

                % iC
                cfg_select = [];
                cfg_select.trials = subiC;
                roi_trl_iC = ft_selectdata(cfg_select,roi_trl);

                % cC
                cfg_select = [];
                cfg_select.trials = cCTrials;
                roi_trl_cC = ft_selectdata(cfg_select,roi_trl);

            elseif ~cCIsSmaller
                
                % Get subsample of trials
                subcC = randsample(cCTrials,niC);

                % cC
                cfg_select = [];
                cfg_select.trials = subcC;
                roi_trl_cC = ft_selectdata(cfg_select,roi_trl);

                % iC
                cfg_select = [];
                cfg_select.trials = iCTrials;
                roi_trl_iC = ft_selectdata(cfg_select,roi_trl);

            else

                % cC
                cfg_select = [];
                cfg_select.trials = cCTrials;
                roi_trl_cC = ft_selectdata(cfg_select,roi_trl);

                % iC
                cfg_select = [];
                cfg_select.trials = iCTrials;
                roi_trl_iC = ft_selectdata(cfg_select,roi_trl);

            end

            % Do FFT      
            cfg_spec            = [];
            cfg_spec.method     = 'mtmfft';
            cfg_spec.taper      = 'hanning';
            cfg_spec.output     = 'fourier';
            cfg_spec.pad        = 1;
            freq_cC             = ft_freqanalysis(cfg_spec, roi_trl_cC);
            freq_iC             = ft_freqanalysis(cfg_spec, roi_trl_iC);

            % Get Granger Causality Spectrum
            cfg                     = [];
            cfg.method              = 'granger';
            cfg.granger.sfmethod    = 'bivariate';
            cfg.granger.conditional = 'no';
            granger_cC              = ft_connectivityanalysis(cfg, freq_cC);
            granger_iC              = ft_connectivityanalysis(cfg, freq_iC);

        end
    end

    save_dir = strcat(root_dir, ['PRJ_Stroop/results/NPG/Stim/',SBJ,'/']);

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    NPG_save_filename = fullfile(save_dir, [SBJ '.mat']);
    save(NPG_save_filename,'granger_cC','granger_iC','granger_con', 'granger_inc','-v7.3');

end