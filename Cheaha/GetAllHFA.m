function GetAllHFA(SBJ, proc_id, an_id,alignment)

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

    %% Remove bad channels from data
    cfg = [];
    cfg.channel = SBJ_vars.ch_lab.ROI;
    data = ft_selectdata(cfg,data);
    
    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);

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
    hfa = extractHFA(data,-1.5,1.5,1000,trial_info,alignment); % extract 70-150 Hz HFA and z-score using boostrapping
    
    hfa.ztrial = hfa.ztrial(:,:,1:5:end);
   
    hfa.newtime = hfa.time(1:5:end);

    hfa.label = data.label;

    if strcmpi(alignment,'S')
        % Save the results for each subject
        save_dir = strcat(root_dir, 'PRJ_Stroop/results/Power/HFA/Stim/');
    
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        power_save_filename = fullfile(save_dir, [SBJ '.mat']);
    else
         % Save the results for each subject
        save_dir = strcat(root_dir, 'PRJ_Stroop/results/Power/HFA/Resp/');
    
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        
        power_save_filename = fullfile(save_dir, [SBJ '.mat']);
    end
    
    save(power_save_filename,'hfa','-v7.3');

end