function GetTrialBlockEncodingElecs(SBJ, proc_id)

    % Set up path stuff for Cheaha
    restoredefaultpath;
    root_dir = '/data/user/anaskhan/';
    ft_dir=[root_dir 'fieldtrip/'];

    % root_dir = '/Users/anaskhan/Desktop/';
    % ft_dir = [root_dir 'Apps/fieldtrip/'];
    % 
    addpath([root_dir 'PRJ_Stroop/scripts']);
    addpath([root_dir 'PRJ_Stroop/scripts/utils']);
    addpath(ft_dir);
    ft_defaults;
   
    %% Load Data
    
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);

    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    %% Trial labels
    trials = fn_create_LMM_design(trial_info,1,1,'power');
   
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
    %% Get block and previous trial encoding electrodes
    hfa = extractHFA(data,-0.5,-0.2,1000,trial_info,'S'); % extract 70-150 Hz HFA and z-score to -0.5:-0.2s using boostrapping
    
    % critical f-value for thresholding
    CritF  = finv(0.95,2,size(hfa.ztrial,1)); % critical f-value for block type (df1 = 2 (MC, EQ, MI = 3 - 1), df2 = # of trials)
    window_size = 0.05; % in seconds


    % hfa.ztrial(~strcmpi(trials.CurrentType,'con'),:,:) = [];
    % trials(~strcmpi(trials.CurrentType,'con'),:) = [];

    hfa.ztrial(strcmpi(trials.PreviousType,'None'),:,:) = [];
    trials(strcmpi(trials.PreviousType,'None'),:) = [];

    % time range 
    trange   = [hfa.time{1}(1) + window_size hfa.time{1}(end) - window_size];
    twin     = round(window_size * data.fsample); % define smoothing window in samples
    
    conditions = {'TrialType','BlockType','pTrialType'};
    
    for condi = 1:numel(conditions)
        switch conditions{condi}
            case 'BlockType'
                w2{condi}.time         = hfa.time{1}(nearest(hfa.time{1},trange(1)):nearest(hfa.time{1},trange(2)));
                w2{condi}.label        = data.label;
                w2{condi}.trial        = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize tiemcourse for exp. variance
                w2{condi}.pval         = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize pval timecourse
                w2{condi}.fstat        = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize f statistic timecourse
                w2{condi}.dimord       = 'chan_time';
                w2{condi}.sig_chans = false(numel(data.label),1);
                w2{condi}.cond = 'BlockType';
                % loop across channels
                for chani = 1:length(w2{condi}.label)
                    
                    fprintf('Channel %d/%d\n', chani, length(data.label));
                    
                    t = 1; % counter for timepoint
                
                    for tidx = nearest(hfa.time{1},trange(1)):nearest(hfa.time{1},trange(2))
            
                        % do the test with factor prediction
                        [p, table] = anovan(squeeze(mean(hfa.ztrial(:,chani,tidx-twin:tidx+twin),3)), ...
                            {trials.BlockType}, 'model', 'linear', 'sstype',2, ...
                            'varnames', {'BlockType'}, 'display', 'off');
            
                        % calc w2 [ω2 = (SSeffect – (dfeffect)(MSerror)) / MSerror + SStotal]
                        w2{condi}.trial(chani,t) = (table{2,2} - (table{2,3} * table{3,5})) / (table{3,5} + table{4,2});
            
                        % get significance
                        w2{condi}.pval(chani,t)  = p; 
                        w2{condi}.fstat(chani,t) = table{2,6};
                                  
                        % trial keeping
                        clear p table
                        t = t+1;
            
                    end %tidx
                    
                    % thresholding
                    fstat = w2{condi}.fstat(chani,:);
                    
                    fstat(fstat < CritF) = 0;
                   
                    % search for cluster for factor prediction
            
                    islands = bwconncomp(fstat);
                    for jj = 1:islands.NumObjects % loop through non-zero traces
                        
                        % check if non-zero elements occur after the verbal response
                        if all(w2{condi}.time(islands.PixelIdxList{jj}) > 0) && numel(islands.PixelIdxList{jj}) > 50 %ceil(sum(w2{condi}.time > 0) * min_pc)
                            w2{condi}.sig_chans(chani) = true;
                            break
                        end
                    end
                    
                    clear fstat islands
                               
                end % chani
            case 'pTrialType'
                w2{condi}.time         = hfa.time{1}(nearest(hfa.time{1},trange(1)):nearest(hfa.time{1},trange(2)));
                w2{condi}.label        = data.label;
                w2{condi}.trial        = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize tiemcourse for exp. variance
                w2{condi}.pval         = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize pval timecourse
                w2{condi}.fstat        = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize f statistic timecourse
                w2{condi}.dimord       = 'chan_time';
                w2{condi}.sig_chans = false(numel(data.label),1);
                w2{condi}.cond = 'PreviousTrial';
                % loop across channels
                for chani = 1:length(w2{condi}.label)
                    
                    fprintf('Channel %d/%d\n', chani, length(data.label));
                    
                    t = 1; % counter for timepoint
                
                    for tidx = nearest(hfa.time{1},trange(1)):nearest(hfa.time{1},trange(2))
            
                        % do the test with factor prediction
                        [p, table] = anovan(squeeze(mean(hfa.ztrial(:,chani,tidx-twin:tidx+twin),3)), ...
                            {trials.CurrentType}, 'model', 'linear', 'sstype',2, ...
                            'varnames', {'PreviousTrial'}, 'display', 'off');
            
                        % calc w2 [ω2 = (SSeffect – (dfeffect)(MSerror)) / MSerror + SStotal]
                        w2{condi}.trial(chani,t) = (table{2,2} - (table{2,3} * table{3,5})) / (table{3,5} + table{4,2});
            
                        % get significance
                        w2{condi}.pval(chani,t)  = p; 
                        w2{condi}.fstat(chani,t) = table{2,6};
                                  
                        % trial keeping
                        clear p table
                        t = t+1;
            
                    end %tidx
                    
                % thresholding
                fstat = w2{condi}.fstat(chani,:);
                
                fstat(fstat < CritF) = 0;
               
                % search for cluster for factor prediction
        
                islands = bwconncomp(fstat);
                
                for jj = 1:islands.NumObjects % loop through non-zero traces
                    
                    % check if non-zero elements occur after the verbal response
                    if all(w2{condi}.time(islands.PixelIdxList{jj}) > 0) && numel(islands.PixelIdxList{jj}) > 50 %ceil(sum(w2{condi}.time > 0) * min_pc)
                        w2{condi}.sig_chans(chani) = true;
                        break
                    end
                end %jj
                
                clear fstat islands
                end %chani

            case 'TrialType'
                w2{condi}.time         = hfa.time{1}(nearest(hfa.time{1},trange(1)):nearest(hfa.time{1},trange(2)));
                w2{condi}.label        = data.label;
                w2{condi}.trial        = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize tiemcourse for exp. variance
                w2{condi}.pval         = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize pval timecourse
                w2{condi}.fstat        = NaN(length(w2{condi}.label), length(w2{condi}.time)); % initialize f statistic timecourse
                w2{condi}.dimord       = 'chan_time';
                w2{condi}.sig_chans = false(numel(data.label),1);
                w2{condi}.cond = 'TrialType';
                % loop across channels
                for chani = 1:length(w2{condi}.label)
                    
                    fprintf('Channel %d/%d\n', chani, length(data.label));
                    
                    t = 1; % counter for timepoint
                
                    for tidx = nearest(hfa.time{1},trange(1)):nearest(hfa.time{1},trange(2))
            
                        % do the test with factor prediction
                        [p, table] = anovan(squeeze(mean(hfa.ztrial(:,chani,tidx-twin:tidx+twin),3)), ...
                            {trials.CurrentType}, 'model', 'linear', 'sstype',2, ...
                            'varnames', {'CurrentType'}, 'display', 'off');
            
                        % calc w2 [ω2 = (SSeffect – (dfeffect)(MSerror)) / MSerror + SStotal]
                        w2{condi}.trial(chani,t) = (table{2,2} - (table{2,3} * table{3,5})) / (table{3,5} + table{4,2});
            
                        % get significance
                        w2{condi}.pval(chani,t)  = p; 
                        w2{condi}.fstat(chani,t) = table{2,6};
                                  
                        % trial keeping
                        clear p table
                        t = t+1;
            
                    end %tidx
                    
                % thresholding
                fstat = w2{condi}.fstat(chani,:);
                
                fstat(fstat < CritF) = 0;
               
                % search for cluster for factor prediction
        
                islands = bwconncomp(fstat);
                
                for jj = 1:islands.NumObjects % loop through non-zero traces
                    
                    % check if non-zero elements occur after the verbal response
                    if all(w2{condi}.time(islands.PixelIdxList{jj}) > 0) && numel(islands.PixelIdxList{jj}) > 50 %ceil(sum(w2{condi}.time > 0) * min_pc)
                        w2{condi}.sig_chans(chani) = true;
                        break
                    end
                end %jj
                
                clear fstat islands
                end %chani
        end %switch/case
    end %condi
        
    % Save the results for each subject
    save_dir = strcat(root_dir, ['PRJ_Stroop/data/' SBJ '/04_proc/']);

    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    save_filename = fullfile(save_dir, [SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat']);

    save(save_filename,'w2','-v7.3');

end