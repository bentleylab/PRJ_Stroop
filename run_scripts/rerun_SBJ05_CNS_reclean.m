%% rerun SBJ05 on all patients
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};
pipeline_id = 'main_ft';

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Step 0 - Processing Variables
eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

%%
for s = 11:numel(SBJs)
    SBJ = SBJs{s};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    %% check for tiny bad
    load([SBJ_vars.dirs.events SBJ '_colin_bad_epochs_preproc.mat']);
    tiny_bad = find(diff(bad_epochs,1,2)<10);
    if ~isempty(tiny_bad)
        warning([SBJ ' = Tiny bad epochs detected:\n']);
        disp(bad_epochs(tiny_bad,:));
        error(SBJ);
    end
    
    %% Load manually corrected trial_info
    ti = {};
    block_lens   = zeros([numel(SBJ_vars.block_name) 1]);
    block_times  = zeros([numel(SBJ_vars.block_name) 1]);
    block_trlcnt = zeros([numel(SBJ_vars.block_name) 1]);
    block_blkcnt = zeros([numel(SBJ_vars.block_name) 1]);
    for b_ix = 1:numel(SBJ_vars.block_name)
        if numel(SBJ_vars.raw_file)>1
            block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
            % Get block length
            tmp = load(strcat(SBJ_vars.dirs.import,SBJ,'_',...
                num2str(proc_vars.resample_freq),'hz',block_suffix,'.mat'));
            block_lens(b_ix) = size(tmp.data.trial{1},2);
            block_times(b_ix) = tmp.data.time{1}(end);
        else
            block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
        end
        ti{b_ix} = load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_manual',block_suffix,'.mat'));
        block_trlcnt(b_ix) = numel(ti{b_ix}.trial_info.trial_n);
        block_blkcnt(b_ix) = max(ti{b_ix}.trial_info.block_n);
        
        % Add run number
        ti{b_ix}.trial_info.run_n = ones(size(ti{b_ix}.trial_info.block_n))*b_ix;
    end
    
    % Combine trial_info structs if necessary
    if numel(ti)<2
        trial_info = ti{1}.trial_info;
    else
        trial_info = ti{1}.trial_info;
        % Add properties of the individual blocks
        trial_info.run_len        = block_lens;
        trial_info.run_time       = block_times;
        trial_info.run_trial_info = ti;     % Keep the original trial_info structs from each run
        for b_ix = 2:numel(SBJ_vars.block_name)
            % Concatenate fields that don't need modification
            %   Not modifying marker_time and onset_time (no idea what those are...)
            trial_info.word          = vertcat(trial_info.word,ti{b_ix}.trial_info.word);
            trial_info.color         = vertcat(trial_info.color,ti{b_ix}.trial_info.color);
            trial_info.trialtype     = vertcat(trial_info.trialtype,ti{b_ix}.trial_info.trialtype);
            trial_info.blocktype     = vertcat(trial_info.blocktype,ti{b_ix}.trial_info.blocktype);
            trial_info.response_time = vertcat(trial_info.response_time,ti{b_ix}.trial_info.response_time);
            trial_info.marker_time   = vertcat(trial_info.marker_time,ti{b_ix}.trial_info.marker_time);
            trial_info.onset_time    = vertcat(trial_info.onset_time,ti{b_ix}.trial_info.onset_time);
            trial_info.condition_n   = vertcat(trial_info.condition_n,ti{b_ix}.trial_info.condition_n);
            trial_info.error         = vertcat(trial_info.error,ti{b_ix}.trial_info.error);
            trial_info.run_n         = vertcat(trial_info.run_n,ti{b_ix}.trial_info.run_n);
            
            % Modify then concatenate counts and indices
            trial_info.block_n = vertcat(trial_info.block_n,ti{b_ix}.trial_info.block_n+sum(block_blkcnt(1:b_ix-1)));
            trial_info.trial_n = vertcat(trial_info.trial_n,ti{b_ix}.trial_info.trial_n+sum(block_trlcnt(1:b_ix-1)));
            trial_info.ignore_trials = horzcat(trial_info.ignore_trials,...
                ti{b_ix}.trial_info.ignore_trials+sum(block_trlcnt(1:b_ix-1)));
            
            trial_info.word_onset = vertcat(trial_info.word_onset,...
                ti{b_ix}.trial_info.word_onset+sum(block_lens(1:b_ix-1)));
            trial_info.resp_onset = vertcat(trial_info.resp_onset,...
                ti{b_ix}.trial_info.resp_onset+sum(block_lens(1:b_ix-1)));
        end
    end
    clear ti block_lens block_times block_trlcnt block_blkcnt
    
    % Toss trials based on behavior and cleaning with Bob
    SBJ05_fname = [SBJ_vars.dirs.events SBJ '_behavior_rejection_results.txt'];
    if exist(SBJ05_fname)
        system(['mv ' SBJ05_fname ' ' SBJ05_fname(1:end-4) '_preCNS.txt']);
    end
    trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,pipeline_id);
    
    %%    
    clear trial_info
    trial_info = trial_info_clean;
    % Document bad trials
    trial_info.bad_trials.variance = SBJ_vars.trial_reject_n';
    trial_info.bad_trials.all = sort([trial_info.bad_trials.all; trial_info.bad_trials.variance]);
    
    % Remove bad trials
    trial_rejected = ismember(trial_info.trial_n,SBJ_vars.trial_reject_n);
    trial_reject_ix = find(trial_rejected);
    trial_info.block_n(trial_reject_ix) = [];
    trial_info.trial_n(trial_reject_ix) = [];
    trial_info.word(trial_reject_ix) = [];
    trial_info.color(trial_reject_ix) = [];
    trial_info.trialtype(trial_reject_ix) = [];
    trial_info.blocktype(trial_reject_ix) = [];
    trial_info.response_time(trial_reject_ix) = [];
    trial_info.marker_time(trial_reject_ix) = [];
    trial_info.onset_time(trial_reject_ix) = [];
    trial_info.word_onset(trial_reject_ix) = [];
    trial_info.resp_onset(trial_reject_ix) = [];
    trial_info.condition_n(trial_reject_ix) = [];
    trial_info.error(trial_reject_ix) = [];
    
    % Make sure no response times are in weird float format
    trial_info.resp_onset = round(trial_info.resp_onset);
    
    final_ti_fname = strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat');
    if exist(final_ti_fname)
        system(['mv ' final_ti_fname ' ' final_ti_fname(1:end-4) '_preCNS.mat']);
    end
    save(final_ti_fname,'trial_info');
    
    clear SBJ SBJ_vars trial_info trial_info_clean
end