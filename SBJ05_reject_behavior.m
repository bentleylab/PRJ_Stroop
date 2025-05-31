function trial_info_clean = SBJ05_reject_behavior(SBJ,trial_info,proc_id)
% Select data of interest and reject all bad trials and channels to get
% final clean dataset for processing
% Criteria:
%   1. Bob bad epoch
%   2. Bad trial (error, interruption, corrupt data, etc.)
%   3. Bad RT (no response, RT outlier, etc.)
% Inputs:
%   SBJ [str]- the dataset to process (e.g., 'IR54')
%   trial_info [struct]- structure with trial information
%       should be SBJ_trial_info_manual.mat from events_dir as saved by SBJ03_RT_manual_adjustments.m
%   evnt_lab [str]- 'S' or 'R', determines the event to which trials are locked
%   trial_lim_s [2x1 float array]- boundaries in SECONDS of trials around events
%       trial_lim(1)- baseline, e.g., -0.5 would be 500 ms before event of interest
%       trial_lim(2)- post-event length, e.g., 2 would be 2000 ms after event of interest
%   RT_std_thresh [int]- # standarad deviation from mean for RT to be tossed as outlier
% Outputs:
%   trial_info_clean_behav [struct]- saves out final version after tossing all bad trials

% Directories
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
elseif exist('/Volumes/hoycw_clust/','dir');root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';
else root_dir='/Users/anaskhan/Desktop/';ft_dir=[root_dir 'Apps/fieldtrip/'];end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

%% Load data
eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' proc_id '_vars.m']);
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));

% Load and convert bad_epochs
colin = load(strcat(SBJ_vars.dirs.events,SBJ,'_colin_bad_epochs_preproc.mat'));
% bob = [];
% for b_ix = 1:numel(SBJ_vars.block_name)
%     if numel(SBJ_vars.raw_file)>1
%         block_suffix = strcat('_',SBJ_vars.block_name{b_ix});
%     else
%         block_suffix = SBJ_vars.block_name{b_ix};   % should just be ''
%     end
%     run_epochs = load(strcat(SBJ_vars.dirs.events,SBJ,'_bob_bad_epochs_preclean',block_suffix,'.mat'));
%     % Convert to analysis_time within it's run
%     %   NOTE: 1 keeps epochs that are only partially overlaping real data
%     %   (i.e., will be trimmed to edges of analysis_time)
%     if ~isempty(run_epochs.bad_epochs)
%         at_epochs  = fn_convert_epochs_full2at(run_epochs.bad_epochs,SBJ_vars.analysis_time{b_ix},...
%             strcat(SBJ_vars.dirs.preproc,SBJ,'_preclean',block_suffix,'.mat'),1);
%         % Convert to combined run time
%         if b_ix>1
%             run_epochs.bad_epochs = run_epochs.bad_epochs+sum(trial_info.run_len(1:b_ix-1));
%         end
%         bob = vertcat(bob,run_epochs.bad_epochs);
%     end
% end

% Parameters
if ~isfield(proc,'RT_std_thresh')
    proc.RT_std_thresh = 3;
end
if ~isfield(proc,'trial_lim_s')
    proc.trial_lim_s = [-0.500 2.5];
end

%% Select channels and events of interest
if strcmp(proc.evnt_lab,'S')
    events = trial_info.word_onset;
elseif strcmp(proc.evnt_lab,'R')
    events = trial_info.resp_onset;
else
    error(stract('ERROR: unknown evnt_lab ',proc.evnt_lab));
end

% Convert trial_lim into samples
trial_lim = proc.trial_lim_s*data.fsample;

%% Reject known artifacts
skip_rt1 = find(trial_info.resp_onset<0);     % couldn't determine RT
skip_rt2 = find(isnan(trial_info.resp_onset));  % no response
if ~isempty(skip_rt2)
    fprintf('WARNING! Found a NaN in trial_info.resp_onset... will be tossed, but shouldnt that be a -1?\n');
end
skip_err = find(trial_info.error==1);           % error
skip_bad = find(trial_info.error<0);          % bad trial

% % Toss epochs that overlap with bad_epochs from Bob
% if ~isempty(bob)
%     skip_bob = fn_find_trials_overlap_epochs(bob,1:size(data.trial{1},2),events,trial_lim);
% else
%     skip_bob = [];
% end
% Toss epochs that overlap with bad_epochs from Colin (extra pass over preprocessed data)
if ~isempty(colin.bad_epochs)
    skip_colin = fn_find_trials_overlap_epochs(colin.bad_epochs,1:size(data.trial{1},2),events,trial_lim);
else
    skip_colin = [];
end

% Find RT outliers
RT_mean = nanmean(trial_info.response_time);
RT_std  = nanstd(trial_info.response_time);
skip_rt_outlier = find(abs(trial_info.response_time-RT_mean)>proc.RT_std_thresh*RT_std);
skip_rt_outlier = setdiff(skip_rt_outlier,skip_rt1);    %don't include RT=-1

% Check against RT bounds, toss late but only warn for early (don't toss)
RT_late = find(trial_info.response_time>proc.rt_bounds(2));
% Toss the trial following the late response (contaminated by
% response/monitoring during stimulus processing)
RT_late_next = RT_late+1;
% Don't toss next trial if it's in a different block (or is the very last trial)
RT_late_next(logical(diff(horzcat(trial_info.block_n(RT_late),trial_info.block_n(RT_late_next)),1,2))) = [];
RT_late_next(RT_late_next>max(trial_info.trial_n)) = [];
skip_rt_outlier = vertcat(skip_rt_outlier, RT_late, RT_late_next);
RT_early = find(trial_info.response_time(trial_info.response_time>0)<proc.rt_bounds(1));

% Compile all a priori bad trials
skip_trial_ix = unique([skip_bad; skip_rt1; skip_rt2; skip_err; skip_colin; skip_rt_outlier]);%skip_bob; 
ok_trial_ix = setdiff(1:length(trial_info.resp_onset),skip_trial_ix);

%% Compile Bad Trials
trial_info_clean = trial_info;
trial_info_clean.evnt_lab = proc.evnt_lab;
trial_info_clean.trial_lim = trial_lim;
trial_info_clean.trial_lim_s = proc.trial_lim_s;
trial_info_clean.RT_std_thresh = proc.RT_std_thresh;

% Document bad trials
trial_info_clean.bad_trials.RT_bad = trial_info.trial_n([skip_rt1 skip_rt2]);
trial_info_clean.bad_trials.RT_out = trial_info.trial_n(skip_rt_outlier);
% trial_info_clean.bad_trials.bob    = trial_info.trial_n(skip_bob);
trial_info_clean.bad_trials.colin  = trial_info.trial_n(skip_colin);
trial_info_clean.bad_trials.error  = trial_info.trial_n(skip_err);
trial_info_clean.bad_trials.bad    = trial_info.trial_n(skip_bad);
trial_info_clean.bad_trials.all    = trial_info.trial_n(skip_trial_ix);

% Remove bad trials
trial_info_clean.block_n(skip_trial_ix)       = [];
trial_info_clean.trial_n(skip_trial_ix)       = [];
trial_info_clean.word(skip_trial_ix)          = [];
trial_info_clean.color(skip_trial_ix)         = [];
trial_info_clean.trialtype(skip_trial_ix)     = [];
trial_info_clean.blocktype(skip_trial_ix)     = [];
trial_info_clean.response_time(skip_trial_ix) = [];
trial_info_clean.marker_time(skip_trial_ix)   = [];
trial_info_clean.onset_time(skip_trial_ix)    = [];
trial_info_clean.word_onset(skip_trial_ix)    = [];
trial_info_clean.resp_onset(skip_trial_ix)    = [];
trial_info_clean.condition_n(skip_trial_ix)   = [];
trial_info_clean.error(skip_trial_ix)         = [];

%prevent formatting errors
trial_info_clean.resp_onset = round(trial_info_clean.resp_onset);

% Print results
fprintf('==============================================================================================\n');
if ~isempty(RT_late)
    fprintf('WARNING! %i RTs > %f sec excluded, plus %i trailing trials!\n',...
        numel(RT_late),proc.rt_bounds(2),numel(RT_late_next));
end
if ~isempty(RT_early)
    fprintf('WARNING! %i RTs < %f sec detected (not excluded)!\n',numel(RT_early),proc.rt_bounds(1));
end
fprintf('Num trials excluded for bad RT     : %i\n',length(skip_rt1)+length(skip_rt2));
fprintf('Num trials excluded for outlier RT : %i\n',length(skip_rt_outlier));
fprintf('Num trials excluded for errors     : %i\n',length(skip_err));
% fprintf('Num trials excluded by PRECLEAN vis: %i\n',length(skip_bob));
fprintf('Num trials excluded by PREPROC vis : %i\n',length(skip_colin));
fprintf('Num trials excluded for other      : %i\n',length(skip_bad));
fprintf('TOTAL TRIALS EXCLUDED A PRIORI     : %i\n',length(skip_trial_ix));
fprintf('TRIALS REMAINING: %i/%i\n',length(trial_info_clean.trial_n),length(trial_info.response_time));
fprintf('==============================================================================================\n');

% Save results
results_filename = [SBJ_vars.dirs.events SBJ '_behavior_rejection_results.txt'];
r_file = fopen(results_filename,'a');
fprintf(r_file,'%s\n',datestr(datetime));
if ~isempty(RT_late)
    fprintf(r_file,'WARNING! %i RTs > %f sec excluded!\n',numel(RT_late),proc.rt_bounds(2));
end
if ~isempty(RT_early)
    fprintf(r_file,'WARNING! %i RTs < %f sec detected (not excluded)!\n',numel(RT_early),proc.rt_bounds(1));
end
fprintf(r_file,'Num trials excluded for bad RT    : %i\n',length(skip_rt1)+length(skip_rt2));
fprintf(r_file,'Num trials excluded for outlier RT: %i\n',length(skip_rt_outlier));
fprintf(r_file,'Num trials excluded for errors    : %i\n',length(skip_err));
% fprintf(r_file,'Num trials excluded by Bob vis    : %i\n',length(skip_bob));
fprintf(r_file,'Num trials excluded for other     : %i\n',length(skip_bad));
fprintf(r_file,'TOTAL TRIALS EXCLUDED A PRIORI    : %i\n',length(skip_trial_ix));
fprintf(r_file,'TRIALS REMAINING: %i/%i\n',length(trial_info_clean.trial_n),length(trial_info.response_time));
fclose(r_file);

end
