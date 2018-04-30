%% EEG04_compute_ERPs
SBJ = 'IR57';
pipeline_id = 'eeg_ft';
ecog_pipeline_id = 'main_ft';

%% Add paths
addpath('/Users/ranaeser/knight/PRJ_Stroop/scripts/');
addpath('/Users/ranaeser/knight/PRJ_Stroop/scripts/utils/');
addpath('/Users/ranaeser/knight/Apps/fieldtrip/');
ft_defaults

%% load the data
eval(['run /Users/ranaeser/knight/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /Users/ranaeser/knight/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

eeg_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_cz_lap.mat');
eog_filename = strcat(SBJ_vars.dirs.preproc,SBJ,'_eog_filt.mat');
load(eeg_filename);
load(eog_filename);

load([SBJ_vars.dirs.preproc SBJ '_trial_info_eeg.mat'])

%% Cut data into trials
if strcmp(proc_vars.event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(proc_vars.event_type,'resp')
    events = trial_info.resp_onset;
else
    error(strcat('ERROR: unknown event_type ',proc_vars.event_type));
end

cz_trials = fn_ft_cut_trials_equal_len(cz_lap,events,...
    trial_info.condition_n',proc_vars.trial_lim_s*cz_lap.fsample);
eog_trials = fn_ft_cut_trials_equal_len(eog,events,...
    trial_info.condition_n',proc_vars.trial_lim_s*eog.fsample);

%% use ft_timelockanalysis to compute the ERPs 

cfg = [];
cz_ERP = ft_timelockanalysis(cfg, cz_trials);
ft_singleplotER(cfg, cz_ERP)

%ERPs based on trial type
con_trials = [];
inc_trials = [];
neu_trials = [];
for ix = 1:numel(trial_info.trialtype);
    if istrue(trial_info.trialtype{ix}=='con');
        con_trials = [con_trials, ix];
    elseif istrue(trial_info.trialtype{ix}=='inc');
        inc_trials = [inc_trials, ix];
    elseif istrue(trial_info.trialtype{ix}=='neu');
        neu_trials = [neu_trials, ix];
    end
end

%congruent trials
cfg = [];
cfg.trials = con_trials;
con = ft_timelockanalysis(cfg, cz_trials);
ft_singleplotER(cfg, con)

%incongruent trials
cfg = [];
cfg.trials = inc_trials;
inc = ft_timelockanalysis(cfg, cz_trials);
ft_singleplotER(cfg, inc)

%neutral trials 
cfg = [];
cfg.trials = neu_trials;
neu = ft_timelockanalysis(cfg, cz_trials);
ft_singleplotER(cfg, neu)


%% The following code allows you to look at the ERP difference waves
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
difference = ft_math(cfg, con, inc);
ft_singleplotER(cfg, difference)