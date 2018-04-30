SBJ = 'IR57';
pipeline_id = 'main_ft';
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

load(['~/PRJ_Stroop/data/' SBJ '/Rana_EEG/' SBJ '_cz_lap.mat']);
load(['~/PRJ_Stroop/data/' SBJ '/Rana_EEG/' SBJ '_eog_filt.mat']);
load(['~/PRJ_Stroop/data/' SBJ '/Rana_EEG/' SBJ '_trial_info_eeg.mat']);
all_ti = load(['~/PRJ_Stroop/data/' SBJ '/03_events/' SBJ '_trial_info_manual.mat']);

events = all_ti.trial_info.word_onset;
cz_trials = fn_ft_cut_trials_equal_len(cz_lap,events,all_ti.trial_info.condition_n',proc_vars.trial_lim_s*cz_lap.fsample);
eog_trials = fn_ft_cut_trials_equal_len(eog,events,all_ti.trial_info.condition_n',proc_vars.trial_lim_s*eog.fsample);
cfga = []; both_trials = ft_appenddata(cfga,cz_trials,eog_trials);

events = trial_info.word_onset;
cz_trials = fn_ft_cut_trials_equal_len(cz_lap,events,trial_info.condition_n',proc_vars.trial_lim_s*cz_lap.fsample);
eog_trials = fn_ft_cut_trials_equal_len(eog,events,trial_info.condition_n',proc_vars.trial_lim_s*eog.fsample);
cfga = []; both_trials = ft_appenddata(cfga,cz_trials,eog_trials);

% all_ti = load([SBJ '_trial_info_manual.mat']);
% final_t_ix = [];
% for t_ix = 1:numel(trial_info.bad_trials.cz_ft_sum)
%     final_t_ix = [final_t_ix find(all_ti.trial_info.trial_n==trial_info.bad_trials.cz_ft_sum(t_ix))];
% end


% trial butterfly
figure; hold on;
plot_data = cz_trials;
for t_ix = 1:numel(plot_data.trial)
    plot(plot_data.trial{t_ix});
end

% single trial Stack
figure; hold on;
plot_data = cz_trials;
stack = NaN(numel(plot_data.trial), size(plot_data.trial{1},2));
for t_ix = 1:numel(plot_data.trial)
    stack(t_ix,:) = plot_data.trial{t_ix};
end
imagesc(stack);

% reject visual
cfg_reject = [];
cfg_reject.method = 'channel';
tmp = ft_rejectvisual(cfg_reject,both_trials);