function [combined] = fn_concat_trials(trials)
%% Concatenate two fieldtrip data structures with continuous data (e.g., two task blocks)

% Initialize matrices
n_trials = numel(trials.trial);
trial_lens = zeros([1 n_trials]);
end_times  = zeros([1 n_trials]);
for t_ix = 1:n_trials
    end_times(t_ix) = trials.time{t_ix}(end);
    trial_lens(t_ix) = numel(trials.time{t_ix});
end
time_step = 1/trials.fsample;

data_concat = zeros([numel(trials.label) sum(trial_lens)]);
time_concat = zeros([1 sum(trial_lens)]);

% Combine data
data_concat(:,1:trial_lens(1)) = trials.trial{1};
time_concat(1,1:trial_lens(1)) = trials.time{1};
for t_ix = 2:n_trials
    data_concat(:,trial_lens(t_ix-1)+1:trial_lens(t_ix-1)+trial_lens(t_ix)) = trials.trial{t_ix};
    time_concat(:,trial_lens(t_ix-1)+1:trial_lens(t_ix-1)+trial_lens(t_ix)) = trials.time{t_ix}+end_times(t_ix-1)+time_step;
end

% Create Fieldtrip struct
combined = trials;
combined.trial = data_concat;
combined.time = time_concat;

end