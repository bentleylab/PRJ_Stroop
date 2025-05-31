%% Photodiode Trace Cleaning Parameters: IR75
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% NOTE: all times and values below are after correcting drift and B8/9 scaling difference

% Set zero/baseline during a block
bsln_val = -231;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 156.0],...% initial offset
    [254.5 256.0],... % blip at end of B1
    [1130.0 1340.001]...% end of task
    };

% Record epochs (in sec) when baseline has shifted
%   this was taken care of in IR69_photo_minimal_processing.m
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
%   B7T36 (inc/medium) was when photodiode came loose (and caused scaling)
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

