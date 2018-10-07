%% Photodiode Trace Cleaning Parameters: IR67
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% NOTE: all times and values below are after correcting drift and B8/9 scaling difference

% Set zero/baseline during a block
bsln_val = -218;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 13.0],...% initial offset (attaching to screen?)
    [112.0 114.0]...% extra blip of a trial at the end of B1
    };

% Record epochs (in sec) when baseline has shifted
%   drift throughout: cftools says -2.958e-05 slope, -297.2 intercept (ignored)
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

