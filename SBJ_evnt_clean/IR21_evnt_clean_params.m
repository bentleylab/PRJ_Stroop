%% Photodiode Trace Cleaning Parameters: IR21
% Mark trials to ignore e.g., interruptions
%   12 blocks of 24 trials in this patient
%   disruption at start of B2 due to screen freeze
ignore_trials = [25 26];

% Set zero/baseline during a block
bsln_val = 18000;

% Record epochs (in sec) with fluctuations that should be set to baseline
%   before B1, b/w B1/B2, first two trials of B2, very end
%   12 blocks, The photodiode data between Block 1 and Block 2 is bad
%   the signal stays high between B1-B2, then squiggles before kicking in on the first trial
bsln_times = {...% The photodiode data between Block 1 and Block 2 is bad
    [0.0 11.0],...% initial offset
    [102.0 160.0],...% extra trial at end of B1
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

