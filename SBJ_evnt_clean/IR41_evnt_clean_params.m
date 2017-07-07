%% Photodiode Trace Cleaning Parameters: IR41
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = 6500;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...% The photodiode data between Block 1 and Block 2 is bad
    [0.0 21.0],...% initial offset
    [122.0 140.0],...% extra trial at end of B1
    [240.0 280.0],...% fixes a dip in baseline b/w B2 and B3
    [606.0 616.0]...% fixes a dip in baseline b/w B5 and B6
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

