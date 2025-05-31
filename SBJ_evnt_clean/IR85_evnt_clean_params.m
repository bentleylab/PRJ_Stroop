%% Photodiode Trace Cleaning Parameters: IR85
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = 1100;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 10.0],...% initial offset
    [761.1 761.8],...% blip between B6 T21-22
    [804.0 1150.0],...% long stretch where one blip is marked
    [1360.0 1370.0],...% long stretch where one blip is marked
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

