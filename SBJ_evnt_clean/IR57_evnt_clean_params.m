%% Photodiode Trace Cleaning Parameters: IR35
% Mark trials to ignore e.g., interruptions
ignore_trials = [1:5];  % see photod comment below

% Set zero/baseline during a block
bsln_val = 2710;

% Record epochs (in sec) with fluctuations that should be set to baseline
%	drop before B1, blip at end of B1
bsln_times = {...
    [0.0 22.75],...% first 4 trials are missing photodiode; onset of 5th is messed up (toss all 5)
    [107.0 109.0]...% extra blip of a trial at the end of B1
    };

% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {[22.96 24.5]};
stim_yval = [30005]; % inc and neu at those points in time
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

