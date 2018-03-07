%% Photodiode Trace Cleaning Parameters: IR39
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = 12280;

% Record epochs (in sec) with fluctuations that should be set to baseline
%	drop before B1, blip at end of B1
bsln_times = {...
    [0.0 12.0],...          % initial offset
    [117.0 118.0],...       % fake trial blip at end of B1
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {...
    };
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {...
    };
stim_yval = []; % inc and neu at those points in time
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

