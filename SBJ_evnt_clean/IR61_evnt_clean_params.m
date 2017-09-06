%% Photodiode Trace Cleaning Parameters: IR61
% Mark trials to ignore e.g., interruptions
ignore_trials = [];  % see photod comment below

% Set zero/baseline during a block
bsln_val = 3324;

% Record epochs (in sec) with fluctuations that should be set to baseline
%	drop before B1, blip at end of B1
%   some noisy stuff in between most blocks, so zeroing that out
bsln_times = {...
    [0.0 10.0],...
    [108.0 152.0],...% extra blip at the end of B1 to B2
    [250.0 283.0],...%B2 to B3
    [382.0 415.0],...%B3 to B4
    [513.0 626.0],...%B4 to B5
    [724.0 745.0],...%B5 to B6
    [844.0 884.0],...%B6 to B7
    [985.0 995.0],...%B7 to B8
    [1100.0 1112.0],...%B8 to B9
    [1213.0 1224.0]...%B9 to end
    };

% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {...
    };
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

