%% Photodiode Trace Cleaning Parameters: IR84
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = -115;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 331.0],...% initial offset
    [1330.0 1380.002]...% offset at end 
    };
% Record epochs (in sec) when baseline has shifted
%   continous, linear drift dwon for first 3 blocks, slope fit on 430:440s = -0.4218 using cftools
bsln_shift_line = -0.0002*(1:size(evnt,2)) - 30;
bsln_shift_line(405001:end) = 0; % stop ~810s in block 5
evnt(1,:) = evnt(1,:) - bsln_shift_line;

bsln_shift_times = {[1 810]};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [115];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

