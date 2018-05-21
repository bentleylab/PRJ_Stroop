%% Photodiode Trace Cleaning Parameters: IR32
% Mark trials to ignore e.g., interruptions
%   B7T4 and T5 for interruption
ignore_trials = [220 221];% these are last trial before nurse interruption and first trial after

% Set zero/baseline during a block
bsln_val = 633;

% Record epochs (in sec) with fluctuations that should be set to baseline
% IR32:
%   noise at front end, B1T1 might have some trouble due to noise at onset
%   blip at end of B1 as usual
yval_bsln = 633; %this is the zero/baseline during a block
bsln_set_times = {[0.0 13.0], [112.0 113.0], [707.5 787.0]}; % in sec
bsln_times = {...
    [0.0 13.0],...% initial offset
    [112.0 113.0],...% extra trial at end of B1
    [707.5 787.0]...% takes care of nurse interruption near beginning of B7
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

