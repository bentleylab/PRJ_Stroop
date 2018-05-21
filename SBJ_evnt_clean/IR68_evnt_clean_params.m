%% Photodiode Trace Cleaning Parameters: IR54
% Mark trials to ignore e.g., interruptions
%   mic died after B5T18 until B6T17, then again B6T29:36
%   B6T17:20 may be too noisy to use...
ignore_trials = [163:196 209:216];

% Set zero/baseline during a block
bsln_val = 17586;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 16.0],...% initial offset (attaching to screen?)
    [115.0 117.0],...% extra blip of a trial at the end of B1
    [568.0 676.0],...% first mic drop out
    [709.0 750.0]...% second mic drop out
    };
% Record epochs (in sec) when baseline has shifted
%   continous, linear drift upwards, slope fit on 500:582s = 0.006421 using cftools
%bsln_shift_line = 0.006421*(1:size(data.trial{1},2));
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

