%% Photodiode Trace Cleaning Parameters: CP24
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = -6.2;% these blocks were demeaned prior to concatenating R1 and R2

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 10.0],...% initial offset (attaching to screen?)
    [108.0 110.0],...% extra blip of a trial at the end of B1
    [1045.0 1059.0],...% initial offset of R2
    [1158.0 1160.0]...% extra blip after R2B1
    };
% Record epochs (in sec) when baseline has shifted
%   continous, linear drift upwards, slope fit on 500:582s = 0.006421 using cftools
%bsln_shift_line = 0.006421*(1:size(data.trial{1},2));
bsln_shift_times = {[1048.001 1383.0]};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [-2.518];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

