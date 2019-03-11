%% Photodiode Trace Cleaning Parameters: CP26
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = 0;% these blocks were demeaned prior to concatenating R1 and R2

% Record epochs (in sec) with fluctuations that should be set to baseline
%   photo is missed sometimes, trying ot get rid of mic noise in baseline to fix it
bsln_times = {...
    [0.0 10.0],...% initial offset (attaching to screen?)
    [196.97 197.2],...%mic in baseline
    [280.15 280.5],...%mic in baseline
    [420.2 420.8],...%mic in baseline
    [490.07 490.6],...%mic in baseline
    [721.9 722.6]...% mic in baseline
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
%   mic appears in the photodiode, so need to flatten stim levels when it causes another event
%   B1T4, B1T29, B1T35, B2T18, B2T29, B3T36, B8T5, B9T2
stim_times = {[19.0 20.0], [87.0 88.0], [103.0 104.0], [164.0 165.0],...
              [193.0 194.0], [319.0 32.0], [786.6 786.8], [885.0 886.2]};
stim_yval = [0.3055 0.3055 0.3055 0.3055 0.3055 0.3055 0.3055 0.3055];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

% photodiode is large, but the difference between stim levels is very small
%   I will fix this here by increasing levels individually
low_idx = evnt(1,:)>0.2 & evnt(1,:)<0.315;
mid_idx = evnt(1,:)>0.315 & evnt(1,:)<0.335;
hi_idx  = evnt(1,:)>0.335;
evnt(1,low_idx) = 1;
evnt(1,mid_idx) = 2;
evnt(1,hi_idx)  = 3;
