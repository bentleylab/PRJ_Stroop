%% Photodiode Trace Cleaning Parameters: IR67
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Large drift in the shift over time:
%   plot(evnt(1,:)); Tools --> Basic Fitting --> linear, show equations:
lin_fit = -0.00022997*[1:size(evnt,2)] - 350.75 ;
photod_ix = strmatch(SBJ_vars.ch_lab.photod,hdr.channel_labels);
evnt(photod_ix,:) = evnt(photod_ix,:)-lin_fit;

% Set zero/baseline during a block
bsln_val = 75.885; % this is after linear correction

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 13.0],...% initial offset (attaching to screen?)
    [112.0 114.0],...% extra blip of a trial at the end of B1
    [1139.0 1142.0]...%end of task when nlx was off
    };
% fix last data point
evnt(photod_ix,end) = bsln_val;

% Record epochs (in sec) when baseline has shifted
%   drift throughout: cftools says -2.958e-05 slope, -297.2 intercept (ignored)
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
%   B7T36 (inc/medium) was when photodiode came loose (and caused scaling)
stim_times = {};
stim_yval = [];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

