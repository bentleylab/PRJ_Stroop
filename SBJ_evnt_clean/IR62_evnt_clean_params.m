%% Photodiode Trace Cleaning Parameters: IR62
% Mark trials to ignore e.g., interruptions
ignore_trials = [];  % see photod comment below

% Set zero/baseline during a block
bsln_val = -366;

% Record epochs (in sec) with fluctuations that should be set to baseline
%	drop before B1, blip at end of B1
%   the baseline fluctuations become consistently up instead of down after B4,
%   which causes the algorithm to mistake the baseline with the smallest level
bsln_times = {...
    [0.0 114.0],...
    [213.0 215.0]...% extra blip of a trial at the end of B1
    };
%    [558.0 863.6]...% gap between B4 and B5 with high fluctuation in baseline

% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {[602.5 1589.2]};
% Record amount of shift in baseline for each epoch 
%   this shift casues mixing of short and medium levels, but event onsets are mostly good
%   since the actual trial typ eis taken from the log file and not the photodidoe level, this is fine
bsln_shift_val = [366.3];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
%   after adjusting the baseline, there are some sloped rises in onsets that cause double trial flagging
stim_times = {...
    [942.69 944.27],...%B5T30
    [951.0 952.5],...%B5T33
    [1202.3 1203.8],...%B7T28
    [1247.2 1248.6],...%B8T3
    [1447.7 1449.2]...%B9T35
    };
stim_yval = [3663, 3663, 3663, 3663, 3663];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

