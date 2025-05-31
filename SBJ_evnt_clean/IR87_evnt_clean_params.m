%% Photodiode Trace Cleaning Parameters: IR87
% Mark trials to ignore e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = 0;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 262.0],...% initial offset
    [1300.0 1700.0]...% offset at end 
    };
% Record epochs (in sec) when baseline has shifted
%   continous, linear drift dwon for first 3 blocks, slope fit on 430:440s = -0.4218 using cftools
time_vec = linspace(0,hdr.length_in_seconds,hdr.n_samples);
line_times = [0 1045; 1045 1300];
line_ix(1,1) = 1;
[~,line_ix(1,2)] = min(abs(time_vec-line_times(1,2)));
[~,line_ix(2,1)] = min(abs(time_vec-line_times(2,1)));
[~,line_ix(2,2)] = min(abs(time_vec-line_times(2,2)));
% ln_time_vec{1} = (line_times(1,1):1/proc.resample_freq:line_times(1,2));
bsln_shift_line1 = -0.115*time_vec(line_ix(1,1):line_ix(1,2)) + 67;
% xvals{2} = (line_times(2,1):1/proc.resample_freq:line_times(2,2));
bsln_shift_line2 = -0.045*time_vec(line_ix(2,1):line_ix(2,2)) - 8;
evnt(1,line_ix(1,1):line_ix(1,2)) = evnt(1,line_ix(1,1):line_ix(1,2)) - bsln_shift_line1;
evnt(1,line_ix(2,1)+1:line_ix(2,2)) = evnt(1,line_ix(2,1)+1:line_ix(2,2)) - bsln_shift_line2(1:end-1);

% This is still too messy (non-linear drift), so binarizing
evnt(1,evnt(1,:)<18) = bsln_val;

bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = [];
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
%   B4T9, B4T21, B7T33, B7T35
% 2nd round: B7T5, B7T17 (oops, corrected B7T6)
% 3rd round: B7T5, B7T7
stim_times = {[662.5 664.2],[695.16 696.8],[1070.2 1071.6],[1076 1077],...
    [997.52 999.17],[1027.16 1028.81],...
    [994.5 996.54],[1000.17 1001.82]};
stim_yval = [41,41,41,30,...
    30,41,...
    38,38];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

