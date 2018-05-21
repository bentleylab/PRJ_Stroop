%% IR57 Clean Photodiode Trace
clear all; close all;
SBJ = 'IR57';
photod_ix = 1;
mic_ix = 2;

%% Load data
SBJ_dir = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/'];
import_dir = [SBJ_dir '01_import/'];
preproc_dir = [SBJ_dir '02_preproc/'];

evnt_filename = strcat(import_dir,SBJ,'_evnt.mat');
load(evnt_filename);
[evnt, hdr] = fn_format_data_ft2KLA(evnt);

%% Plot event channels
plot(linspace(0,hdr.length_in_seconds,hdr.n_samples), evnt);

%% Manually Mark photodiode trace
% Mark trials to ignore
% e.g., interruptions
ignore_trials = [1:5];

% Set zero/baseline during a block
bsln_val = 2710;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 22.75],...% first 4 trials are missing photodiode; onset of 5th is messed up (toss all 5)
    [107.0 109.0]...% extra blip of a trial at the end of B1
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = {};
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% Record within trial corrections
% wavy stimulus shade on T22 and T31 in B5
stim_times = {[22.96 24.5]};
stim_yval = [30005];
if length(stim_times)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

%% Correct baseline shift
for shift_ix = 1:length(bsln_shift_times)
    epoch_idx = floor(bsln_shift_times{shift_ix}(1)*hdr.sample_rate):floor(bsln_shift_times{shift_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = evnt(photod_ix,epoch_idx) - bsln_shift_val(shift_ix);
end
% zero out drifts
for zero_ix = 1:length(bsln_times)
    epoch_idx = floor(bsln_times{zero_ix}(1)*hdr.sample_rate):floor(bsln_times{zero_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = bsln_val;
end

% level out stimulus periods
for stim_ix = 1:length(stim_times)
    epoch_idx = floor(stim_times{stim_ix}(1)*hdr.sample_rate):floor(stim_times{stim_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = stim_yval(stim_ix);
end

%% Save corrected data
out_filename = [preproc_dir SBJ '_evnt_clean.mat'];
save(out_filename, 'evnt', 'hdr', 'photod_ix', 'mic_ix', 'ignore_trials');
