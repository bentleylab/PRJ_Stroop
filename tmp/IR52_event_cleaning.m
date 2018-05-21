%% IR52 Clean Photodiode Trace
clear all; close all;
SBJ = 'IR52';
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

%% Manually correct photodiode trace
% Mark trials to ignore
% e.g., interruptions
ignore_trials = [];

% Set zero/baseline during a block
bsln_val = 476;

% Record epochs (in sec) with fluctuations that should be set to baseline
bsln_times = {...
    [0.0 5.0],...% initial offset (attaching to screen?
    [105.0 107.0]...% extra blip of a trial at the end of B1
    };
% Record epochs (in sec) when baseline has shifted
bsln_shift_times = {};
% Record amount of shift in baseline for each epoch 
bsln_shift_val = {};
if length(bsln_shift_times)~=length(bsln_shift_val)
    error('Number of epochs and values for baseline shift periods do not match.');
end

% correct baseline shift
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

% Record within trial corrections
% % fix double event detections on T22 and T31 in B5
% extra_times = [570.0 571.5 594.7 596.25];
% extra_yval = [5013 3582]; % inc and neu at those points in time
% data_evnt(photo_channel_n,floor(extra_times(1)*header_evnt.sample_rate):floor(extra_times(2)*header_evnt.sample_rate)) = extra_yval(1);
% data_evnt(photo_channel_n,floor(extra_times(3)*header_evnt.sample_rate):floor(extra_times(4)*header_evnt.sample_rate)) = extra_yval(2);

%% Save corrected data
out_filename = [preproc_dir SBJ '_evnt_clean.mat'];
save(out_filename, 'evnt', 'hdr', 'photod_ix', 'mic_ix', 'ignore_trials');
