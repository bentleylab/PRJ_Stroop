%% compare photodiode trigger for IR62 SUR Target Time
% Import helper functions to Matlab path
helper_function_dir_name = '/home/knight/hoycw/PRJ_Stroop/scripts/utils/';
addpath(genpath(helper_function_dir_name));


% PLot the two photodiode signals
load ~/PRJ_Stroop/data/IR62/00_raw/NLX_DC1_trigger.mat
time = linspace(0,numel(photodiode.time{1})/photodiode.original_fsample,numel(photodiode.time{1}));
% plot(time,photodiode.trial{1})

cfgl = [];cfgl.continuous = 'yes';cfgl.channel='DC*';
cfgl.dataset = '~/PRJ_Stroop/data/IR62/00_raw/2017071315_0004.besa';
macro = ft_preprocessing(cfgl);
figure;plot(macro.time{1},macro.trial{1}(3,:))


%% Adjust the NLX signal
% Downsample for ease of processing
new_srate = 1000;
decimate_v = floor(photodiode.original_fsample/new_srate);
photo_dec = photodiode.trial{1}(1:decimate_v:end);
plot(photo_dec);

% Set zero/baseline during a block
bsln_val = -32767;% this is the bottom baseline, so setting it to that not the middle (noise level = 0)

% Record epochs (in sec) with fluctuations that should be set to baseline
%   NOTE: bottom to top of one oscillation/block is ~52ms!
%   for simplicity, I'm setting everything to the front bottom corner of the
%   first trial in each block (the first point from which it rises)
bsln_ix = {...
    [1 118270],...%mess before B1, blip at beginning of B1
    [216000 240800],...% extra blip of a trial at the end of B1, gap before B2
    [339000 351300],...%B2 to B3
    [450000 463500],...%B3 to B4
    [562000 866600],...%B4 to B5
    [965000 976700],...%B5 to B6
    [1074000 1131100],...%B6 to B7
    [1230000 1244600],...%B7 to B8
    [1343000 1357700],...%B8 to B9
    [1456000 1581016]...%B9 to end
    };
% Record within trial corrections
%   after adjusting the baseline, there are some sloped rises in onsets that cause double trial flagging
stim_ix = {...
    [138000 139000],...
    [146500 147200],...
    [154375 155312],...
    [168125 169062],...
    [170625 171875],...
    [184375 185625],...
    [187500 188125],...
    [190312 190938],...
    [198000 199000],...
    [206500 207200],...
    [212000 212700],...
    [214500 215500]...    
    };
stim_yval = repmat(32767,size(stim_ix));
if length(stim_ix)~=length(stim_yval)
    error('Number of epochs and values for stimulus correction periods do not match.');
end

%% Apply corrections
photo_fix = photo_dec;
% Correct baseline shift
for shift_n = 1:length(bsln_ix)
    epoch_idx = floor(bsln_ix{shift_n}(1)):floor(bsln_ix{shift_n}(2));
    epoch_idx(epoch_idx<1) = [];
    photo_fix(epoch_idx) = bsln_val;
end

% level out stimulus periods
for stim_n = 1:length(stim_ix)
    epoch_idx = floor(stim_ix{stim_n}(1)):floor(stim_ix{stim_n}(2));
    epoch_idx(epoch_idx<1) = [];
    photo_fix(epoch_idx) = stim_yval(stim_n);
end

%% Find events
photo = photo_fix/1000; %rescale for ease
photo = photo + min(photo); % bring baseline to zero
min_event_length = 1.3*new_srate;

[~, ~, shades] = read_photodiode(photo, min_event_length, 1);  %only 1 shade

shades = [diff(shades) 0]; % Add a point because diff removes one
word_onset = find(shades>0)'; % 1 to 2,3,4 is word onset. Transpose to make column vector
% word_onsets = word_onsets*decimate_v; % convert back to ms
fprintf('\t\tFound %d trials in photodiode channel\n', length(word_onset));

%% Plot events
% Plot photodiode data
% plot_photo = data_photo_orig - min(data_photo_orig);
plot_photo = photo / (max(photo)-min(photo));
plot_photo = plot_photo + 0.25;

% Plot identified events
figure;
plot(plot_photo, 'k'); hold on;
% Plot word onsets
for word_n = 1:length(word_onset)
%     plot([word_onsets(word_n) word_onsets(word_n)],[1.30 1.40],'r','LineWidth',2);
    plot([word_onset(word_n) word_onset(word_n)],[-0.35 0.35],'r','LineWidth',2);
end

% Plot fit against NK/macro photodiode events
load('~/PRJ_Stroop/data/IR62/03_events/IR62_trial_info_auto.mat');
if trial_info.sample_rate~=new_srate
    error('WARNING! sample rates are not equal between NLX/micro and NK/macro word_onsets');
end
% trial_len_NLX = diff(word_onset);
% trial_len_NK  = diff(trial_info.word_onset);

figure;
subplot(2,1,1);
plot(word_onsets,'k');
hold on;
plot(trial_info.word_onset,'r');
legend({'NLX','NK'});
title('Time between trial onsets');

subplot(2,1,2);
plot(word_onsets-trial_info.word_onset);
legend('NLX-NK');
title('Difference between NK and NLX trial onsets');

%% Build full trial_info for NLX
NLX_trial_info = trial_info;
NLX_trial_info.response_time = [];
NLX_trial_info.resp_onset    = [];
NLX_trial_info.marker_time   = [];
NLX_trial_info.rt_window     = [];
NLX_trial_info.ignore_trials = [];
NLX_trial_info.sample_rate   = new_srate;
NLX_trial_info.orig_srate    = photodiode.original_fsample;
NLX_trial_info.word_onset    = word_onset;
NLX_trial_info.word_onset_orig_srate = word_onset*photodiode.original_fsample;

%% Save NLX micro results
save('~/PRJ_Stroop/data/IR62/03_events/IR62_NLX_trial_info_auto.mat','NLX_trial_info');