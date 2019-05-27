if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(genpath('/Users/colinhoy/Code/Apps/wave_clus/'));
addpath(genpath('/Users/colinhoy/Code/Apps/UR_NLX2MAT_releaseDec2015/'));
addpath(ft_dir);
ft_defaults

%% SBJ variables
SBJ = 'IR69';
block_suffix = '_R1';

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

SBJ_evnt_clean_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_evnt_clean/' SBJ '_evnt_clean_params',block_suffix,'.m'];
eval(SBJ_evnt_clean_cmd);


%% Read photodiode
photo_label = 'Photo1_0010';
inverted = 1;
nlx_dir = [root_dir 'PRJ_Stroop/data/IR69/00_raw/SU_2018-02-07_13-17-17/'];
photo = ft_read_neuralynx_interp({[nlx_dir 'photo/' photo_label '.ncs']});
photo.label = {'Photo1'};
if inverted
    photo.trial{1} = photo.trial{1}*-1;
end

%% Minimal processing for CPPC 2018 Abstract
% Plot original
plot(photo.time{1},photo.trial{1})

% Correct linear drift, fit beginning from 0:106s (before drop to task baseline)
cfgs = [];
cfgs.latency= [0 106];
photo_baseline = ft_selectdata(cfgs,photo);
photo_bsln_data = photo_baseline.trial{1};
% cftool
% fit and remove linear model: f(x) = p1*x + p2, 
p1 = 3.014e-05;
p2 = 2659;
drift = p1*[0:numel(photo.trial{1})-1] + p2;
photo.trial{1} = photo.trial{1}-drift;

% Scale last 2 blocks to match first 7
%   Photodiode came loose in last trial of B7
%   B7T36 was inc (middle), but drop makes level close to neu (shortest)
% B7T36 is last trial before drop in size: 841-842s
% B8T1 is first of new scale: 859-860s
init_stim = [find(photo.time{1}==841) find(photo.time{1}==842)];
next_stim = [find(photo.time{1}==859) find(photo.time{1}==860)];
max_init = min(photo.trial{1}(init_stim(1):init_stim(2)));  % used min instead of mean because mean/mean ratio made it a little too big
max_next = mean(photo.trial{1}(next_stim(1):next_stim(2)));

next_times = [844.05 1070];
next_ix = [find(photo.time{1}==next_times(1)) find(photo.time{1}==next_times(2))];
photo.trial{1}(next_ix(1):next_ix(2)) = photo.trial{1}(next_ix(1):next_ix(2))*max_init/max_next;
% This gets the proportions about right, but the whole chunk of data is too
% low and needs a baseline boost
%   bsln between B6/7 = 740:745
%   bsln between B7/8 = 847:855
init_bsln_ix = [find(photo.time{1}==740) find(photo.time{1}==745)];
next_bsln_ix = [find(photo.time{1}==847) find(photo.time{1}==855)];
init_bsln = mean(photo.trial{1}(init_bsln_ix(1):init_bsln_ix(2)));
next_bsln = mean(photo.trial{1}(next_bsln_ix(1):next_bsln_ix(2)));
photo.trial{1}(next_ix(1):next_ix(2)) = photo.trial{1}(next_ix(1):next_ix(2))+init_bsln-next_bsln;

% Skipping cut to analysis time for now... (don't have to adjust spike times)

%% Create correction times and values in a separate file in ~/PRJ_Stroop/scripts/SBJ_evnt_clean/
% Convert to KLA format
[evnt, hdr] = fn_format_data_ft2KLA(photo);
photod_ix = 1;

% Correct baseline shift
for shift_ix = 1:length(bsln_shift_times)
    epoch_idx = floor(bsln_shift_times{shift_ix}(1)*hdr.sample_rate):floor(bsln_shift_times{shift_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = evnt(photod_ix,epoch_idx) - bsln_shift_val(shift_ix);
end
% zero out drifts
for zero_ix = 1:length(bsln_times)
    if any(bsln_times{zero_ix}==-1) % Cheating and just making it zero to the end manually... won't be a problem when I use analysis_time.
        evnt(photod_ix,floor(bsln_times{zero_ix}(1)*hdr.sample_rate):end) = bsln_val;
    else
        epoch_idx = floor(bsln_times{zero_ix}(1)*hdr.sample_rate):floor(bsln_times{zero_ix}(2)*hdr.sample_rate);
        epoch_idx(epoch_idx<1) = [];
        evnt(photod_ix,epoch_idx) = bsln_val;
    end
end

% level out stimulus periods
for stim_ix = 1:length(stim_times)
    epoch_idx = floor(stim_times{stim_ix}(1)*hdr.sample_rate):floor(stim_times{stim_ix}(2)*hdr.sample_rate);
    epoch_idx(epoch_idx<1) = [];
    evnt(photod_ix,epoch_idx) = stim_yval(stim_ix);
end

%% Cut to analysis time
tmp_photo = ft_read_neuralynx_interp({[SBJ_vars.dirs.SU 'photo/' SBJ_vars.ch_lab.photod{1} SBJ_vars.ch_lab.suffix '.ncs']});
tmp_photo.trial{1} = evnt;
cfgs = []; cfgs.latency = SBJ_vars.analysis_time{b_ix}{1};
tmp_photo = ft_selectdata(cfgs,tmp_photo);
tmp_photo.time{1} = tmp_photo.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);

%% Read, downsample, and save mic
mic = ft_read_neuralynx_interp({[SBJ_vars.dirs.SU 'mic/' SBJ_vars.ch_lab.mic{1} SBJ_vars.ch_lab.suffix '.ncs']});
mic.label = {'mic'};

% Cut to analysis time
mic = ft_selectdata(cfgs,mic);
mic.time{1} = mic.time{1}-SBJ_vars.analysis_time{b_ix}{1}(1);

% Downsample
if proc_vars.mic_resample
    cfg_dsmp = [];
    cfg_dsmp.resamplefs = photo.fsample;
    mic_dsmp = ft_resampledata(cfg_dsmp,mic);
end

%% Process Microphone data
mic_data = mic.trial{1};
%rescale to prevent clipping, add 0.05 fudge factor
mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);

%% Concatenate photo and mic
% Check that time vectors are close enough
dif = tmp_photo.time{1}-mic_dsmp.time{1};
if max(dif)>0.000001
    error('time vectors not aligned, check that!');
elseif ~isequal(photo.time{1},mic_dsmp.time{1})
    warning(['Mic and photo time vectors still not equal, but differences is small: ' num2str(max(dif))]);
    mic_dsmp.time{1} = tmp_photo.time{1};
end

% Append
evnt = vertcat(tmp_photo.trial{1},mic_dsmp.trial{1});
photod_ix = 1;
mic_ix = 2;
ignore_trials = [];
% cfga = [];
% cfga.keepsampleinfo = 'no';
% evnt = ft_appenddata(cfga,photo,mic_dsmp);
% evnt.fsample = photo.fsample;

%% Save data out
evnt_out_filename = strcat(SBJ_vars.dirs.import,SBJ,'_evnt',block_suffix,'.mat');
save(evnt_out_filename, '-v7.3', 'evnt');

% Save Microphone data separately as .wav for listening
mic_data_filename = strcat(SBJ_vars.dirs.import,SBJ,'_mic_recording',block_suffix,'.wav');
audiowrite(mic_data_filename,mic_data_rescale,evnt.fsample);

%% Save corrected data
out_filename = [SBJ_vars.dirs.preproc SBJ '_evnt_clean',block_suffix,'.mat'];
save(out_filename, 'evnt', 'hdr', 'ignore_trials', 'photod_ix');%, 'mic_ix'
