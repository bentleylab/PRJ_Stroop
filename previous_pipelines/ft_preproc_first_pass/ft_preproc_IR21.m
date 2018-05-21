%% Extract data with fieldtrip pipeline
clear all; close all;

addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

SBJ = 'IR21';%input('Enter subject ID:  ','s');
raw_filename = '2015063014_0002.edf';
% Channel Selection
probe_names = {'ROF*','RAH*','RC*','LIN*','LAH*'};
bad_ch = {...
    'RAH6','RAH7',...%epileptic
    'RC13','RC14','RC15','RC16',...%out of brain
    'DC02','DC03'....% Not real data
    'E','LSh ','LLE','RSh','V1','V2','V3','V4','V5','V6','REF',...% Not real data
    'EKG','C16','EDF Annotations'...
    };
eeg_ch = {};
evnt_ch = {...
    'DC01',...%photodiode
    'DC04',...%microphone
    };
mic_ch = {'DC04'};
% analysis_time = {[65 1375]; %A00 IR21

% View PSDs
plot_psd = 'no';       % 'all', '1by1', or 'no'

% Filtering parameters --> FOR PLOTTING PURPOSES ONLY! Not saved this way
resamp_yn   = 'yes';
demean_yn   = 'yes';
dft_yn      = 'yes';
bs_yn       = 'yes';                % Overrides notch
hp_yn       = 'yes';
lp_yn       = 'yes';
resamp_freq = 1000;
% nothing at 60, most have samll 300 bump and maybe 180
% LIN has 100, 120, 180, 200, 240, 300 Hz, even 360
notch_freqs = [60 100 120 180 200 240 300]; 
bs_width    = 2;
hp_freq     = 0.5;
hp_order    = 4;            % Leaving blank causes instability error, 1 or 2 works 
lp_freq     = 300;
bs_freq_lim = NaN([length(notch_freqs) 2]);
for f_ix = 1:length(notch_freqs)
    bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(notch_freqs(f_ix),bs_width);
end

% basic plotting parameters
% plot_it = 0;
% cfg_plot = [];
% cfg_plot.viewmode = 'vertical'; %trace per line
% cfg_plot.continuous = 'yes'; % not trials
% cfg_plot.plotlabels = 'yes'; %to know what to cut
% cfg_plot.blocksize = 10; % good for bob, can adjust

SBJ_dir = fullfile('/home/knight/hoycw/PRJ_Stroop/data/',SBJ);
raw_dir = fullfile(SBJ_dir, '00_raw/');
raw_filename = fullfile(raw_dir,raw_filename);
import_dir = fullfile(SBJ_dir,'01_import/');
preproc_dir = fullfile(SBJ_dir,'02_preproc/');
if ~exist(import_dir,'dir')
    mkdir(import_dir);
end
if ~exist(preproc_dir,'dir')
    mkdir(preproc_dir);
end

%% Load and preprocess the data
bad_ch     = fn_ch_lab_add_prefix(bad_ch,'POL ');
bad_ch     = fn_ch_lab_add_suffix(bad_ch,'-Ref');
bad_ch_neg = fn_ch_lab_negate(bad_ch);
eeg_ch     = fn_ch_lab_add_prefix(eeg_ch,'POL ');
eeg_ch     = fn_ch_lab_add_suffix(eeg_ch,'-Ref');
eeg_ch_neg = fn_ch_lab_negate(eeg_ch);
evnt_ch     = fn_ch_lab_add_prefix(evnt_ch,'POL ');
evnt_ch     = fn_ch_lab_add_suffix(evnt_ch,'-Ref');
event_ch_neg = fn_ch_lab_negate(evnt_ch);
cfg            = [];
cfg.dataset    = raw_filename;
cfg.continuous = 'yes';
cfg.channel    = {'all',bad_ch_neg{:},eeg_ch_neg{:},event_ch_neg{:},'-EDF Annotations'};
data = ft_preprocessing(cfg);   % just load the data, don't process it

% % EEG data
% cfg.channel = eeg_ch;
% eeg = ft_preprocessing(cfg);   % just load the data, don't process it
% output_filename = strcat(import_dir,SBJ,'_ft_eeg.mat');
% save(output_filename, '-v7.3', 'eeg');
% clear eeg

% Pull out event data
cfg.channel = evnt_ch;
evnt = ft_preprocessing(cfg);   % just load the data, don't process it
% Remove extra characters from channel labels
for channel_n = 1:length(evnt.label)
    evnt.label{channel_n} = strrep(evnt.label{channel_n},'POL ','');
    evnt.label{channel_n} = strrep(evnt.label{channel_n},'-Ref','');
end
output_filename = strcat(import_dir,SBJ,'_ft_evnt.mat');
save(output_filename, '-v7.3', 'evnt');

% Save Microphone data separately as .wav for listening
cfg = [];
cfg.channel = mic_ch;
mic_data = ft_selectdata(cfg,evnt);
mic_data = mic_data.trial{1};
%rescale to prevent clipping, add 0.05 fudge factor
mic_data_rescale = mic_data./(max(abs(mic_data))+0.05);
mic_data_filename = strcat(import_dir,SBJ,'_mic_recording.wav');
audiowrite(mic_data_filename,mic_data_rescale,evnt.fsample);
clear evnt mic_data

% Remove extra characters from channel labels
for channel_n = 1:length(data.label)
    data.label{channel_n} = strrep(data.label{channel_n},'POL ','');
    data.label{channel_n} = strrep(data.label{channel_n},'-Ref','');
end

%% Resample data
if strcmp(resamp_yn,'yes')
    cfg = [];
    cfg.resamplefs = resamp_freq;
    cfg.detrend = 'no';
    data = ft_resampledata(cfg, data);
end
output_filename = strcat(import_dir,SBJ,'_ft_resamp_',num2str(resamp_freq),'.mat');
save(output_filename, '-v7.3', 'data');

%% View resampled data
if strcmp(plot_psd,'all')
    fn_plot_PSD_all(data.trial{1},data.label,data.fsample,1);
elseif strcmp(plot_psd,'1by1')
    fn_plot_PSD_1by1(data.trial{1},data.label,data.fsample);
end

%% Preprocess data
cfg           = [];
cfg.demean    = demean_yn;
cfg.dftfilter = dft_yn; % line noise removal using discrete fourier transform
cfg.dftfreq   = notch_freqs;
cfg.bsfilter  = bs_yn;
cfg.bsfreq    = bs_freq_lim;
cfg.lpfilter  = lp_yn;
cfg.lpfreq    = lp_freq;
cfg.hpfilter  = hp_yn;
cfg.hpfreq    = hp_freq;
cfg.hpfiltord = hp_order;
data = ft_preprocessing(cfg,data);

%% Rerefrence
for d = 1:numel(probe_names)
    cfg = [];
    cfg.channel = ft_channelselection(probe_names{d}, data.label);
    cfg.montage.labelold = cfg.channel;
    cfg.montage.labelnew = strcat(cfg.channel(1:end-1),'-',cfg.channel(2:end));
    cfg.montage.tra = eye(numel(cfg.channel)-1,numel(cfg.channel));
    for r = 1:size(cfg.montage.tra, 1)
        for c = 1:size(cfg.montage.tra, 2)
            if cfg.montage.tra(r,c) == 1; cfg.montage.tra(r,c+1) = -1; end
        end
    end
    data_reref{d} = ft_preprocessing(cfg, data);
end

%Concatenate together again
cfg = [];
% cfg.appendsens = 'yes';
data = ft_appenddata(cfg,data_reref{:});

%% Save data
output_filename = strcat(preproc_dir,SBJ,'_ft_preproc.mat');
save(output_filename, '-v7.3', 'data');

%% Preprocess Event Data
% % % Manually correct bad event traces
% % Fix IR35: drop before B1, blip at end of B1
% % photo_channel_n = 1;
% yval = 2600; %this is the zero/baseline during a block
% set_times = [0.0 11.5 110.5 130.0]; % in sec
% data_evnt(photo_channel_n,1:floor(set_times(2)*header_evnt.sample_rate)) = yval;
% data_evnt(photo_channel_n,floor(set_times(3)*header_evnt.sample_rate):floor(set_times(4)*header_evnt.sample_rate)) = yval;
% save(preproc_filename, 'header_ecog', 'data_ecog', 'header_evnt', 'data_evnt', 'data_id');
% % 
% %save
