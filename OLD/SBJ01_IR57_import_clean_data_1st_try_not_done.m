%% Extract data with fieldtrip pipeline
clear all; close all;

addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

SBJ = 'IR57';%input('Enter subject ID:  ','s');
raw_filename = '2017032319_0023.besa';
% Channel Selection
ch_labels.probe_names = {'RSM','RAC','ROF','RIN','RTI','RAM','RHH','RTH',...
                         'LSMA','LAC','LOF','LIN','LTI','LAM','LHH','LTH'};
ch_labels.all = {
    'DC01' 'DC02' 'DC03' 'DC04'...
    'RSM1' 'RSM2' 'RSM3' 'RSM4' 'RSM5' 'RSM6' 'RSM7' 'RSM8' 'RSM9' 'RSM10'...
    'RAC1' 'RAC2' 'RAC3' 'RAC4' 'RAC5' 'RAC6' 'RAC7' 'RAC8' 'RAC9' 'RAC10'...
    'ROF1' 'ROF2' 'ROF3' 'ROF4' 'ROF5' 'ROF6' 'ROF7' 'ROF8' 'ROF9' 'ROF10'...
    'RIN1' 'RIN2' 'RIN3' 'RIN4' 'RIN5' 'RIN6' 'RIN7' 'RIN8' 'RIN9' 'RIN10'...
    'RTI1' 'RTI2' 'RTI3' 'RTI4' 'RTI5' 'RTI6' 'RTI7' 'RTI8' 'RTI9' 'RTI10'...
    'RAM1' 'RAM2' 'RAM3' 'RAM4' 'RAM5' 'RAM6' 'RAM7' 'RAM8' 'RAM9' 'RAM10'...
    'EKG' 'XREF' 'E'...
    'RHH1' 'RHH2' 'RHH3' 'RHH4' 'RHH5' 'RHH6' 'RHH7' 'RHH8' 'RHH9' 'RHH10'...
    'RTH1' 'RTH2' 'RTH3' 'RTH4' 'RTH5' 'RTH6' 'RTH7' 'RTH8' 'RTH9' 'RTH10'...
    'LSMA1' 'LSMA2' 'LSMA3' 'LSMA4' 'LSMA5' 'LSMA6' 'LSMA7' 'LSMA8' 'LSMA9' 'LSMA10'...
    'LAC1' 'LAC2' 'LAC3' 'LAC4' 'LAC5' 'LAC6' 'LAC7' 'LAC8' 'LAC9' 'LAC10'...
    'LOF1' 'LOF2' 'LOF3' 'LOF4' 'LOF5' 'LOF6' 'LOF7' 'LOF8' 'LOF9' 'LOF10'...
    'LIN1' 'LIN2' 'LIN3' 'LIN4' 'LIN5' 'LIN6' 'LIN7' 'LIN8' 'LIN9' 'LIN10'...
    'LTI1' 'LTI2' 'LTI3' 'LTI4' 'LTI5' 'LTI6' 'LTI7' 'LTI8' 'LTI9' 'LTI10'...
    'LAM1' 'LAM2' 'LAM3' 'LAM4' 'LAM5' 'LAM6' 'LAM7' 'LAM8' 'LAM9' 'LAM10'...
    'LHH1' 'LHH2' 'LHH3' 'LHH4' 'LHH5' 'LHH6' 'LHH7' 'LHH8' 'LHH9' 'LHH10'...
    'LTH1' 'LTH2' 'LTH3' 'LTH4' 'LTH5' 'LTH6' 'LTH7' 'LTH8' 'LTH9' 'LTH10'...
    'FPZ' 'CZ' 'OZ' 'C3' 'C4' 'Z' 'FP1' 'FP2' 'T3' 'T4' 'O1' 'O2'...
    'LUC' 'LLC' 'RUC' 'RLC'...
    };
ch_labels.bad = {...
%     'LHH1','LHH2','LHH3','RHH1','RHH2','RHH3','LTH3','LTH4',...%epileptic
%     'LAC4',...%crazy strong noise
%     'RIN10','LAM8','LAM9','LAM10','LPC9','LPC10','LHH9','LHH10',...%out of brain
%     'LTH1','LTH10','ROF8','ROF9','ROF10','RTH8','RTH9','RTH10',...%out of brain
%     'RHH7','RHH8','RHH9','RHH10','RAM9','RAM10','LIN9','LIN10',...%out of brain
%     'LOF1','LOF8','LOF9','LOF10','LAC11','LAC12',...%out of brain
%     'NULL','NULL-1','NILL','NULL-2','DC03','DC04'....% Not real data
%     'E','LSH','LLE','RSH','V1','V2','V3','V4','V5','V6','xREF',...% Not real data
%     'EKG'...
    };
ch_labels.eeg = {'FPZ' 'CZ' 'OZ' 'C3' 'C4' 'FP1' 'FP2' 'T3' 'T4' 'O1' 'O2'};
ch_labels.evnt = {...
    'DC01',...%photodiode
    'DC02',...%microphone
    };
ch_labels.mic = {'DC02'};

analysis_time = {[105 1140]};       % from A00 IR35

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
% Big bumps at 60, 100, (small)120, 180, 200, 240, 300 Hz, some other small ones
% LAC4 has HUGE noise! evrything is way bigger, new peak at 150
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

SBJ_dir = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/'];
raw_dir = [SBJ_dir '00_raw/'];
raw_filename = [raw_dir raw_filename];
import_dir = [SBJ_dir '01_import/'];
preproc_dir = [SBJ_dir '02_preproc/'];
if ~exist(import_dir,'dir')
    mkdir(import_dir);
end
if ~exist(preproc_dir,'dir')
    mkdir(preproc_dir);
end

%% Load and preprocess the data
bad_ch_neg = fn_ch_lab_negate(bad_ch);
eeg_ch_neg = fn_ch_lab_negate(eeg_ch);
event_ch_neg = fn_ch_lab_negate(evnt_ch);
cfg            = [];
cfg.dataset    = raw_filename;
cfg.continuous = 'yes';
cfg.channel    = {'all',bad_ch_neg{:},eeg_ch_neg{:},event_ch_neg{:}};
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
    data.label{channel_n} = strrep(data.label{channel_n},'POL','');
    data.label{channel_n} = strrep(data.label{channel_n},'Ref','');
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
