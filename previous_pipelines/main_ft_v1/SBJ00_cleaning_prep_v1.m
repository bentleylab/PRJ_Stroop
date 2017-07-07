%% SBJ00_cleaning_prep
% Load raw data, extract channel labels, downsample, bandstop line noise,
% and save the result for visual cleaning inspection.
function SBJ00_cleaning_prep(SBJ,raw_file,plot_psd,block_prefix)
% Inputs:
%   SBJ [str]- name of the subject
%   raw_file [str]- name of the raw data file (not including path)
%   plot_psd [str]- name of function to plot+save PSDs of raw channels
%       'all'- plots all PSDs on one figure
%       '1by1'- plots each PSD on a separate figure
%       any other string- nothing is plotted or saved
%   block_prefix [str]- name of the block, e.g., '_B1'
%       used if there are more than one raw data files
% Outputs:
%   SBJ_preclean.mat- the resampled and bandstop filtered data for visual cleaning
%   SBJ_raw_labels.mat- the raw labels for all channels, to inform channel selection
%   raw_psds/- PSDs of channels as determined by plot_psd

addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

% plot_psd      = '1by1';         % 'all','1by1','none'
psd_fig_type  = 'jpg';
resamp_it     = 1;
resamp_freq   = 1000;
filter_it     = 1;
notch_freqs   = [60 120 180 240 300];
bs_width      = 2;
bs_freq_lim   = NaN([length(notch_freqs) 2]);
for f_ix = 1:length(notch_freqs)
    bs_freq_lim(f_ix,:) = fn_freq_lim_from_CFBW(notch_freqs(f_ix),bs_width);
end

%% Processing
SBJ_dir = ['/home/knight/hoycw/PRJ_Stroop/data/' SBJ '/'];
raw_dir = [SBJ_dir '00_raw/'];
fprintf('============== Processing %s, %s ==============\n',SBJ,raw_file);

%% Set up directories
raw_filename = [raw_dir raw_file];
% Assume only one block for now...
%         if length(SBJ_files{SBJ_ix}{2})>1
%             block_prefix = strcat('_B',num2str(file_ix));
psd_dir = strcat(SBJ_dir,'01_import/raw_psds/');
if ~exist(psd_dir,'dir')
    mkdir(psd_dir);
end
out_dir = strcat(SBJ_dir,'02_preproc/');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

%% Load the data
cfg = [];
cfg.dataset = raw_filename;
cfg.continuous = 'yes';
cfg.channel = 'all';
data = ft_preprocessing(cfg);

%% Check noise profile
if strcmp(plot_psd,'all')
    fprintf('============== Plotting PSDs %s, %s ==============\n',SBJ,raw_file);
    fn_plot_PSD_all_save(data.trial{1},data.label,data.fsample,...
        strcat(psd_dir,SBJ,'_raw_psd',block_prefix,'_all.',psd_fig_type));
elseif strcmp(plot_psd,'1by1')
    fprintf('============== Plotting PSDs %s, %s ==============\n',SBJ,raw_file);
    fn_plot_PSD_1by1_save(data.trial{1},data.label,data.fsample,...
        strcat(psd_dir,SBJ,'_raw_psd',block_prefix),psd_fig_type);
end

%% Resample data to speed things up
if (resamp_it) && (data.fsample ~= resamp_freq)
    cfg_resamp = [];
    cfg_resamp.resamplefs = resamp_freq;
    cfg_resamp.detrend = 'yes';
    data = ft_resampledata(cfg_resamp, data);
end

%% Filter data for ease of viewing
fprintf('============== Bandstop filtering %s, %s ==============\n',SBJ,raw_file);
if (filter_it)% && (data_resamp.cfg.dftfreq ~= notch_freqs)
    cfg_bs = [];
    cfg_bs.continuous = 'yes';
    cfg_bs.bsfilter = 'yes';
    cfg_bs.bsfreq = bs_freq_lim;
    cfg_bs.demean = 'yes';
    data = ft_preprocessing(cfg_bs, data);
end

%% Save data
% Save preclean data
out_filename = strcat(out_dir,SBJ,'_preclean',block_prefix,'.mat');
fprintf('============== Saving %s, %s ==============\n',out_filename,raw_file);
save(out_filename, 'data');

% Save data labels
raw_labels = data.label;
label_filename = strcat(SBJ_dir,'01_import/',SBJ,'_raw_labels',block_prefix,'.mat');
% if raw_labels exit, load, compare, save if different
if ~exist(label_filename)
    save(label_filename,'raw_labels');
else
    error(strcat('ERROR: raw_labels filename already exists: ',label_filename));
end

end