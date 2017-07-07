%% Stroop analysis - frequency analyses
clear all; %close all
addpath(genpath('/home/knight/hoycw/Apps/fieldtrip/'));
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
%addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));
ft_defaults

% Analysis Parameters
SBJ       = 'IR35';
data_id   = strcat(SBJ,'_LAC_ft_KLA');
cond_lab  = {'con', 'neu', 'inc'}; % conditions to average
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [1000, 1000];            % buffers are for cutting time series and then plotting
bsln_it   = 0;
bsln_lim  = [-200, 0];             % relative to stim
filt_it   = 1;                     % 0=no low pass; 1=pre-averaging lp; 2=post-averaging lp
filt_freq = [70 150];                    % low pass frequency
env_it    = 0;                      % take the envelope [0 1]
post_filt_it = 0;
post_filt_freq = 10;

% Plotting parameters
save_fig      = 0;
cond_colors   = {'c', 'k', 'r'};    % colors for cond_lab plotting
block_colors  = {'r', 'g', 'b'};    % colors for [mcon, same, minc]
prop_con_lab  = {'con-mcon', 'con-same', 'con-minc', 'neu-mcon', 'neu-same',...
    'neu-minc', 'inc-mcon', 'inc-same', 'inc-minc'};
fig_type      = 'eps';

% Process parameters
epoch_id = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];
bsln_id  = [num2str(bsln_lim(1)) '.' num2str(bsln_lim(2))];
if filt_it == 0
    bp_id = 'bp.none';
elseif filt_it == 1
    bp_id = ['bp' num2str(filt_freq(1)) '.' num2str(filt_freq(2))];
elseif filt_it == 2
    bp_id = ['bp.avgs.' num2str(filt_freq(1)) '.' num2str(filt_freq(2))];
else
    error('bandpass flag not in [0, 1, 2]');
end
if env_it ==1
    env_id = 'env_';
else
    env_id = '';
end

SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/ERPs/',SBJ,'/');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Load data
load(strcat(SBJ_dir,'04_proc/',data_id,'.mat'));
% need only ok_epochs index, the rest is all in 04_proc file (variables seem same as of 10/10/16...)
load(strcat(SBJ_dir,'03_events/',SBJ,'_trial_info_full.mat'));

% need only ok_epochs index, the rest is all in 04_proc file (variables seem same as of 10/10/16...)
bob_bad_epochs = load(strcat(SBJ_dir,'03_events/',SBJ,'_bob_bad_epochs.mat'),'bad_epochs');
bob_bad_epochs = bob_bad_epochs.bad_epochs;
%         if strcmp(data_id,'IR39_RAC_BP_data') || strcmp(data_id,'IR39_ft_preproc')%!!! SHITTY WORK AROUND!!!
%             KLA_ok_epochs = load(strcat(SBJ_dir,'06_epochs/IR39_RAC_WM_',epoch_id,'.mat'),'ok_epochs');
%         else
%             KLA_ok_epochs = load(strcat(SBJ_dir,'06_epochs/',data_id,'_',epoch_id,'.mat'),'ok_epochs');
%         end
%         KLA_ok_epochs = KLA_ok_epochs.ok_epochs;

% Correct all info to only good trials
if header_ecog.sample_rate~=1000
    error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!\n')
end
% Toss epochs that overlap with bad_epochs from Bob
bad_samples = [];
time_epoch_idx = fn_epoch_cuts_datapad(1:size(data_ecog,2),trial_info.word_onset,...
    trial_info.resp_onset,buff_lim);
for bad_epoch_ix = 1:size(bob_bad_epochs,1)
    bad_samples = [bad_samples bob_bad_epochs(bad_epoch_ix,1):bob_bad_epochs(bad_epoch_ix,2)];
end
bob_ok_epochs = [];
for epoch_ix = 1:size(time_epoch_idx,1)
    if isempty(intersect(time_epoch_idx(epoch_ix,:),bad_samples));
        bob_ok_epochs = [bob_ok_epochs epoch_ix];
    end
end
ok_epochs = bob_ok_epochs;%intersect(bob_ok_epochs,KLA_ok_epochs);
% load(strcat(SBJ_dir,'06_epochs/',data_id,'_',epoch_id,'.mat'),'ok_epochs');

%% Calculate Behavioral Timing
% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
% trial_type = NaN(size(trial_info.condition_n));
no_RT_trials = find(isnan(trial_info.response_time));
ok_epochs    = setdiff(ok_epochs,no_RT_trials);
n_trials_all = length(trial_info.response_time);
word_onset = trial_info.word_onset(ok_epochs);
resp_onset = trial_info.resp_onset(ok_epochs);
cond_n     = trial_info.condition_n(ok_epochs)'; %convert to column vector to match
RTs        = round(1000*trial_info.response_time(ok_epochs)); % converts sec to ms
max_RT     = max(RTs);

%% Convert to fieldtrip
raw_data = [];
raw_data.label = header_ecog.channel_labels;
raw_data.fsample = header_ecog.sample_rate;
raw_data.trial{1} = data_ecog;
raw_data.time{1}  = linspace(0, size(data_ecog,2)/header_ecog.sample_rate, size(data_ecog,2));

%% define trials
cfg = [];
% cfg.dataset = strcat(SBJ_dir,'04_proc/',data_id,'_data_only.mat');
% Define trials as my whole period of interest (buffer to buffer, no offset)
cfg.trl = [word_onset-buff_lim(1), ...        % trial onset
    resp_onset+buff_lim(2),...    % trial offset
    repmat(0,length(ok_epochs),1),...         %baseline before each trial
    cond_n];                 % trial type
% cfg.continuous = 'yes';
trl_data = ft_redefinetrial(cfg, raw_data);

%% Overwrite ft timing variables
for t_ix = 1:length(RTs)
    % Overwrite/redefine time to reflect my trial structure (in sec)
    trl_data.time{t_ix} = trl_data.time{t_ix}-(buff_lim(1)/trl_data.fsample);
    % Overwrite sampleinfo to avoid overlap between trials (fool ft into
    % thinking it's continuous data)
    t_len = length(trl_data.trial{t_ix});
    if t_ix ==1
        trl_data.sampleinfo(t_ix,:) = [1 t_len+1];
    else
        last_end = trl_data.sampleinfo(t_ix-1,2);
        trl_data.sampleinfo(t_ix,:) = [last_end+1 last_end+t_len+1];
    end
end


%% frequency analysis
foi_center = [2 2.5 3 3.5 4 4.5 5 5.5 6 7 8 9 10 11 12 13 14 15 16 18 20 25 30 40 50 60 80 100 125 150 175 200 250];%2.^[2:1/4:8]; % Center frequencies
octave = 3/4;              % Frequency resolution
foi_min = 2^(-octave/2)*foi_center;
foi_max = 2^(octave/2)*foi_center;
foi = (foi_min+foi_max)/2;
delta_freq = foi_max-foi_min;
delta_time = 0.5;
n_taper_all = max(1,round(delta_freq.*delta_time-1));
foi_center  = round(foi_center*10)/10;
delta_freq_true = (n_taper_all+1)./delta_time; % total badn width around

tfr = {};
for cond_ix = 1:2
    cfg = [];
    cfg.output       = 'pow';
    cfg.channel      = 'LAC*';
    cfg.trials       = !!!find(trial_info.condition_n==3);
    cfg.method       = 'mtmconvol';%'wavelet';
    % cfg.width        = 7;
    % cfg.gwidth       = 3;
    cfg.taper        = 'dpss';
    cfg.tapsmofrq    = delta_freq_true./2; %ft wants half bandwidth around the foi
    cfg.keeptapers   = 'no';
    cfg.pad          = 'maxperlen'; %add time on either side of window
    cfg.padtype      = 'zero';
    cfg.foi          = foi_center;%2:5:150;                         % analysis 2 to 30 Hz in steps of 2 Hz
    % cfg.t_ftimwin    = ones(1,length(cfg.tapsmofrq))*delta_time;
    cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
    cfg.toi          = 'all';%-buff_lim(1):0.1:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
    cfg.keeptrials   = 'yes';
    tfr{cond_ix} = ft_freqanalysis(cfg, trl_data);
end

%zscore afterwards
%ft_baselinecorrect
%% Single Plots
cfg_s = [];
cfg_s.latency = [-0.5 2];
tfr_plot = ft_selectdata(cfg_s, tfr);
cfg_lay = [];
cfg_lay.layout = 'vertical';
layout = ft_prepare_layout(cfg_lay,tfr_plot);
conditions = {'con','inc'};
event = 'stim';
for cond_ix = 1:length(conditions)
%     subplot(length(tfr_plot.label),1,cond_ix);
    out_file = [data_id '_TFR_' event '_' conditions{cond_ix}];
%         '_Bob.ep' epoch_id '_bsln' bsln_id env_id smooth_id sig_id y_scale_id];
%     if n_chan > 2
    fig_height = 1;
    vis_fig='on';
%     else
%         fig_height = n_chan/3;
%     end
    figure('Name',out_file);%,'units','normalized',...
%         'outerposition',[0 0 1 fig_height],'Visible',vis_fig);
    cfg = [];
    cfg.trials = logical(fn_condition_index(cond_lab{cond_ix}, cond_n));
    cfg.baseline = [-.2 0];%bsln_lim; %should be in sec
    cfg.baselinetype = 'db';
    cfg.channel = 'all';
    cfg.colorbar = 'yes';
    cfg.zlim = [-3 3];%'maxabs'; %
    cfg.parameter = 'powspctrm';
    cfg.layout = layout;
    % cfg.maskparameter = 'trial';
    % cfg.maskstyle = 'outline';
    cfg.showlabels   = 'yes';
    ft_multiplotTFR(cfg, tfr_plot);
end

%% Single Plots
% cfg_s = [];
% cfg_s.latency = [-0.5 2];
% tfr_plot = ft_selectdata(cfg_s, tfr);
% figure;
% for ch_ix = 1:length(tfr.label)
%     subplot(length(tfr.label),1,ch_ix);
%     cfg = [];
%     cfg.baseline = [-.2 0];%bsln_lim; %should be in sec
%     cfg.baselinetype = 'db';
%     cfg.channel = tfr.label{ch_ix};
%     cfg.colorbar = 'yes';
%     cfg.zlim = [-3 3];%'maxabs'; %
%     cfg.parameter = 'powspctrm';
%     % cfg.maskparameter = 'trial';
%     % cfg.maskstyle = 'outline';
%     cfg.showlabels   = 'yes';
%     ft_singleplotTFR(cfg, tfr_plot);
% end

%%
% cfg.padding = 1;
% cfg.padtype = 'data';
% cfg.continuous = 'yes';
% cfg.channel = 'all';
% % cfg.hpfreq = 0.8; %what kris already did, so hopefully nothing changes
% % cfg.lpfreq = 200;
% cfg.method = 'channel';
% 
% % Remove bad trials
% bad_epochs = [1:n_trials_all];
% cfg.trl(ok_epochs,:) = [];
