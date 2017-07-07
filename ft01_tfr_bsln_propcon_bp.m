%% Stroop analysis - frequency analyses
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
SBJ       = 'IR31';
data_id   = strcat(SBJ,'_LAC_WM');
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

save_it = 0;


% Plotting parameters
save_fig      = 0;
plot_prop_con = 0;
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
load(strcat(SBJ_dir,'06_epochs/',data_id,'_',epoch_id,'.mat'),'ok_epochs');

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

cfg = [];
cfg.output       = 'pow';
cfg.channel      = 'all';
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
tfr = ft_freqanalysis(cfg, trl_data);

%zscore afterwards
%ft_baselinecorrect
%% 
cfg_s = [];
cfg_s.latency = [-0.5 2];
tfr_plot = ft_selectdata(cfg_s, tfr);
figure;
for ch_ix = 1:length(tfr.label)
    subplot(2,1,ch_ix);
    cfg = [];
    cfg.baseline = [-.2 0];%bsln_lim; %should be in sec
    cfg.baselinetype = 'db';
    cfg.channel = tfr.label{ch_ix};
    cfg.colorbar = 'yes';
    cfg.zlim = 'maxabs'; %[-5 5];
    cfg.parameter = 'powspctrm';
    % cfg.maskparameter = 'trial';
    % cfg.maskstyle = 'outline';
    cfg.showlabels   = 'yes';
    ft_singleplotTFR(cfg, tfr_plot);
end

%%
cfg.padding = 1;
cfg.padtype = 'data';
cfg.continuous = 'yes';
cfg.channel = 'all';
% cfg.hpfreq = 0.8; %what kris already did, so hopefully nothing changes
% cfg.lpfreq = 200;
cfg.method = 'channel';

% Remove bad trials
bad_epochs = [1:n_trials_all];
cfg.trl(ok_epochs,:) = [];
%%
ft_line_noise_freqs = line_noise_freqs;
cfg = [];
cfg.dftfilter = 'yes';
cfg.dftfreq = ft_line_noise_freqs;
data = ft_preprocessing(cfg, data);

% Convert back to mat...
data_ecog = data.trial{:};
header_ecog.line_noise_freqs = ft_line_noise_freqs;







%%
% Correct all info to only good trials
n_chan = size(data_ecog,1);
n_trials = length(ok_epochs);
if header_ecog.sample_rate~=1000
    error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!\n')
end

%% Bandpass if bp_it ==1
if filt_it == 1
    for ch_ix = 1:size(data_ecog,1)
        data_ecog(ch_ix,:) = eegfilt(data_ecog(ch_ix,:), header_ecog.sample_rate, filt_freq(1), []);
        data_ecog(ch_ix,:) = eegfilt(data_ecog(ch_ix,:), header_ecog.sample_rate, [], filt_freq(2));
    end
end
if env_it == 1
    for ch_ix = 1:size(data_ecog,1)
        data_ecog(ch_ix,:) = abs(hilbert(data_ecog(ch_ix,:)));
    end
end
if post_filt_it == 1
    for ch_ix = 1:size(data_ecog,1)
        data_ecog(ch_ix,:) = eegfilt(data_ecog(ch_ix,:), header_ecog.sample_rate, [], post_filt_freq);
    end
end
   
%% Format trials
% Cut Trials
trials_stim = NaN([n_chan,n_trials,buff_lim(1)+max_RT+buff_lim(2)+1]);
trials_resp = NaN([n_chan,n_trials,buff_lim(1)+buff_lim(2)+1]);
for ch = 1:n_chan
    for t = 1:n_trials
        trial = data_ecog(ch,word_onset(t)-buff_lim(1):word_onset(t)+max_RT+buff_lim(2));
        % Baseline the data
        bsln_mean = nanmean(trial(buff_lim(1)+bsln_lim(1):buff_lim(1)+bsln_lim(2)));
        bsln_std  = nanstd(trial(buff_lim(1)+bsln_lim(1):buff_lim(1)+bsln_lim(2)));
        trial_norm = (trial-bsln_mean)/bsln_std;
        
        % Cut into trials
        trials_stim(ch,t,:) = trial_norm;
        if ~isnan(resp_onset(t))
            trials_resp(ch,t,:) = trial_norm(RTs(t):RTs(t)+buff_lim(1)+buff_lim(2));
        end
    end
end

%% ERPs by condition
for ix = 1:length(cond_lab)
    % Get binary condition index
    eval(['cond_idx.' cond_lab{ix} ' = fn_condition_index(cond_lab{ix}, cond_n);']);
    % if plot_prop_con, get that label
    % Average trials matching condition index
    eval(['mean_stim.' cond_lab{ix} ' = squeeze(nanmean(trials_stim(:,cond_idx.' cond_lab{ix} '==1,:),2));']);
    eval(['mean_resp.' cond_lab{ix} ' = squeeze(nanmean(trials_resp(:,cond_idx.' cond_lab{ix} '==1,:),2));']);
    % Low pass if required
    if filt_it == 2
        eval(['mean_stim.' cond_lab{ix} ' = eegfilt(mean_stim.' cond_lab{ix} ...
            ', header_ecog.sample_rate, [], lp_freq);']);
        eval(['mean_resp.' cond_lab{ix} ' = eegfilt(mean_resp.' cond_lab{ix} ...
            ', header_ecog.sample_rate, [], lp_freq);']);
    end        
    
end

% mean_diff_i_c_stim = mean_inc_stim-mean_con_stim;
% % mean_diff_mi_mc_stim = mean_minc_stim-mean_mcon_stim;
% 
% mean_diff_i_c_resp = mean_inc_resp-mean_con_resp;
% % mean_diff_mi_mc_resp = mean_minc_resp-mean_mcon_resp;

%% Plot ERPs
x_step = 250;
x_lab_stim = -buff_lim(1):x_step:max_RT+buff_lim(2);
x_lab_resp = -buff_lim(1):x_step:buff_lim(2);

% Plot Condition ERPs
figure;
if n_chan > 2
    fig_height = 1;
else
    fig_height = n_chan/3;
end
set(gcf,'units','normalized','outerposition',[0 0 1 fig_height]);
for ch = 1:n_chan
    % Plot Stim Locked
    subplot(n_chan,2,ch*2-1); hold on;
    for c_ix = 1:length(cond_lab)
        if n_chan==1
            eval(['plot(mean_stim.' cond_lab{c_ix} ',cond_colors{c_ix});']);
        else
            eval(['plot(mean_stim.' cond_lab{c_ix} '(ch,:),cond_colors{c_ix});']);
        end
    end
%     plot(mean_diff_i_c_stim(ch,:),'g');
    line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
    line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
    % Axes and Labels
    ax = gca;
    ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
    ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
    ax.XTickLabel = x_lab_stim;
    for lab_ix = 1:length(cond_lab)
        uscore = '_';       %!!!FIX ME too lazy to figure out single quoptes in side single suqotes
        cond_lab_legend{lab_ix} = eval(['[cond_lab{lab_ix} uscore num2str(size(cond_idx.' cond_lab{lab_ix} '))];']);
    end
    title(strcat('Stim Locked: ',header_ecog.channel_labels(ch)));
    if ch==1
        legend(cond_lab_legend{:},'stim','mean RT','Location','southeast');%,'inc-con'
    end
    
    % Plot RT locked
    subplot(n_chan,2,ch*2); hold on;
    for c_ix = 1:length(cond_lab)
        if n_chan==1
            eval(['plot(mean_resp.' cond_lab{c_ix} ',cond_colors{c_ix});']);
        else
            eval(['plot(mean_resp.' cond_lab{c_ix} '(ch,:),cond_colors{c_ix});']);
        end
    end
    %     plot(mean_diff_i_c_resp(ch,:),'g');
    line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_lim(1)+buff_lim(2)];
    ax.XTick = 0:x_step:buff_lim(1)+buff_lim(2);
    ax.XTickLabel = x_lab_resp;
    title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
    if ch==1
        legend(cond_lab{:},'resp','Location','southeast');%,'inc-con'
    end
end

% Save plots
out_filename = [fig_dir data_id '_ERPs_ep' epoch_id '_bsln' bsln_id '_' env_id bp_id '.' fig_type];
if save_fig ==1
    fprintf('Saving %s\n',out_filename);
    eval(['export_fig ' out_filename]);
end
% saveas(gcf,out_filename);

%% Plot Condition ERPs split by proportion congruency
% !!!
if plot_prop_con == 1
    for ch = 1:n_chan
        figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);
        % Plot Stim Locked
        % Congruent
        subplot(3,2,1); hold on;
        if n_chan==1
            plot(mean_con_mcon_stim,'r');
            plot(mean_con_same_stim,'g');
            plot(mean_con_minc_stim,'b');
            plot(mean_con_stim,'k','LineWidth',3);
        else
            plot(mean_con_mcon_stim(ch,:),'r');
            plot(mean_con_same_stim(ch,:),'g');
            plot(mean_con_minc_stim(ch,:),'b');
            plot(mean_con_stim(ch,:),'k','LineWidth',3);
        end
        line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
        line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
        ax = gca;
        ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
        ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
        ax.XTickLabel = x_lab_stim;
        title(strcat('Stim Locked, Con: ',header_ecog.channel_labels(ch)));
        legend('mcon','same','minc','avg','stim','mean RT','Location','southeast');%,'inc-con'
        
        % Neutral
        subplot(3,2,3); hold on;
        if n_chan==1
            plot(mean_neu_mcon_stim,'r');
            plot(mean_neu_same_stim,'g');
            plot(mean_neu_minc_stim,'b');
            plot(mean_neu_stim,'k','LineWidth',3);
        else
            plot(mean_neu_mcon_stim(ch,:),'r');
            plot(mean_neu_same_stim(ch,:),'g');
            plot(mean_neu_minc_stim(ch,:),'b');
            plot(mean_neu_stim(ch,:),'k','LineWidth',3);
        end
        line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
        line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
        ax = gca;
        ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
        ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
        ax.XTickLabel = x_lab_stim;
        title(strcat('Stim Locked, Neu: ',header_ecog.channel_labels(ch)));
        %     legend('mcon','same','minc','avg','stim','mean RT');%,'inc-con'
        
        % Incongruent
        subplot(3,2,5); hold on;
        if n_chan==1
            plot(mean_inc_mcon_stim,'r');
            plot(mean_inc_same_stim,'g');
            plot(mean_inc_minc_stim,'b');
            plot(mean_inc_stim,'k','LineWidth',3);
        else
            plot(mean_inc_mcon_stim(ch,:),'r');
            plot(mean_inc_same_stim(ch,:),'g');
            plot(mean_inc_minc_stim(ch,:),'b');
            plot(mean_inc_stim(ch,:),'k','LineWidth',3);
        end
        %     plot(mean_diff_i_c_stim(ch,:),'g');
        line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
        line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
        ax = gca;
        ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
        ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
        ax.XTickLabel = x_lab_stim;
        title(strcat('Stim Locked, Inc: ',header_ecog.channel_labels(ch)));
        %     legend('mcon','same','minc','avg','stim','mean RT');%,'inc-con'
        
        % Plot RT locked
        % Congruent
        subplot(3,2,2); hold on;
        if n_chan==1
            plot(mean_con_mcon_resp,'r');
            plot(mean_con_same_resp,'g');
            plot(mean_con_minc_resp,'b');
            plot(mean_con_resp,'k','LineWidth',3);
        else
            plot(mean_con_mcon_resp(ch,:),'r');
            plot(mean_con_same_resp(ch,:),'g');
            plot(mean_con_minc_resp(ch,:),'b');
            plot(mean_con_resp(ch,:),'k','LineWidth',3);
        end
        line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
        ax = gca;
        ax.XLim = [0,buff_lim(1)+buff_lim(2)];
        ax.XTick = 0:x_step:buff_lim(1)+buff_lim(2);
        ax.XTickLabel = x_lab_resp;
        title(strcat('Resp Locked, Con: ',header_ecog.channel_labels(ch)));
        legend('mcon','same','minc','avg','resp','Location','southeast');%,'inc-con'
        
        % Neutral
        subplot(3,2,4); hold on;
        if n_chan==1
            plot(mean_neu_mcon_resp,'r');
            plot(mean_neu_same_resp,'g');
            plot(mean_neu_minc_resp,'b');
            plot(mean_neu_resp,'k','LineWidth',3);
        else
            plot(mean_neu_mcon_resp(ch,:),'r');
            plot(mean_neu_same_resp(ch,:),'g');
            plot(mean_neu_minc_resp(ch,:),'b');
            plot(mean_neu_resp(ch,:),'k','LineWidth',3);
        end
        line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
        ax = gca;
        ax.XLim = [0,buff_lim(1)+buff_lim(2)];
        ax.XTick = 0:x_step:buff_lim(1)+buff_lim(2);
        ax.XTickLabel = x_lab_resp;
        title(strcat('Resp Locked, Neu: ',header_ecog.channel_labels(ch)));
        %     legend('mcon','same','minc','avg','resp');%,'inc-con'
        
        % Incongruent
        subplot(3,2,6); hold on;
        if n_chan==1
            plot(mean_inc_mcon_resp,'r');
            plot(mean_inc_same_resp,'g');
            plot(mean_inc_minc_resp,'b');
            plot(mean_inc_resp,'k','LineWidth',3);
        else
            plot(mean_inc_mcon_resp(ch,:),'r');
            plot(mean_inc_same_resp(ch,:),'g');
            plot(mean_inc_minc_resp(ch,:),'b');
            plot(mean_inc_resp(ch,:),'k','LineWidth',3);
        end
        line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
        ax = gca;
        ax.XLim = [0,buff_lim(1)+buff_lim(2)];
        ax.XTick = 0:x_step:buff_lim(1)+buff_lim(2);
        ax.XTickLabel = x_lab_resp;
        title(strcat('Resp Locked, Inc: ',header_ecog.channel_labels(ch)));
        %     legend('mcon','same','minc','avg','resp');%,'inc-con'
        
        out_filename = [fig_dir data_id '_ERPs_' num2str(ch) '_ep' epoch_id ...
            '_propcon_bsln' bsln_id '_' bp_id '.' fig_type];
        fprintf('Saving %s\n',out_filename);
        eval(['export_fig ' out_filename]);
        %  saveas(gcf,out_filename);
    end
end