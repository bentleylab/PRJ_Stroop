%% Stroop analysis - power time series
clear all; close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));

% Dataset Parameters
SBJ = 'S02IR21';
ROI = 'RAC';
conditions = {'con', 'neu', 'inc'};
cond_colors = {'c', 'k','r'};

% Filtering Parameters
fbands = {'delta','theta','theta-alpha','alpha','betaL','betaH','HG'};

bsln       = 'demean';      % 'zscore', 'demean', 'none'
buff_pre   = 750;           %buffers are for cutting time series and then plotting
buff_post  = 750;
bsln_start = -200;          %relative to stim
bsln_end   = 0;             %relative to stim

% Plotting parameters
fig_type = 'eps';
fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/pow_ts/',SBJ,'/');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Load data, get basic parameters
load(strcat('/home/knight/kla/Desktop/stroop/',SBJ,'/03_fulldata/',SBJ,'_',ROI,'.mat'));

n_ch = size(data_ecog,1);
n_trl = size(trial_info.condition_n,2);
srate = header_ecog.sample_rate;
if srate~=1000
    error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!')
end
RTs = round(1000*trial_info.response_time); % converts sec to ms
max_RT = max(RTs);

%% Filter full time series
data_env = NaN([length(fbands), n_ch, size(data_ecog,2)]);

for ch_ix = 1:n_ch
    curr_ch = header_ecog.channel_labels{ch_ix};
    disp(strcat('Calculating envelopes for elec=',curr_ch,...
        ' (',num2str(ch_ix),'/',num2str(n_ch),')'));
    
    for fband = 1:length(fbands)
        % Filter elec data
        bp_lim = func_bp_lim(fbands{fband});
        disp(strcat(fbands{fband},'[lo,hi] = ',num2str(bp_lim(1)),'-',num2str(bp_lim(2))));
        ch_bp = func_EEGlab_bandpass(data_ecog(ch_ix,:), srate, bp_lim(1), bp_lim(2));
%         ch_hp = eegfilt(data_ecog(ch_ix,:), srate, bp_lim(1), []);
%         ch_bp = eegfilt(ch_hp, srate, [], bp_lim(2));
        data_env(fband,ch_ix,:) = abs(hilbert(ch_bp));
    end
end
%% Cut trials, baseline normalization
env_trl_stim = NaN([length(fbands), n_ch, n_trl, buff_pre+max_RT+buff_post+1]);
env_trl_resp = NaN([length(fbands), n_ch, n_trl, buff_pre+buff_post+1]);
for fband = 1:length(fbands)
    for ch = 1:n_ch
        for t = 1:n_trl
            trl = data_env(fband,ch,...
                trial_info.word_onset(t)-buff_pre:trial_info.word_onset(t)+max_RT+buff_post);
            % Baseline the data
            bsln_mean = nanmean(trl(buff_pre+bsln_start:buff_pre+bsln_end));
            bsln_std  = nanstd(trl(buff_pre+bsln_start:buff_pre+bsln_end));
            if strcmp(bsln, 'zscore')
                trl_norm = (trl-bsln_mean)/bsln_std;
            elseif strcmp(bsln, 'demean')
                trl_norm = trl-bsln_mean;
            else
                trl_norm = trl;
            end
            % Cut into trials
            env_trl_stim(fband,ch,t,:) = trl_norm;
            if ~isnan(trial_info.resp_onset(t))
                env_trl_resp(fband,ch,t,:) = trl_norm(RTs(t):RTs(t)+buff_pre+buff_post);
            end
        end
    end
end

%% Plotting Power Time Series
% Get condition for each trial
trl_cond_ix = func_kla_cond_ix(trial_info, conditions);

% Average power time seires per condition
for cond = 1:length(conditions)
    eval(['avg_env_stim.' conditions{cond} ' = nanmean(env_trl_stim(:,:,trl_cond_ix. ' ...
        conditions{cond} ',:),3);']);
    eval(['avg_env_resp.' conditions{cond} ' = nanmean(env_trl_resp(:,:,trl_cond_ix. ' ...
        conditions{cond} ',:),3);']);
end

% Plotting parameters
x_step = 250;
x_lab_stim = -buff_pre:x_step:max_RT+buff_post;
x_lab_resp = -buff_pre:x_step:buff_post;

% Plot Power times series by condition, one figure per frequency band
for fband = 1:length(fbands)
    figure;suptitle(strcat(SBJ,',',ROI,',',fbands{fband}));
    if n_ch > 2
        fig_height = 1;
    else
        fig_height = n_ch/3;
    end
    set(gcf,'units','normalized','outerposition',[0 0 1 fig_height]);
    
    for ch = 1:n_ch
        % Plot Stim Locked
        subplot(n_ch,2,ch*2-1); hold on;
        for cond = 1:length(conditions)
            eval(['plot(squeeze(avg_env_stim.' conditions{cond} '(fband,ch,:)),cond_colors{cond});']);
        end
        %     plot(mean_diff_i_c_stim(ch,:),'g');
        line([buff_pre, buff_pre],ylim);%, 'k--');
        line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
        ax = gca;
        ax.XLim = [0,buff_pre+max_RT+buff_post];
        ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
        ax.XTickLabel = x_lab_stim;
        title(strcat('Stim Locked: ',header_ecog.channel_labels(ch)));
        if ch==1
            legend({conditions{:},'stim','mean RT'},'Location','southeast');%,'inc-con'
        end
        
        % Plot RT locked
        subplot(n_ch,2,ch*2); hold on;
        for cond = 1:length(conditions)
            eval(['plot(squeeze(avg_env_resp.' conditions{cond} '(fband,ch,:)),cond_colors{cond});']);
        end
        %     plot(mean_diff_i_c_resp(ch,:),'g');
        line([buff_pre, buff_pre],ylim);%, 'k--');
        ax = gca;
        ax.XLim = [0,buff_pre+buff_post];
        ax.XTick = 0:x_step:buff_pre+buff_post;
        ax.XTickLabel = x_lab_resp;
        title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
        if ch==1
            legend({conditions{:},'resp'},'Location','southeast');%,'inc-con'
        end
    end
    
    % Save plots
    out_filename = [fig_dir SBJ '_' ROI '_pow_' fbands{fband} '_kla_bsln_' bsln '_' num2str(bsln_start) ...
        '.' num2str(bsln_end) '.' fig_type];
    eval(['export_fig ' out_filename]);
%     saveas(gcf,out_filename);
end


