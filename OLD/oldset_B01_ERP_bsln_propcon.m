%% Stroop analysis - frequency analyses
clear all; close all
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));

SBJ       = 'IR21';
data_id   = strcat(SBJ,'_RC_WM');
epoch_lim = [-200 2000];
fig_type  = 'eps';

epoch_id = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];
SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/ERPs/',SBJ,'/');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Load data
load(strcat(SBJ_dir,'04_proc/',data_id,'.mat'));
% need only ok_epochs index, the rest is all in 04_proc file (versions seem same as of 10/10/16...)
load(strcat(SBJ_dir,'06_epochs/',data_id,'_',epoch_id,'.mat'),'ok_epochs');

% Correct all info to only good trials
n_chan = size(data_ecog,1);
n_trials = length(ok_epochs);
if header_ecog.sample_rate~=1000
    error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!\n')
end
% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
% trial_type = NaN(size(trial_info.condition_n));
word_onset = trial_info.word_onset(ok_epochs);
resp_onset = trial_info.resp_onset(ok_epochs);
cond_n     = trial_info.condition_n(ok_epochs);
RTs        = round(1000*trial_info.response_time(ok_epochs)); % converts sec to ms
max_RT     = max(RTs);


%% Format trials
% Parameters
buff_pre   = 500;   %buffers are for cutting time series and then plotting
buff_post  = 500;
bsln_start = -200;  %relative to stim
bsln_end   = 0;     %relative to stim

% Cut Trials
trials_stim = NaN([n_chan,n_trials,buff_pre+max_RT+buff_post+1]);
trials_resp = NaN([n_chan,n_trials,buff_pre+buff_post+1]);
for ch = 1:n_chan
    for t = 1:n_trials
        trial = data_ecog(ch,word_onset(t)-buff_pre:word_onset(t)+max_RT+buff_post);
        % Baseline the data
        bsln_mean = nanmean(trial(buff_pre+bsln_start:buff_pre+bsln_end));
        bsln_std  = nanstd(trial(buff_pre+bsln_start:buff_pre+bsln_end));
        trial_norm = (trial-bsln_mean)/bsln_std;
        
        % Cut into trials
        trials_stim(ch,t,:) = trial_norm;
        if ~isnan(resp_onset(t))
            trials_resp(ch,t,:) = trial_norm(RTs(t):RTs(t)+buff_pre+buff_post);
        end
    end
end

%% ERPs by condition
con_idx = [cond_n==1]+[cond_n==2]+[cond_n==3];
neu_idx = [cond_n==4]+[cond_n==5]+[cond_n==6];
inc_idx = [cond_n==7]+[cond_n==8]+[cond_n==9];

mcon_idx = [cond_n==3]+[cond_n==6]+[cond_n==9];
same_idx = [cond_n==2]+[cond_n==5]+[cond_n==8];
minc_idx = [cond_n==1]+[cond_n==4]+[cond_n==7];

% Calculate stimulus locked ERP
mean_con_stim = squeeze(nanmean(trials_stim(:,con_idx==1,:),2));
mean_neu_stim = squeeze(nanmean(trials_stim(:,neu_idx==1,:),2));
mean_inc_stim = squeeze(nanmean(trials_stim(:,inc_idx==1,:),2));
mean_diff_i_c_stim = mean_inc_stim-mean_con_stim;

mean_con_mcon_stim = squeeze(nanmean(trials_stim(:,cond_n==3,:),2));
mean_con_same_stim = squeeze(nanmean(trials_stim(:,cond_n==2,:),2));
mean_con_minc_stim = squeeze(nanmean(trials_stim(:,cond_n==1,:),2));
mean_neu_mcon_stim = squeeze(nanmean(trials_stim(:,cond_n==6,:),2));
mean_neu_same_stim = squeeze(nanmean(trials_stim(:,cond_n==5,:),2));
mean_neu_minc_stim = squeeze(nanmean(trials_stim(:,cond_n==4,:),2));
mean_inc_mcon_stim = squeeze(nanmean(trials_stim(:,cond_n==9,:),2));
mean_inc_same_stim = squeeze(nanmean(trials_stim(:,cond_n==8,:),2));
mean_inc_minc_stim = squeeze(nanmean(trials_stim(:,cond_n==7,:),2));
% mean_diff_mi_mc_stim = mean_minc_stim-mean_mcon_stim;

% Calculate RT locked ERP
mean_con_resp = squeeze(nanmean(trials_resp(:,con_idx==1,:),2));
mean_neu_resp = squeeze(nanmean(trials_resp(:,neu_idx==1,:),2));
mean_inc_resp = squeeze(nanmean(trials_resp(:,inc_idx==1,:),2));
mean_diff_i_c_resp = mean_inc_resp-mean_con_resp;

mean_con_mcon_resp = squeeze(nanmean(trials_resp(:,cond_n==3,:),2));
mean_con_same_resp = squeeze(nanmean(trials_resp(:,cond_n==2,:),2));
mean_con_minc_resp = squeeze(nanmean(trials_resp(:,cond_n==1,:),2));
mean_neu_mcon_resp = squeeze(nanmean(trials_resp(:,cond_n==6,:),2));
mean_neu_same_resp = squeeze(nanmean(trials_resp(:,cond_n==5,:),2));
mean_neu_minc_resp = squeeze(nanmean(trials_resp(:,cond_n==4,:),2));
mean_inc_mcon_resp = squeeze(nanmean(trials_resp(:,cond_n==9,:),2));
mean_inc_same_resp = squeeze(nanmean(trials_resp(:,cond_n==8,:),2));
mean_inc_minc_resp = squeeze(nanmean(trials_resp(:,cond_n==7,:),2));
% mean_mcon_resp = squeeze(nanmean(trials_resp(:,mcon_idx==1,:),2));
% mean_same_resp = squeeze(nanmean(trials_resp(:,same_idx==1,:),2));
% mean_minc_resp = squeeze(nanmean(trials_resp(:,minc_idx==1,:),2));
% mean_diff_mi_mc_resp = mean_minc_resp-mean_mcon_resp;

%% Plot ERPs
x_step = 250;
x_lab_stim = -buff_pre:x_step:max_RT+buff_post;
x_lab_resp = -buff_pre:x_step:buff_post;

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
    if n_chan==1
        plot(mean_con_stim,'b');
        plot(mean_neu_stim,'k');
        plot(mean_inc_stim,'r');
    else
        plot(mean_con_stim(ch,:),'b');
        plot(mean_neu_stim(ch,:),'k');
        plot(mean_inc_stim(ch,:),'r');
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
        legend('con','neu','inc','stim','mean RT','Location','southeast');%,'inc-con'
    end
    
    % Plot RT locked
    subplot(n_chan,2,ch*2); hold on;
    if n_chan==1
        plot(mean_con_resp,'b');
        plot(mean_neu_resp,'k');
        plot(mean_inc_resp,'r');
    else
        plot(mean_con_resp(ch,:),'b');
        plot(mean_neu_resp(ch,:),'k');
        plot(mean_inc_resp(ch,:),'r');
    end
    %     plot(mean_diff_i_c_resp(ch,:),'g');
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
    ax.XTickLabel = x_lab_resp;
    title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
    if ch==1
        legend('con','neu','inc','resp','Location','southeast');%,'inc-con'
    end
end

% Save plots
out_filename = [fig_dir data_id '_ERPs_kla_bsln_' num2str(bsln_start) ...
    '.' num2str(bsln_end) '.' fig_type];
eval(['export_fig ' out_filename]);
% saveas(gcf,out_filename);

%% Plot Condition ERPs split by proportion congruency
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
    line([buff_pre, buff_pre],ylim);%, 'k--');
    line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
    ax = gca;
    ax.XLim = [0,buff_pre+max_RT+buff_post];
    ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
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
    line([buff_pre, buff_pre],ylim);%, 'k--');
    line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
    ax = gca;
    ax.XLim = [0,buff_pre+max_RT+buff_post];
    ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
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
    line([buff_pre, buff_pre],ylim);%, 'k--');
    line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
    ax = gca;
    ax.XLim = [0,buff_pre+max_RT+buff_post];
    ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
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
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
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
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
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
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
    ax.XTickLabel = x_lab_resp;
    title(strcat('Resp Locked, Inc: ',header_ecog.channel_labels(ch)));
%     legend('mcon','same','minc','avg','resp');%,'inc-con'
    
    out_filename = [fig_dir data_id '_ERPs_' num2str(ch) '_propcon_kla_bsln_' ...
        num2str(bsln_start) '.' num2str(bsln_end) '.' fig_type];
    eval(['export_fig ' out_filename]);
%  saveas(gcf,out_filename);
end
