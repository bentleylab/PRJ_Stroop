%% Stroop analysis - Simple ERPs, no baseline

SBJ = 'S02IR21';
ROI = 'RAC';
load(strcat('/home/knight/kla/Desktop/stroop/',SBJ,'/03_fulldata/',SBJ,'_',ROI,'.mat'));

n_chan = size(data_ecog,1);
n_trials = size(trial_info.condition_n,2);
if header_ecog.sample_rate~=1000
    error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!')
end
% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
% trial_type = NaN(size(trial_info.condition_n));
RTs = round(1000*trial_info.response_time); % converts sec to ms
max_RT = max(RTs);

%% Format trials
% Parameters
buff_pre   = 500;   %buffers are for cutting time series and then plotting
buff_post  = 500;
bsln_start = -100;  %relative to stim
bsln_end   = 0;     %relative to stim

% Cut Trials
trials_stim = NaN([n_chan,n_trials,buff_pre+max_RT+buff_post+1]);
trials_resp = NaN([n_chan,n_trials,buff_pre+buff_post+1]);
for ch = 1:n_chan
    for t = 1:n_trials
        trial = data_ecog(ch,trial_info.word_onset(t)-buff_pre:trial_info.word_onset(t)+max_RT+buff_post);
        % Baseline the data
        bsln_mean = nanmean(trial(buff_pre+bsln_start:buff_pre+bsln_end));
        bsln_std  = nanstd(trial(buff_pre+bsln_start:buff_pre+bsln_end));
        trial_norm = (trial-bsln_mean)/bsln_std;
        
        % Cut into trials
        trials_stim(ch,t,:) = trial_norm;
        if ~isnan(trial_info.resp_onset(t))
            trials_resp(ch,t,:) = trial_norm(RTs(t):RTs(t)+buff_pre+buff_post);
        end
    end
    
    
end

%% ERPs by condition
con_idx = [trial_info.condition_n==1]+[trial_info.condition_n==2]+[trial_info.condition_n==3];
neu_idx = [trial_info.condition_n==4]+[trial_info.condition_n==5]+[trial_info.condition_n==6];
inc_idx = [trial_info.condition_n==7]+[trial_info.condition_n==8]+[trial_info.condition_n==9];

mcon_idx = [trial_info.condition_n==3]+[trial_info.condition_n==6]+[trial_info.condition_n==9];
same_idx = [trial_info.condition_n==2]+[trial_info.condition_n==5]+[trial_info.condition_n==8];
minc_idx = [trial_info.condition_n==1]+[trial_info.condition_n==4]+[trial_info.condition_n==7];

% Calculate stimulus locked ERP
mean_con_stim = squeeze(nanmean(trials_stim(:,con_idx==1,:),2));
mean_neu_stim = squeeze(nanmean(trials_stim(:,neu_idx==1,:),2));
mean_inc_stim = squeeze(nanmean(trials_stim(:,inc_idx==1,:),2));
mean_diff_i_c_stim = mean_inc_stim-mean_con_stim;

mean_con_mcon_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==3,:),2));
mean_con_same_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==2,:),2));
mean_con_minc_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==1,:),2));
mean_neu_mcon_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==6,:),2));
mean_neu_same_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==5,:),2));
mean_neu_minc_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==4,:),2));
mean_inc_mcon_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==9,:),2));
mean_inc_same_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==8,:),2));
mean_inc_minc_stim = squeeze(nanmean(trials_stim(:,trial_info.condition_n==7,:),2));
% mean_diff_mi_mc_stim = mean_minc_stim-mean_mcon_stim;

% Calculate RT locked ERP
mean_con_resp = squeeze(nanmean(trials_resp(:,con_idx==1,:),2));
mean_neu_resp = squeeze(nanmean(trials_resp(:,neu_idx==1,:),2));
mean_inc_resp = squeeze(nanmean(trials_resp(:,inc_idx==1,:),2));
mean_diff_i_c_resp = mean_inc_resp-mean_con_resp;

mean_con_mcon_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==3,:),2));
mean_con_same_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==2,:),2));
mean_con_minc_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==1,:),2));
mean_neu_mcon_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==6,:),2));
mean_neu_same_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==5,:),2));
mean_neu_minc_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==4,:),2));
mean_inc_mcon_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==9,:),2));
mean_inc_same_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==8,:),2));
mean_inc_minc_resp = squeeze(nanmean(trials_resp(:,trial_info.condition_n==7,:),2));
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
for ch = 1:n_chan
    % Plot Stim Locked
    subplot(n_chan,2,ch*2-1); hold on;
    plot(mean_con_stim(ch,:),'r');
    plot(mean_neu_stim(ch,:),'b');
    plot(mean_inc_stim(ch,:),'k');
%     plot(mean_diff_i_c_stim(ch,:),'g');
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+max_RT+buff_post];
    ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
    ax.XTickLabel = x_lab_stim;
    title(strcat('Stim Locked: ',header_ecog.channel_labels(ch)));
    legend('con','neu','inc','stim');%,'inc-con'
    
    % Plot RT locked
    subplot(n_chan,2,ch*2); hold on;
    plot(mean_con_resp(ch,:),'r');
    plot(mean_neu_resp(ch,:),'b');
    plot(mean_inc_resp(ch,:),'k');
%     plot(mean_diff_i_c_resp(ch,:),'g');
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
    ax.XTickLabel = x_lab_resp;
    title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
    legend('con','neu','inc','resp');%,'inc-con'
end

% Plot Condition ERPs split by proportion congruency
figure;
for ch = 1:n_chan
    % Plot Stim Locked
    subplot(n_chan,2,ch*2-1); hold on;
    plot(mean_con_stim(ch,:),'r');
    plot(mean_neu_stim(ch,:),'b');
    plot(mean_inc_stim(ch,:),'k');
%     plot(mean_diff_i_c_stim(ch,:),'g');
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+max_RT+buff_post];
    ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
    ax.XTickLabel = x_lab_stim;
    title(strcat('Stim Locked: ',header_ecog.channel_labels(ch)));
    legend('con','neu','inc','stim');%,'inc-con'
    
    % Plot RT locked
    subplot(n_chan,2,ch*2); hold on;
    plot(mean_con_resp(ch,:),'r');
    plot(mean_neu_resp(ch,:),'b');
    plot(mean_inc_resp(ch,:),'k');
%     plot(mean_diff_i_c_resp(ch,:),'g');
    line([buff_pre, buff_pre],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
    ax.XTickLabel = x_lab_resp;
    title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
    legend('con','neu','inc','resp');%,'inc-con'
end
