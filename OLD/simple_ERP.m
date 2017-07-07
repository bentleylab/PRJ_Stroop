%% Stroop analysis - Simple ERPs, no baseline

SBJ = 'S02IR21';
ROI = 'RAC';
load(strcat('/home/knight/kla/Desktop/stroop/',SBJ,'/03_fulldata/',SBJ,'_',ROI,'.mat'));

n_chan = size(data_ecog,1);
n_trials = size(trial_info.condition_n,2);
% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
% trial_type = NaN(size(trial_info.condition_n));
max_RT = max(trial_info.response_time); % in sec

%% Cut trials
buff_pre = 1000;
buff_post = 1500;
trials_stim = NaN([n_chan,n_trials,buff_pre+buff_post+1]);
trials_resp = NaN([n_chan,n_trials,buff_pre+buff_post+1]);
for ch = 1:n_chan
    for t = 1:n_trials
        trials_stim(ch,t,:) = data_ecog(ch,trial_info.word_onset(t)-buff_pre:trial_info.word_onset(t)+buff_post);
        if ~isnan(trial_info.resp_onset(t))
            trials_resp(ch,t,:) = data_ecog(ch,trial_info.resp_onset(t)-buff_pre:trial_info.resp_onset(t)+buff_post);
        end
    end
end

%% ERPs by condition
con_idx = [trial_info.condition_n==1]+[trial_info.condition_n==2]+[trial_info.condition_n==3];
neu_idx = [trial_info.condition_n==4]+[trial_info.condition_n==5]+[trial_info.condition_n==6];
inc_idx = [trial_info.condition_n==7]+[trial_info.condition_n==8]+[trial_info.condition_n==9];

mean_con_stim = squeeze(nanmean(trials_stim(:,con_idx==1,:),2));
mean_neu_stim = squeeze(nanmean(trials_stim(:,neu_idx==1,:),2));
mean_inc_stim = squeeze(nanmean(trials_stim(:,inc_idx==1,:),2));
mean_diff_i_c_stim = mean_inc_stim-mean_con_stim;

mean_con_resp = squeeze(nanmean(trials_resp(:,con_idx==1,:),2));
mean_neu_resp = squeeze(nanmean(trials_resp(:,neu_idx==1,:),2));
mean_inc_resp = squeeze(nanmean(trials_resp(:,inc_idx==1,:),2));
mean_diff_i_c_resp = mean_inc_resp-mean_con_resp;

%% Plot ERPs
x_step = 250;
x_labels = -buff_pre:x_step:buff_post+1;

figure;
for ch = 1:n_chan
    subplot(n_chan,2,ch*2-1); hold on;
    plot(mean_con_stim(ch,:),'r');
    plot(mean_neu_stim(ch,:),'b');
    plot(mean_inc_stim(ch,:),'k');
%     plot(mean_diff_i_c_stim(ch,:),'g');
    line([1000, 1000],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
    ax.XTickLabel = x_labels;
    title(strcat('Stim Locked: ',header_ecog.channel_labels(ch)));
    legend('con','neu','inc','stim');%,'inc-con'
    
    subplot(n_chan,2,ch*2); hold on;
    plot(mean_con_resp(ch,:),'r');
    plot(mean_neu_resp(ch,:),'b');
    plot(mean_inc_resp(ch,:),'k');
%     plot(mean_diff_i_c_resp(ch,:),'g');
    line([1000, 1000],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,buff_pre+buff_post];
    ax.XTick = 0:x_step:buff_pre+buff_post;
    ax.XTickLabel = x_labels;
    title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
    legend('con','neu','inc','resp');%,'inc-con'
end