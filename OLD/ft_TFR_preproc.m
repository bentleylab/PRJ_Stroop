%% Stroop analysis - frequency analyses
clear all; close all
addpath('/home/knight/hoycw/Apps/fieldtrip-20160927');

SBJ = 'S02IR21';
ROI = 'RAC';

fig_type = 'eps';
fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/ERPs/',SBJ,'/');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

load(strcat('/home/knight/kla/Desktop/stroop/',SBJ,'/03_fulldata/',SBJ,'_',ROI,'.mat'));
if header_ecog.sample_rate~=1000
    error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!')
end

n_ch         = size(data_ecog,1);
RTs_all      = trial_info.response_time;
good_trl_idx = ~isnan(RTs_all);
n_trl        = sum(good_trl_idx);
RTs          = round(1000*RTs_all(good_trl_idx)); % converts sec to ms
max_RT       = max(RTs);
stim_onsets  = trial_info.word_onset(good_trl_idx);
trial_type   = trial_info.condition_n(good_trl_idx);


% Condition Indices (trial_info.condition_n):
% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
con_idx = trial_type==1 | trial_type==2 | trial_type==3;
neu_idx = trial_type==4 | trial_type==5 | trial_type==6;
inc_idx = trial_type==7 | trial_type==8 | trial_type==9;

mcon_idx = trial_type==3 | trial_type==6 | trial_type==9;
same_idx = trial_type==2 | trial_type==5 | trial_type==8;
minc_idx = trial_type==1 | trial_type==4 | trial_type==7;

%% Format trials
%CWH: option for resp locked trial is just create a trl based only on the
%   RTs, but the problem will be baselining...
%   stim locked trails can either be all the same length based on the max
%   RT or defined individually, which seems like more of a pain in the ass
%   (isn't ft supposed to take care of the paiN??? apparently not...)

% trl_stim.label = header_ecog.channel_labels';   % should be column
% trl_stim.fsample = header_ecog.sample_rate;
% % how do I insert the event labels? Which variable?

% Parameters
buff_pre   = 500;   %buffers are for cutting time series and then plotting
buff_post  = 500;
bsln_start = -200;  %relative to stim
bsln_end   = 0;     %relative to stim

cfg = [];
cfg.dataset = '/home/knight/hoycw/PRJ_Stroop/prcs_data/IR21/IR21_RAC.mat';
% cfg.trialfun = 'ft_trialfun_general';
% [trl, event] = ft_definetrial();
cfg.trl = [trial_info.word_onset(good_trl_idx), ...
    trial_info.resp_onset(good_trl_idx)+buff_post, repmat(buff_pre,sum(good_trl_idx),1)];
cfg.headerfile = '/home/knight/hoycw/PRJ_Stroop/prcs_data/IR21/IR21_RAC_hdr.mat';
cfg.padding = 1;
cfg.padtype = 'data';
cfg.continuous = 'yes';
cfg.channel = 'all';
% cfg.hpfreq = 0.8; %what kris already did, so hopefully nothing changes
% cfg.lpfreq = 200;
cfg.method = 'channel';

[data] = ft_preprocessing(cfg);
% Cut Trials
trl_stim.trial = {};%NaN([n_ch,n_trl,buff_pre+max_RT+buff_post+1]);
trl_resp.trial = {};%NaN([n_ch,n_trl,buff_pre+buff_post+1]);
for t = 1:n_trl
    if ~isnan(RTs(t)) %only analyze trials with responses
        trl = data_ecog(:,stim_onsets(t)-buff_pre:stim_onsets(t)+RTs(t)+buff_post);
        % Baseline the data
        bsln_mean = nanmean(trl(:,buff_pre+bsln_start:buff_pre+bsln_end),2);
        bsln_std  = nanstd(trl(:,buff_pre+bsln_start:buff_pre+bsln_end),0,2); %0 means normalize by n-1
        for ch = 1:n_ch
            trl_norm(ch,:) = (trl(ch,:)-bsln_mean(ch))/bsln_std(ch);
        end
        
        % Cut into trials
        trl_stim.trial{t} = trl_norm;
        trl_stim.time{t}  = 1:size(trl_norm,2);
        trl_resp.trial{t} = trl_norm(RTs(t):RTs(t)+buff_pre+buff_post);
        trl_resp.time{t}  = 1:size(trl_norm,2);
    end
    clear trl trl_norm
end

%% Frequency Analysis
% PSDs by condition
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'pow';
cfg.keeptrials = 'yes';
cfg.trials = con_idx;
cfg.foilim = [1 200];
cfg.tapsmofrq = 2;
cfg.toi = [buff_pre:size(;
% cfg.pad    = 100;

[psd_s] = ft_freqanalysis(cfg, trl_stim);
[psd_r] = ft_freqanalysis(cfg, trl_resp);

% % average TFR, locked to stim and resp
% func_envelopes_oneFband_resamp1K(PrcsDataDir,SUBJECT_TASK{SbjTask},FBAND);
% 
% % phase concentration, stim and resp
% 
% %% ERPs by condition
% 
% % Calculate stimulus locked ERP
% mean_con_stim = squeeze(nanmean(trl_stim(:,con_idx==1,:),2));
% mean_neu_stim = squeeze(nanmean(trl_stim(:,neu_idx==1,:),2));
% mean_inc_stim = squeeze(nanmean(trl_stim(:,inc_idx==1,:),2));
% mean_diff_i_c_stim = mean_inc_stim-mean_con_stim;
% 
% mean_con_mcon_stim = squeeze(nanmean(trl_stim(:,trial_type==3,:),2));
% mean_con_same_stim = squeeze(nanmean(trl_stim(:,trial_type==2,:),2));
% mean_con_minc_stim = squeeze(nanmean(trl_stim(:,trial_type==1,:),2));
% mean_neu_mcon_stim = squeeze(nanmean(trl_stim(:,trial_type==6,:),2));
% mean_neu_same_stim = squeeze(nanmean(trl_stim(:,trial_type==5,:),2));
% mean_neu_minc_stim = squeeze(nanmean(trl_stim(:,trial_type==4,:),2));
% mean_inc_mcon_stim = squeeze(nanmean(trl_stim(:,trial_type==9,:),2));
% mean_inc_same_stim = squeeze(nanmean(trl_stim(:,trial_type==8,:),2));
% mean_inc_minc_stim = squeeze(nanmean(trl_stim(:,trial_type==7,:),2));
% % mean_diff_mi_mc_stim = mean_minc_stim-mean_mcon_stim;
% 
% % Calculate RT locked ERP
% mean_con_resp = squeeze(nanmean(trl_resp(:,con_idx==1,:),2));
% mean_neu_resp = squeeze(nanmean(trl_resp(:,neu_idx==1,:),2));
% mean_inc_resp = squeeze(nanmean(trl_resp(:,inc_idx==1,:),2));
% mean_diff_i_c_resp = mean_inc_resp-mean_con_resp;
% 
% mean_con_mcon_resp = squeeze(nanmean(trl_resp(:,trial_type==3,:),2));
% mean_con_same_resp = squeeze(nanmean(trl_resp(:,trial_type==2,:),2));
% mean_con_minc_resp = squeeze(nanmean(trl_resp(:,trial_type==1,:),2));
% mean_neu_mcon_resp = squeeze(nanmean(trl_resp(:,trial_type==6,:),2));
% mean_neu_same_resp = squeeze(nanmean(trl_resp(:,trial_type==5,:),2));
% mean_neu_minc_resp = squeeze(nanmean(trl_resp(:,trial_type==4,:),2));
% mean_inc_mcon_resp = squeeze(nanmean(trl_resp(:,trial_type==9,:),2));
% mean_inc_same_resp = squeeze(nanmean(trl_resp(:,trial_type==8,:),2));
% mean_inc_minc_resp = squeeze(nanmean(trl_resp(:,trial_type==7,:),2));
% % mean_mcon_resp = squeeze(nanmean(trials_resp(:,mcon_idx==1,:),2));
% % mean_same_resp = squeeze(nanmean(trials_resp(:,same_idx==1,:),2));
% % mean_minc_resp = squeeze(nanmean(trials_resp(:,minc_idx==1,:),2));
% % mean_diff_mi_mc_resp = mean_minc_resp-mean_mcon_resp;
% 
% %% Plot ERPs
% x_step = 250;
% x_lab_stim = -buff_pre:x_step:max_RT+buff_post;
% x_lab_resp = -buff_pre:x_step:buff_post;
% 
% % Plot Condition ERPs
% figure;
% if n_ch > 2
%     fig_height = 1;
% else
%     fig_height = n_ch/3;
% end
% set(gcf,'units','normalized','outerposition',[0 0 1 fig_height]);
% for ch = 1:n_ch
%     % Plot Stim Locked
%     subplot(n_ch,2,ch*2-1); hold on;
%     if n_ch==1
%         plot(mean_con_stim,'b');
%         plot(mean_neu_stim,'k');
%         plot(mean_inc_stim,'r');
%     else
%         plot(mean_con_stim(ch,:),'b');
%         plot(mean_neu_stim(ch,:),'k');
%         plot(mean_inc_stim(ch,:),'r');
%     end
% %     plot(mean_diff_i_c_stim(ch,:),'g');
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+max_RT+buff_post];
%     ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
%     ax.XTickLabel = x_lab_stim;
%     title(strcat('Stim Locked: ',header_ecog.channel_labels(ch)));
%     if ch==1
%         legend('con','neu','inc','stim','mean RT','Location','southeast');%,'inc-con'
%     end
%     
%     % Plot RT locked
%     subplot(n_ch,2,ch*2); hold on;
%     if n_ch==1
%         plot(mean_con_resp,'b');
%         plot(mean_neu_resp,'k');
%         plot(mean_inc_resp,'r');
%     else
%         plot(mean_con_resp(ch,:),'b');
%         plot(mean_neu_resp(ch,:),'k');
%         plot(mean_inc_resp(ch,:),'r');
%     end
%     %     plot(mean_diff_i_c_resp(ch,:),'g');
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+buff_post];
%     ax.XTick = 0:x_step:buff_pre+buff_post;
%     ax.XTickLabel = x_lab_resp;
%     title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
%     if ch==1
%         legend('con','neu','inc','resp','Location','southeast');%,'inc-con'
%     end
% end
% 
% out_filename = [fig_dir SBJ '_ERPs_' ROI '_kla_bsln_' num2str(bsln_start) ...
%     '.' num2str(bsln_end) '.' fig_type];
% saveas(gcf,out_filename);
% 
% %% Plot Condition ERPs split by proportion congruency
% for ch = 1:n_ch
%     figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%     % Plot Stim Locked
%     % Congruent
%     subplot(3,2,1); hold on;
%     if n_ch==1
%         plot(mean_con_mcon_stim,'r');
%         plot(mean_con_same_stim,'g');
%         plot(mean_con_minc_stim,'b');
%         plot(mean_con_stim,'k','LineWidth',3);
%     else
%         plot(mean_con_mcon_stim(ch,:),'r');
%         plot(mean_con_same_stim(ch,:),'g');
%         plot(mean_con_minc_stim(ch,:),'b');
%         plot(mean_con_stim(ch,:),'k','LineWidth',3);
%     end
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+max_RT+buff_post];
%     ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
%     ax.XTickLabel = x_lab_stim;
%     title(strcat('Stim Locked, Con: ',header_ecog.channel_labels(ch)));
%     legend('mcon','same','minc','avg','stim','mean RT','Location','southeast');%,'inc-con'
% 
%     % Neutral
%     subplot(3,2,3); hold on;
%     if n_ch==1
%         plot(mean_neu_mcon_stim,'r');
%         plot(mean_neu_same_stim,'g');
%         plot(mean_neu_minc_stim,'b');
%         plot(mean_neu_stim,'k','LineWidth',3);
%     else
%         plot(mean_neu_mcon_stim(ch,:),'r');
%         plot(mean_neu_same_stim(ch,:),'g');
%         plot(mean_neu_minc_stim(ch,:),'b');
%         plot(mean_neu_stim(ch,:),'k','LineWidth',3);
%     end
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+max_RT+buff_post];
%     ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
%     ax.XTickLabel = x_lab_stim;
%     title(strcat('Stim Locked, Neu: ',header_ecog.channel_labels(ch)));
% %     legend('mcon','same','minc','avg','stim','mean RT');%,'inc-con'
% 
%     % Incongruent
%     subplot(3,2,5); hold on;
%     if n_ch==1
%         plot(mean_inc_mcon_stim,'r');
%         plot(mean_inc_same_stim,'g');
%         plot(mean_inc_minc_stim,'b');
%         plot(mean_inc_stim,'k','LineWidth',3);
%     else
%         plot(mean_inc_mcon_stim(ch,:),'r');
%         plot(mean_inc_same_stim(ch,:),'g');
%         plot(mean_inc_minc_stim(ch,:),'b');
%         plot(mean_inc_stim(ch,:),'k','LineWidth',3);
%     end
%     %     plot(mean_diff_i_c_stim(ch,:),'g');
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     line([buff_pre+nanmean(RTs), buff_pre+nanmean(RTs)],ylim,'LineStyle','--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+max_RT+buff_post];
%     ax.XTick = 0:x_step:buff_pre+max_RT+buff_post;
%     ax.XTickLabel = x_lab_stim;
%     title(strcat('Stim Locked, Inc: ',header_ecog.channel_labels(ch)));
% %     legend('mcon','same','minc','avg','stim','mean RT');%,'inc-con'
%     
%     % Plot RT locked
%     % Congruent
%     subplot(3,2,2); hold on;
%     if n_ch==1
%         plot(mean_con_mcon_resp,'r');
%         plot(mean_con_same_resp,'g');
%         plot(mean_con_minc_resp,'b');
%         plot(mean_con_resp,'k','LineWidth',3);
%     else
%         plot(mean_con_mcon_resp(ch,:),'r');
%         plot(mean_con_same_resp(ch,:),'g');
%         plot(mean_con_minc_resp(ch,:),'b');
%         plot(mean_con_resp(ch,:),'k','LineWidth',3);
%     end
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+buff_post];
%     ax.XTick = 0:x_step:buff_pre+buff_post;
%     ax.XTickLabel = x_lab_resp;
%     title(strcat('Resp Locked, Con: ',header_ecog.channel_labels(ch)));
%     legend('mcon','same','minc','avg','resp','Location','southeast');%,'inc-con'
%     
%     % Neutral
%     subplot(3,2,4); hold on;
%     if n_ch==1
%         plot(mean_neu_mcon_resp,'r');
%         plot(mean_neu_same_resp,'g');
%         plot(mean_neu_minc_resp,'b');
%         plot(mean_neu_resp,'k','LineWidth',3);
%     else
%         plot(mean_neu_mcon_resp(ch,:),'r');
%         plot(mean_neu_same_resp(ch,:),'g');
%         plot(mean_neu_minc_resp(ch,:),'b');
%         plot(mean_neu_resp(ch,:),'k','LineWidth',3);
%     end
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+buff_post];
%     ax.XTick = 0:x_step:buff_pre+buff_post;
%     ax.XTickLabel = x_lab_resp;
%     title(strcat('Resp Locked, Neu: ',header_ecog.channel_labels(ch)));
% %     legend('mcon','same','minc','avg','resp');%,'inc-con'
% 
%     % Incongruent
%     subplot(3,2,6); hold on;
%     if n_ch==1
%         plot(mean_inc_mcon_resp,'r');
%         plot(mean_inc_same_resp,'g');
%         plot(mean_inc_minc_resp,'b');
%         plot(mean_inc_resp,'k','LineWidth',3);
%     else
%         plot(mean_inc_mcon_resp(ch,:),'r');
%         plot(mean_inc_same_resp(ch,:),'g');
%         plot(mean_inc_minc_resp(ch,:),'b');
%         plot(mean_inc_resp(ch,:),'k','LineWidth',3);
%     end
%     line([buff_pre, buff_pre],ylim);%, 'k--');
%     ax = gca;
%     ax.XLim = [0,buff_pre+buff_post];
%     ax.XTick = 0:x_step:buff_pre+buff_post;
%     ax.XTickLabel = x_lab_resp;
%     title(strcat('Resp Locked, Inc: ',header_ecog.channel_labels(ch)));
% %     legend('mcon','same','minc','avg','resp');%,'inc-con'
%     
%     out_filename = [fig_dir SBJ '_ERPs_' ROI num2str(ch) '_propcon_kla_bsln_' ...
%         num2str(bsln_start) '.' num2str(bsln_end) '.' fig_type];
%     saveas(gcf,out_filename);
% end
