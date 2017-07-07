%% Stroop analysis - frequency analyses
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
SBJ       = 'IR21';
data_id   = strcat(SBJ,'_RC_WM');
cond_lab  = {'con', 'neu', 'inc'}; % conditions to average
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [500, 1000];            % buffers are for cutting time series and then plotting
bsln_it   = 0;
bsln_lim  = [-200, 0];             % relative to stim
filt_it     = 1;                     % 0=no low pass; 1=pre-averaging lp; 2=post-averaging lp
filt_freq   = 10;                    % low pass frequency

% Plotting parameters
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
    lp_id = 'lp.none';
elseif filt_it == 1
    lp_id = ['lp' num2str(filt_freq) '.data'];
elseif filt_it == 2
    lp_id = ['lp' num2str(filt_freq) '.avgs'];
else
    error('low pass flag not in [0, 1, 2]');
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

%% Low pass if lp_it ==1
if filt_it == 1
    for ch_ix = 1:size(data_ecog,1)
        data_ecog(ch_ix,:) = eegfilt(data_ecog(ch_ix,:), header_ecog.sample_rate, [], filt_freq);
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
        if bsln_it == 1
            bsln_mean = nanmean(trial(buff_lim(1)+bsln_lim(1):buff_lim(1)+bsln_lim(2)));
            bsln_std  = nanstd(trial(buff_lim(1)+bsln_lim(1):buff_lim(1)+bsln_lim(2)));
            trial_norm = (trial-bsln_mean)/bsln_std;
        else
            trial_norm= trial;
        end
        
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
    if n_chan==1
        for c_ix = 1:length(cond_lab)
            eval(['plot(mean_stim.' cond_lab{c_ix} ',cond_colors{c_ix});']);
        end
    else
        for c_ix = 1:length(cond_lab)
            eval(['plot(mean_stim.' cond_lab{c_ix} '(ch,:),cond_colors{c_ix});']);
        end
    end
%     plot(mean_diff_i_c_stim(ch,:),'g');
    line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
    line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
    ax = gca;
    ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
    ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
    ax.XTickLabel = x_lab_stim;
    title(strcat('Stim Locked: ',header_ecog.channel_labels(ch)));
    if ch==1
        legend(cond_lab{:},'stim','mean RT','Location','southeast');%,'inc-con'
    end
    
    % Plot RT locked
    subplot(n_chan,2,ch*2); hold on;
    if n_chan==1
        for c_ix = 1:length(cond_lab)
            eval(['plot(mean_resp.' cond_lab{c_ix} ',cond_colors{c_ix});']);
        end
    else
        for c_ix = 1:length(cond_lab)
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
out_filename = [fig_dir data_id '_ERPs_ep' epoch_id '_bsln' bsln_id '_' lp_id '.' fig_type];
fprintf('Saving %s\n',out_filename);
eval(['export_fig ' out_filename]);
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
            '_propcon_bsln' bsln_id '_' lp_id '.' fig_type];
        fprintf('Saving %s\n',out_filename);
        eval(['export_fig ' out_filename]);
        %  saveas(gcf,out_filename);
    end
end