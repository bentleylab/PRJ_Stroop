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