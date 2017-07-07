%% Stroop analysis - High Gamma Power
!!! make this script a single plot for S R with a gap between
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
SBJ       = 'IR39';
data_id   = strcat(SBJ,'_RAC_WM');
cond_lab  = {'con', 'neu', 'inc'}; % conditions to average
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [500, 500];            % buffers are for cutting time series and then plotting

% Analysis parameters
hg_type      = 'multiband';            % 'wideband', 'multiband'
bsln_type    = 'demean';               % 'zscore', 'demean', 'none'
bsln_lim     = [200, 0];               % ms before and after event
env_it       = 1;                      % take the envelope [0 1]
smooth_it    = 1;                      % smooth the result before averaging (0/1)
smooth_freq  = 10;
trial_type   = 'nanpad';               % 'nanpad', 'datapad'
avg_channels = 1;                      % average across channels

% Plotting parameters
save_fig      = 1;
plot_lim_s    = [250 1000];           % time window to plot around stimulus
plot_lim_r    = [500 500];            % time window to plot around response
x_step        = 250;                  % step of x tick marks
plot_prop_con = 0;
cond_colors   = {'c', 'k', 'r'};    % colors for cond_lab plotting
block_colors  = {'r', 'g', 'b'};    % colors for [mcon, same, minc]
prop_con_lab  = {'con-mcon', 'con-same', 'con-minc', 'neu-mcon', 'neu-same',...
    'neu-minc', 'inc-mcon', 'inc-same', 'inc-minc'};
fig_type      = 'eps';

%% Process parameters
% High Gamma Extraction
if strcmp(hg_type,'wideband')
    hg_cf = [110];
    hg_bw = 80;
elseif strcmp(hg_type,'multiband')
    hg_cf = [70:10:150];
    hg_bw = 10;
else
    error(strcat('Unknown high gamma extraction type: ',hg_type));
end
for f_ix = 1:length(hg_cf)
    hg_lim{f_ix} = [hg_cf(f_ix)-hg_bw/2 hg_cf(f_ix)+hg_bw/2];
end

% Filtering
if env_it ==1
    env_id = 'env_';
else
    env_id = '';
end
if smooth_it == 0
    smooth_id = '';
elseif smooth_it == 1
    smooth_id = strcat('sm',num2str(smooth_freq));
% elseif smooth_it == 2
%     smooth_id = ['_sm.avgs.' num2str(smooth_freq)];
else
    error('smooth_it not in [0,1]');
end

% Epoching/Averaging
epoch_id = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];
bsln_id  = [num2str(bsln_lim(1)) '.' num2str(bsln_lim(2))];
if strcmp(trial_type,'nanpad')
    trl_id = 'nan';
elseif strcmp(trial_type,'datapad')
    trl_id = 'data';
else
    error(strcat('Unknown trial type: ',trial_type));
end
if avg_channels == 1
    grand_id = '_ch.avg';
else
    grand_id = '';
end

SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/power/',SBJ,'/');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Load data
load(strcat(SBJ_dir,'04_proc/',data_id,'.mat'));
% need only ok_epochs index, the rest is all in 04_proc file (variables seem same as of 10/10/16...)
load(strcat(SBJ_dir,'06_epochs/',data_id,'_',epoch_id,'.mat'),'ok_epochs');

% Correct all info to only good trials
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

%% High Gamma Extraction
% Filter data
data_filt = NaN([size(data_ecog,1) length(hg_cf) size(data_ecog,2)]);
for ch_ix = 1:size(data_ecog,1)
    for cf_ix = 1:length(hg_cf)
        % Filter data
        data_filt(ch_ix,cf_ix,:) = eegfilt(data_ecog(ch_ix,:), ...
            header_ecog.sample_rate, hg_lim{cf_ix}(1), []);
        data_filt(ch_ix,cf_ix,:) = eegfilt(data_filt(ch_ix,cf_ix,:), ...
            header_ecog.sample_rate, [], hg_lim{cf_ix}(2));
        % Envelope
        if env_it ==1
            data_filt(ch_ix,cf_ix,:) = abs(hilbert(data_filt(ch_ix,cf_ix,:)));
        end
        % Baseline?
    end
end
% Average to get single HG time series
data_hg = mean(data_filt,2);

% Smooth
if smooth_it == 1
    for ch_ix = 1:size(data_ecog,1)
        data_hg(ch_ix,:) = eegfilt(data_hg(ch_ix,:), header_ecog.sample_rate, [], smooth_freq);
    end
end
   
%% Format trials
% Cut Trials
trials_stim = NaN([size(data_hg,1),n_trials,buff_lim(1)+max_RT+buff_lim(2)+1]);
trials_resp = NaN([size(data_hg,1),n_trials,buff_lim(1)+buff_lim(2)+1]);
for ch = 1:size(data_hg,1)
    % Cut trials
    if strcmp(trial_type,'nanpad')
        trials_stim(ch,:,:) = fn_epoch_cuts_nanpad(data_hg(ch,:),word_onset,resp_onset,buff_lim);
        trials_resp(ch,:,:) = fn_epoch_cuts_nanpad(data_hg(ch,:),resp_onset,resp_onset,buff_lim);
    elseif strcmp(trial_type,'datapad')
        trials_stim(ch,:,:) = fn_epoch_cuts_datapad(data_hg(ch,:),word_onset,resp_onset,buff_lim);
        trials_resp(ch,:,:) = fn_epoch_cuts_datapad(data_hg(ch,:),resp_onset,resp_onset,buff_lim);
    else
        error('Bad trial_type');
    end
    
    % Baseline Trials
    trials_stim(ch,:,:) = fn_bsln_epochs(squeeze(trials_stim(ch,:,:)),...
        repmat(buff_lim(1),n_trials),bsln_lim,bsln_type);
    trials_resp(ch,:,:) = fn_bsln_epochs(squeeze(trials_resp(ch,:,:)),...
        repmat(buff_lim(1),n_trials),bsln_lim,bsln_type);
end

%% ERPs by condition
for ix = 1:length(cond_lab)
    % Get binary condition index
    eval(['cond_idx.' cond_lab{ix} ' = fn_condition_index(cond_lab{ix}, cond_n);']);
    % if plot_prop_con, get that label
    
    % Average trials matching condition index
    eval(['mean_stim.' cond_lab{ix} ' = squeeze(nanmean(trials_stim(:,cond_idx.' cond_lab{ix} '==1,:),2));']);
    eval(['mean_resp.' cond_lab{ix} ' = squeeze(nanmean(trials_resp(:,cond_idx.' cond_lab{ix} '==1,:),2));']);
    
    % Average channels if required
    if avg_channels == 1
        eval(['mean_stim.' cond_lab{ix} ' = squeeze(nanmean(mean_stim.' cond_lab{ix} ',1));']);
        eval(['mean_resp.' cond_lab{ix} ' = squeeze(nanmean(mean_resp.' cond_lab{ix} ',1));']);
    end
    
    % Low pass if required
    if smooth_it == 2
        eval(['mean_stim.' cond_lab{ix} ' = eegfilt(mean_stim.' cond_lab{ix} ...
            ', header_ecog.sample_rate, [], smooth_freq);']);
        eval(['mean_resp.' cond_lab{ix} ' = eegfilt(mean_resp.' cond_lab{ix} ...
            ', header_ecog.sample_rate, [], smooth_freq);']);
    end

end

% mean_diff_i_c_stim = mean_inc_stim-mean_con_stim;
% % mean_diff_mi_mc_stim = mean_minc_stim-mean_mcon_stim;
% 
% mean_diff_i_c_resp = mean_inc_resp-mean_con_resp;
% % mean_diff_mi_mc_resp = mean_minc_resp-mean_mcon_resp;

%% Plot ERPs
x_lab_stim   = -plot_lim_s(1):x_step:plot_lim_s(2);
x_lab_resp   = -plot_lim_r(1):x_step:plot_lim_r(2);
stim_win_on  = buff_lim(1)-plot_lim_s(1)+1;
resp_win_on  = buff_lim(1)-plot_lim_r(1)+1;
stim_win_off = buff_lim(1)+plot_lim_s(2)+1;
resp_win_off = buff_lim(1)+plot_lim_r(2)+1;

% Plot Condition ERPs
out_file = [data_id '_HG' grand_id '_ep' epoch_id '_' trl_id '_bsln' bsln_id '_' env_id smooth_id];
figure;suptitle(strrep(out_file,'_','\_'));
n_plots = eval(['size(mean_stim.' cond_lab{1} ',1)']);
if n_plots > 2
    fig_height = 1;
else
    fig_height = n_plots/3;
end
set(gcf,'units','normalized','outerposition',[0 0 1 fig_height]);
for ch = 1:n_plots
    % Plot Stim Locked
    subplot(n_plots,2,ch*2-1); hold on;
    for c_ix = 1:length(cond_lab)
        if n_plots==1
            eval(['plot(mean_stim.' cond_lab{c_ix} '(stim_onset:stim_offset),cond_colors{c_ix});']);
        else
            eval(['plot(mean_stim.' cond_lab{c_ix} '(ch,stim_onset:stim_offset),cond_colors{c_ix});']);
        end
    end
%     plot(mean_diff_i_c_stim(ch,:),'g');
    line([plot_lim_s(1) plot_lim_s(1)],ylim);%, 'k--');
    line([plot_lim_s(1)+nanmean(RTs), plot_lim_s(1)+nanmean(RTs)],ylim,'LineStyle','--');
    % Axes and Labels
    ax = gca;
    ax.XLim = [0,sum(plot_lim_s)];
    ax.XTick = 0:x_step:sum(plot_lim_s);
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
    subplot(n_plots,2,ch*2); hold on;
    for c_ix = 1:length(cond_lab)
        if n_plots==1
            eval(['plot(mean_resp.' cond_lab{c_ix} '(resp_onset:resp_offset,cond_colors{c_ix});']);
        else
            eval(['plot(mean_resp.' cond_lab{c_ix} '(ch,resp_onset:resp_offset),cond_colors{c_ix});']);
        end
    end
    %     plot(mean_diff_i_c_resp(ch,:),'g');
    line([plot_lim_r(1) plot_lim_r(1)],ylim);%, 'k--');
    ax = gca;
    ax.XLim = [0,sum(plot_lim_r)];
    ax.XTick = 0:x_step:sum(plot_lim_r);
    ax.XTickLabel = x_lab_resp;
    title(strcat('Resp Locked: ',header_ecog.channel_labels(ch)));
    if ch==1
        legend(cond_lab{:},'resp','Location','southeast');%,'inc-con'
    end
end

% Save plots
out_filename = [fig_dir out_file '.' fig_type];
if save_fig ==1
    fprintf('Saving %s\n',out_filename);
    eval(['export_fig ' out_filename]);
end
% saveas(gcf,out_filename);

%% Plot Condition ERPs split by proportion congruency
% % !!!
% if plot_prop_con == 1
%     for ch = 1:n_chan
%         figure;set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%         % Plot Stim Locked
%         % Congruent
%         subplot(3,2,1); hold on;
%         if n_chan==1
%             plot(mean_con_mcon_stim,'r');
%             plot(mean_con_same_stim,'g');
%             plot(mean_con_minc_stim,'b');
%             plot(mean_con_stim,'k','LineWidth',3);
%         else
%             plot(mean_con_mcon_stim(ch,:),'r');
%             plot(mean_con_same_stim(ch,:),'g');
%             plot(mean_con_minc_stim(ch,:),'b');
%             plot(mean_con_stim(ch,:),'k','LineWidth',3);
%         end
%         line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
%         line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
%         ax = gca;
%         ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
%         ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
%         ax.XTickLabel = x_lab_stim;
%         title(strcat('Stim Locked, Con: ',header_ecog.channel_labels(ch)));
%         legend('mcon','same','minc','avg','stim','mean RT','Location','southeast');%,'inc-con'
%         
%         % Neutral
%         subplot(3,2,3); hold on;
%         if n_chan==1
%             plot(mean_neu_mcon_stim,'r');
%             plot(mean_neu_same_stim,'g');
%             plot(mean_neu_minc_stim,'b');
%             plot(mean_neu_stim,'k','LineWidth',3);
%         else
%             plot(mean_neu_mcon_stim(ch,:),'r');
%             plot(mean_neu_same_stim(ch,:),'g');
%             plot(mean_neu_minc_stim(ch,:),'b');
%             plot(mean_neu_stim(ch,:),'k','LineWidth',3);
%         end
%         line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
%         line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
%         ax = gca;
%         ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
%         ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
%         ax.XTickLabel = x_lab_stim;
%         title(strcat('Stim Locked, Neu: ',header_ecog.channel_labels(ch)));
%         %     legend('mcon','same','minc','avg','stim','mean RT');%,'inc-con'
%         
%         % Incongruent
%         subplot(3,2,5); hold on;
%         if n_chan==1
%             plot(mean_inc_mcon_stim,'r');
%             plot(mean_inc_same_stim,'g');
%             plot(mean_inc_minc_stim,'b');
%             plot(mean_inc_stim,'k','LineWidth',3);
%         else
%             plot(mean_inc_mcon_stim(ch,:),'r');
%             plot(mean_inc_same_stim(ch,:),'g');
%             plot(mean_inc_minc_stim(ch,:),'b');
%             plot(mean_inc_stim(ch,:),'k','LineWidth',3);
%         end
%         %     plot(mean_diff_i_c_stim(ch,:),'g');
%         line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
%         line([buff_lim(1)+nanmean(RTs), buff_lim(1)+nanmean(RTs)],ylim,'LineStyle','--');
%         ax = gca;
%         ax.XLim = [0,buff_lim(1)+max_RT+buff_lim(2)];
%         ax.XTick = 0:x_step:buff_lim(1)+max_RT+buff_lim(2);
%         ax.XTickLabel = x_lab_stim;
%         title(strcat('Stim Locked, Inc: ',header_ecog.channel_labels(ch)));
%         %     legend('mcon','same','minc','avg','stim','mean RT');%,'inc-con'
%         
%         % Plot RT locked
%         % Congruent
%         subplot(3,2,2); hold on;
%         if n_chan==1
%             plot(mean_con_mcon_resp,'r');
%             plot(mean_con_same_resp,'g');
%             plot(mean_con_minc_resp,'b');
%             plot(mean_con_resp,'k','LineWidth',3);
%         else
%             plot(mean_con_mcon_resp(ch,:),'r');
%             plot(mean_con_same_resp(ch,:),'g');
%             plot(mean_con_minc_resp(ch,:),'b');
%             plot(mean_con_resp(ch,:),'k','LineWidth',3);
%         end
%         line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
%         ax = gca;
%         ax.XLim = [0,buff_lim(1)+buff_lim(2)];
%         ax.XTick = 0:x_step:buff_lim(1)+buff_lim(2);
%         ax.XTickLabel = x_lab_resp;
%         title(strcat('Resp Locked, Con: ',header_ecog.channel_labels(ch)));
%         legend('mcon','same','minc','avg','resp','Location','southeast');%,'inc-con'
%         
%         % Neutral
%         subplot(3,2,4); hold on;
%         if n_chan==1
%             plot(mean_neu_mcon_resp,'r');
%             plot(mean_neu_same_resp,'g');
%             plot(mean_neu_minc_resp,'b');
%             plot(mean_neu_resp,'k','LineWidth',3);
%         else
%             plot(mean_neu_mcon_resp(ch,:),'r');
%             plot(mean_neu_same_resp(ch,:),'g');
%             plot(mean_neu_minc_resp(ch,:),'b');
%             plot(mean_neu_resp(ch,:),'k','LineWidth',3);
%         end
%         line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
%         ax = gca;
%         ax.XLim = [0,buff_lim(1)+buff_lim(2)];
%         ax.XTick = 0:x_step:buff_lim(1)+buff_lim(2);
%         ax.XTickLabel = x_lab_resp;
%         title(strcat('Resp Locked, Neu: ',header_ecog.channel_labels(ch)));
%         %     legend('mcon','same','minc','avg','resp');%,'inc-con'
%         
%         % Incongruent
%         subplot(3,2,6); hold on;
%         if n_chan==1
%             plot(mean_inc_mcon_resp,'r');
%             plot(mean_inc_same_resp,'g');
%             plot(mean_inc_minc_resp,'b');
%             plot(mean_inc_resp,'k','LineWidth',3);
%         else
%             plot(mean_inc_mcon_resp(ch,:),'r');
%             plot(mean_inc_same_resp(ch,:),'g');
%             plot(mean_inc_minc_resp(ch,:),'b');
%             plot(mean_inc_resp(ch,:),'k','LineWidth',3);
%         end
%         line([buff_lim(1), buff_lim(1)],ylim);%, 'k--');
%         ax = gca;
%         ax.XLim = [0,buff_lim(1)+buff_lim(2)];
%         ax.XTick = 0:x_step:buff_lim(1)+buff_lim(2);
%         ax.XTickLabel = x_lab_resp;
%         title(strcat('Resp Locked, Inc: ',header_ecog.channel_labels(ch)));
%         %     legend('mcon','same','minc','avg','resp');%,'inc-con'
%         
%         out_filename = [fig_dir data_id '_ERPs_' num2str(ch) '_ep' epoch_id ...
%             '_propcon_bsln' bsln_id '_' bp_id '.' fig_type];
%         fprintf('Saving %s\n',out_filename);
%         eval(['export_fig ' out_filename]);
%         %  saveas(gcf,out_filename);
%     end
% end