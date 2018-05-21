%% RT Behavioral analysis
clear all; close all
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));

% Analysis Parameters
SBJ       = 'IR35';
data_id   = strcat(SBJ,'_LAC_WM');
late_RT_cut = 600;      %ms window before next stim to eliminate

% Plotting parameters
save_fig      = 1;
vis           = 'on';
n_bins        = 50;
line_w        = 2;
trial_lab     = {'con', 'neu', 'inc'};
block_lab     = {'mcon', 'same', 'minc'};
trial_colors   = {'b', 'k', 'r'};    % colors for cond_lab plotting
block_colors  = {repmat(0.8,3,1), repmat(0.5,3,1), repmat(0.2,3,1)};    % colors for [mcon, same, minc]
prop_con_lab  = {'con_mcon', 'con_same', 'con_minc', 'neu_mcon', 'neu_same',...
    'neu_minc', 'inc_mcon', 'inc_same', 'inc_minc'};
fig_type      = 'eps';

% Process parameters
SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/RTs/',SBJ,'/');
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Load data
load(strcat(SBJ_dir,'06_epochs/',data_id,'_clean.mat'),'trial_info');

% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
% trial_type = NaN(size(trial_info.condition_n));
RTs       = round(1000*trial_info.response_time); % converts sec to ms
max_RT    = max(RTs);
n_edges = linspace(min(RTs),max_RT,n_bins);

% Check for RTs overlapping with stim onset or baseline
late_RT_idx = [];
for t_ix = 2:length(trial_info.resp_onset)
    if trial_info.word_onset(t_ix)-late_RT_cut*trial_info.sample_rate/1000 <= trial_info.resp_onset(t_ix-1)
        late_RT_idx = [late_RT_idx(:) t_ix];
    end
end
fprintf('%i late trials detected\n',length(late_RT_idx));


%% Gratton Effects: Trial Type
cC_idx = fn_sequence_index('con', 'con', trial_info.condition_n, trial_info.block_n);
iC_idx = fn_sequence_index('inc', 'con', trial_info.condition_n, trial_info.block_n);
[~,C_pval] = ttest2(RTs(cC_idx==1), RTs(iC_idx==1));
cC_lab = ['cC (n=' num2str(sum(cC_idx)) ')'];
iC_lab = ['iC (n=' num2str(sum(iC_idx)) ')'];

cI_idx = fn_sequence_index('con', 'inc', trial_info.condition_n, trial_info.block_n);
iI_idx = fn_sequence_index('inc', 'inc', trial_info.condition_n, trial_info.block_n);
[~,I_pval] = ttest2(RTs(cI_idx==1), RTs(iI_idx==1));
cI_lab = ['cI (n=' num2str(sum(cI_idx)) ')'];
iI_lab = ['iI (n=' num2str(sum(iI_idx)) ')'];
% % Function Check
% for ix = 1:length(RTs)
%     if cC_idx(ix)==1
%         fprintf('%i: %s, %s, %s\n', ix,...
%             trial_info.condition_types{trial_info.condition_n(ix-1)},...
%             trial_info.condition_types{trial_info.condition_n(ix)},...
%             trial_info.condition_types{trial_info.condition_n(ix+1)});
%     end
% end

% Plot Boxplots
fig_name = strcat(data_id,'_RT_con_seq');
figure('Name',fig_name,'units','normalized','outerposition',[0 0 1 1],'Visible',vis);hold on;
subplot(1,2,1);hold on;
max_trl_num = max([sum(cC_idx) sum(iC_idx)]);
RT_plot = NaN([max_trl_num 2]);
RT_plot(1:sum(cC_idx),1) = RTs(cC_idx==1);
RT_plot(1:sum(iC_idx),2) = RTs(iC_idx==1);
boxplot(RT_plot,'Colors','br','Labels',{cC_lab,iC_lab},'Notch','on');
ylabel('Reaction Time (ms)');
title(strcat('Congruent Trials: p=',num2str(C_pval)));

subplot(1,2,2);hold on;
max_trl_num = max([sum(cI_idx) sum(iI_idx)]);
RT_plot = NaN([max_trl_num 2]);
RT_plot(1:sum(cI_idx),1) = RTs(cI_idx==1);
RT_plot(1:sum(iI_idx),2) = RTs(iI_idx==1);
boxplot(RT_plot,'Colors','br','Labels',{cI_lab,iI_lab},'Notch','on');
ylabel('Reaction Time (ms)');
title(strcat('Incongruent Trials: p=',num2str(I_pval)));

fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
if save_fig ==1
    fprintf('Saving %s\n',fig_filename);
    eval(['export_fig ' fig_filename]);
end

%% Gratton Effects: Within Trial Type Across Blocks
fig_name = strcat(data_id,'_RT_con_seq_by_block');
figure('Name',fig_name,'units','normalized','outerposition',[0 0 1 1],'Visible',vis);
for b_ix = 1:length(block_lab)
    cC_idx = fn_sequence_index(['con_' block_lab{b_ix}], ['con_' block_lab{b_ix}], trial_info.condition_n, trial_info.block_n);
    iC_idx = fn_sequence_index(['inc_' block_lab{b_ix}], ['con_' block_lab{b_ix}], trial_info.condition_n, trial_info.block_n);
    [~,C_pval] = ttest2(RTs(cC_idx==1), RTs(iC_idx==1));
    cC_lab = ['cC (n=' num2str(sum(cC_idx)) ')'];
    iC_lab = ['iC (n=' num2str(sum(iC_idx)) ')'];

    cI_idx = fn_sequence_index(['con_' block_lab{b_ix}], ['inc_' block_lab{b_ix}], trial_info.condition_n, trial_info.block_n);
    iI_idx = fn_sequence_index(['inc_' block_lab{b_ix}], ['inc_' block_lab{b_ix}], trial_info.condition_n, trial_info.block_n);
    [~,I_pval] = ttest2(RTs(cI_idx==1), RTs(iI_idx==1));
    cI_lab = ['cI (n=' num2str(sum(cI_idx)) ')'];
    iI_lab = ['iI (n=' num2str(sum(iI_idx)) ')'];

    % Plot Boxplots
    subplot(length(block_lab),2,2*b_ix-1);hold on;
    max_trl_num = max([sum(cC_idx) sum(iC_idx)]);
    RT_plot = NaN([max_trl_num 2]);
    RT_plot(1:sum(cC_idx),1) = RTs(cC_idx==1);
    RT_plot(1:sum(iC_idx),2) = RTs(iC_idx==1);
    boxplot(RT_plot,'Colors','br','Labels',{cC_lab,iC_lab});%,'Notch','on');
    ylabel('Reaction Time (ms)');
    title(strcat('"',block_lab{b_ix},'" Congruent Trials: p=',num2str(C_pval)));
    
    subplot(length(block_lab),2,2*b_ix);hold on;
    max_trl_num = max([sum(cI_idx) sum(iI_idx)]);
    RT_plot = NaN([max_trl_num 2]);
    RT_plot(1:sum(cI_idx),1) = RTs(cI_idx==1);
    RT_plot(1:sum(iI_idx),2) = RTs(iI_idx==1);
    boxplot(RT_plot,'Colors','br','Labels',{cI_lab,iI_lab});%,'Notch','on');
    ylabel('Reaction Time (ms)');
    title(strcat('"',block_lab{b_ix},'" Incongruent Trials: p=',num2str(I_pval)));
    
    clear cC_idx iC_idx cI_idx iI_idx C_pval I_pval cC_lab iC_lab cI_lab iI_lab
end

fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
if save_fig ==1
    fprintf('Saving %s\n',fig_filename);
    eval(['export_fig ' fig_filename]);
end

%% Gratton Effects: Error Likelihood by Proportion congruency



%%
% % Split RTs by condition type
% RT_cond = {};
% hist_cond = {};
% max_hist = 0;
% for lab_ix = 1:length(prop_con_lab)
%     eval(['RT_cond.' prop_con_lab{lab_ix} ' = fn_condition_index([prop_con_lab{lab_ix}], trial_info.condition_n);']);
%     eval(['hist_cond.' prop_con_lab{lab_ix} ' = histogram(RTs(RT_cond.' prop_con_lab{lab_ix} '),n_bins);']);
%     if max_hist < eval(['max(hist_cond.' prop_con_lab{lab_ix} '.Values)'])
%         eval(['max_hist = max(hist_cond.' prop_con_lab{lab_ix} '.Values);']);
%     end
% end

% %% Histograms per condition
% % Trial Type
% fig_name = strcat(data_id,'_RT_hist_trial_type');
% figure('Name',fig_name,'units','normalized','outerposition',[0 0 1 1],'Visible',vis);hold on;
% trial_RTs{1} = RTs(fn_condition_index([trial_lab{1}], trial_info.condition_n)==1);
% trial_RTs{2} = RTs(fn_condition_index([trial_lab{2}], trial_info.condition_n)==1);
% trial_RTs{3} = RTs(fn_condition_index([trial_lab{3}], trial_info.condition_n)==1);
% for lab_ix = 1:length(trial_lab)
%     histogram([trial_RTs{lab_ix}],n_edges,'FaceColor',trial_colors(lab_ix));
% end
% ylimits = ylim;
% for lab_ix = 1:length(trial_lab)
%     line([mean([trial_RTs{lab_ix}]) mean([trial_RTs{lab_ix}])], ylimits,...
%         'Color', [trial_colors{lab_ix}], 'LineWidth', line_w);
% end
% [~,pval] = ttest2([trial_RTs{1}],[trial_RTs{3}]);%,'Alpha',alpha);
% title(strcat('RT Histogram by Trial Type: p=',num2str(pval)));
% legend(trial_lab);
% 
% fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
% if save_fig ==1
%     fprintf('Saving %s\n',fig_filename);
%     eval(['export_fig ' fig_filename]);
% end
% 
% %% Compare trial types within each block type
% fig_name = strcat(data_id,'_RT_hist_trial_type_by_block');
% figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.6 1],'Visible',vis);
% ylimits = {};
% pval = {};
% for b_ix = 1:length(block_lab)
%     subplot(3,1,b_ix);hold on;
%     for lab_ix = 1:length(trial_lab)
%         block_RTs{lab_ix} = RTs(fn_condition_index([trial_lab{lab_ix} '_' block_lab{b_ix}],...
%             trial_info.condition_n)==1);
%         histogram([block_RTs{lab_ix}],n_edges,'FaceColor',trial_colors(lab_ix));
%     end
%     [~,pval{b_ix}] = ttest2([block_RTs{1}],[block_RTs{3}]);%,'Alpha',alpha);
%     ylimits{b_ix} = ylim;
% end
% max_ylim = max([ylimits{:}]);
% for b_ix = 1:length(block_lab)
%     subplot(3,1,b_ix);
%     ylim([0 max_ylim]);
%     for lab_ix = 1:length(trial_lab)
%         line([mean([block_RTs{lab_ix}]) mean([block_RTs{lab_ix}])], ylim,...
%             'Color', [trial_colors{lab_ix}], 'LineWidth', line_w);
%     end
%     legend(trial_lab);
%     title(strcat('RT Histogram for "',block_lab(b_ix),'" Blocks: p=',num2str(pval{b_ix})));
% end
% 
% fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
% if save_fig ==1
%     fprintf('Saving %s\n',fig_filename);
%     eval(['export_fig ' fig_filename]);
% end
% 
% %% Compare block effects within trial type
% fig_name = strcat(data_id,'_RT_hist_trial_type_prop_con');
% figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.6 1],'Visible',vis);
% for t_ix = 1:length(trial_lab)
%     subplot(3,1,t_ix);hold on;
%     trl_blk_RTs = {};
%     for b_ix = 1:length(block_lab)
%         trl_blk_RTs{b_ix} = RTs(fn_condition_index([trial_lab{t_ix} '_' block_lab{b_ix}],...
%             trial_info.condition_n)==1);
%         histogram([trl_blk_RTs{b_ix}],n_edges,'FaceColor',...
%             trial_colors(b_ix));%,'EdgeColor',trial_colors(b_ix));
%     end
%     for b_ix = 1:length(block_lab)
%         line([mean([trl_blk_RTs{b_ix}]) mean([trl_blk_RTs{b_ix}])], ylim, ...
%             'Color',[trial_colors{b_ix}], 'LineWidth', line_w);
%     end
%     fprintf('Trial type: %s\n',trial_lab{t_ix});
%     fprintf('Mean RT mcon = %f\n',mean([trl_blk_RTs{1}]));
%     fprintf('Mean RT same = %f\n',mean([trl_blk_RTs{2}]));
%     fprintf('Mean RT minc = %f\n',mean([trl_blk_RTs{3}]));
%     [~,pval] = ttest2([trl_blk_RTs{1}],[trl_blk_RTs{3}]);%,'Alpha',alpha);
%     legend(block_lab);
%     title(strcat('"',trial_lab(t_ix), '" RTs across block conditions: p=',num2str(pval)));
%     clear pval
% end
% 
% fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
% if save_fig ==1
%     fprintf('Saving %s\n',fig_filename);
%     eval(['export_fig ' fig_filename]);
% end
% 

