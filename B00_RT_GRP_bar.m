%% RT Behavioral analysis- draw box plots per SBJ and for Group
clear all; close all
% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
%addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));

% Analysis Parameters
conditions = 'CNI';
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};

late_RT_cut = 0.6;      % window in sec before next stim to eliminate
% Plotting parameters
save_fig      = 1;
vis           = 'on';
n_bins        = 50;
line_w        = 2;
[cond_lab, cond_colors, cond_styles] = fn_condition_label_styles(conditions);
% block_lab     = {'mcon', 'same', 'minc'};
% trial_colors   = {'b', 'k', 'r'};    % colors for cond_lab plotting
% block_colors  = {repmat(0.8,3,1), repmat(0.5,3,1), repmat(0.2,3,1)};    % colors for [mcon, same, minc]
% prop_con_lab  = {'con_mcon', 'con_same', 'con_minc', 'neu_mcon', 'neu_same',...
%     'neu_minc', 'inc_mcon', 'inc_same', 'inc_minc'};
fig_type      = 'png';

% Process parameters
fig_dir  = strcat([root_dir 'PRJ_Stroop/results/RTs/GRP/']);
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Load data
RTs       = cell([numel(SBJs) numel(cond_lab)]);
RT_means  = NaN([numel(SBJs) numel(cond_lab)]);
RT_vars   = NaN([numel(SBJs) numel(cond_lab)]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    
    % con = 1-3, neu = 4-6, inc = 7-9
    % within those: same, mic, 
    for cond_ix = 1:numel(cond_lab)
        RTs{sbj_ix, cond_ix} = trial_info.response_time(logical(fn_condition_index(cond_lab{cond_ix},...
                                                                        trial_info.condition_n)));
        RT_means(sbj_ix,cond_ix) = mean(RTs{sbj_ix,cond_ix});
        RT_vars(sbj_ix,cond_ix)  = std(RTs{sbj_ix,cond_ix});
    end
    
    % Check for RTs overlapping with stim onset or baseline
    late_RT_count = 0;
    for t_ix = 2:length(trial_info.resp_onset)
        if trial_info.word_onset(t_ix)-late_RT_cut*trial_info.sample_rate <= trial_info.resp_onset(t_ix-1)
            late_RT_count = late_RT_count + 1;
        end
    end
    fprintf('%s: %i late trials detected\n',SBJ,length(late_RT_count));
    clear trial_info
end

% % Add means across groups
% for cond_ix = 1:numel(cond_lab)
%     RTs{end-numel(cond_lab)+cond_ix} = RT_means(:,cond_ix);
%     RT_groups(end-numel(cond_lab)+cond_ix,:) = [numel(SBJs)+1, cond_ix];
% end
%     
% % Reformat RTs to matrix with 
% max_size = max(cellfun(@numel,RTs));
% RT_mat = NaN([max_size numel(RTs)]);
% for col_ix = 1:numel(RTs)
%     RT_mat(1:numel(RTs{col_ix})) = RTs{col_ix};
% end

%% Histograms per condition
% Trial Type
fig_name = strcat('GRP_RT_box_trial_type');
figure('Name',fig_name,'units','normalized','outerposition',[0 0 1 1],'Visible',vis);hold on

% Plot Group means
b = [];
for cond_ix = 1:numel(cond_lab)
    b(cond_ix) = bar(cond_ix,mean(RT_means(:,cond_ix)),cond_colors{cond_ix});
end
% Plot Individuals
for sbj_ix = 1:numel(SBJs)
    !!! stopped here!!!
    % Plot individual subject percentages as scatters on top
    s = scatter(scat_offsets+bar_offsets(groi_ix)+an_ix,...
        scat_vars{an_ix}(:,groi_ix)./elec_g_count(:,groi_ix),50,'k*');
end
if strcmp(event,'stim')
    legend([b{1},s],groi_list{:},'Individuals','Location','northwest');
else
    legend([b{1},s],groi_list{:},'Individuals','Location','northeast');
end


trial_RTs{1} = RTs(fn_condition_index([cond_lab{1}], trial_info.condition_n)==1);
trial_RTs{2} = RTs(fn_condition_index([cond_lab{2}], trial_info.condition_n)==1);
trial_RTs{3} = RTs(fn_condition_index([cond_lab{3}], trial_info.condition_n)==1);
for lab_ix = 1:length(cond_lab)
    histogram([trial_RTs{lab_ix}],n_edges,'FaceColor',trial_colors(lab_ix));
end
ylimits = ylim;
for lab_ix = 1:length(cond_lab)
    line([mean([trial_RTs{lab_ix}]) mean([trial_RTs{lab_ix}])], ylimits,...
        'Color', [trial_colors{lab_ix}], 'LineWidth', line_w);
end
[~,pval] = ttest2([trial_RTs{1}],[trial_RTs{3}]);%,'Alpha',alpha);
title(strcat('RT Histogram by Trial Type: p=',num2str(pval)));
trial_legend = {};
for lab_ix = 1:length(cond_lab)
    trial_legend{lab_ix} = [cond_lab{lab_ix} '-' num2str(length(trial_RTs{lab_ix}))];
end
legend(trial_legend);

fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
if save_fig ==1
    fprintf('Saving %s\n',fig_filename);
    eval(['export_fig ' fig_filename]);
end

%% Compare trial types within each block type
fig_name = strcat(SBJ,'_RT_hist_trial_type_by_block');
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.6 1],'Visible',vis);
ylimits = {};
pval = {};
for b_ix = 1:length(block_lab)
    subplot(3,1,b_ix);hold on;
    for lab_ix = 1:length(cond_lab)
        block_RTs{b_ix,lab_ix} = RTs(fn_condition_index([cond_lab{lab_ix} '_' block_lab{b_ix}],...
            trial_info.condition_n)==1);
        histogram([block_RTs{b_ix,lab_ix}],n_edges,'FaceColor',trial_colors(lab_ix));
    end
    [~,pval{b_ix}] = ttest2([block_RTs{b_ix,1}],[block_RTs{b_ix,3}]);%,'Alpha',alpha);
    ylimits{b_ix} = ylim;
end
max_ylim = max([ylimits{:}]);
for b_ix = 1:length(block_lab)
    subplot(3,1,b_ix);
    ylim([0 max_ylim]);
    trial_legend = {};
    for lab_ix = 1:length(cond_lab)
        line([mean([block_RTs{b_ix,lab_ix}]) mean([block_RTs{b_ix,lab_ix}])], ylim,...
            'Color', [trial_colors{lab_ix}], 'LineWidth', line_w);
        trial_legend{lab_ix} = [cond_lab{lab_ix} '-' num2str(length(block_RTs{b_ix,lab_ix}))];
    end
    legend(trial_legend);
    title(strcat('RT Histogram for "',block_lab(b_ix),'" Blocks: p=',num2str(pval{b_ix})));
end

fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
if save_fig ==1
    fprintf('Saving %s\n',fig_filename);
    eval(['export_fig ' fig_filename]);
end

%% Compare block effects within trial type
fig_name = strcat(SBJ,'_RT_hist_trial_type_prop_con');
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.6 1],'Visible',vis);
for t_ix = 1:length(cond_lab)
    subplot(3,1,t_ix);hold on;
    trl_blk_RTs = {};
    for b_ix = 1:length(block_lab)
        trl_blk_RTs{b_ix} = RTs(fn_condition_index([cond_lab{t_ix} '_' block_lab{b_ix}],...
            trial_info.condition_n)==1);
        histogram([trl_blk_RTs{b_ix}],n_edges,'FaceColor',...
            trial_colors(b_ix));%,'EdgeColor',trial_colors(b_ix));
    end
    block_legend = {};
    for b_ix = 1:length(block_lab)
        line([mean([trl_blk_RTs{b_ix}]) mean([trl_blk_RTs{b_ix}])], ylim, ...
            'Color',[trial_colors{b_ix}], 'LineWidth', line_w);
        block_legend{b_ix} = [block_lab{b_ix} '-' num2str(length(trl_blk_RTs{b_ix}))];
    end
    fprintf('Trial type: %s\n',cond_lab{t_ix});
    fprintf('Mean RT mcon = %f\n',mean([trl_blk_RTs{1}]));
    fprintf('Mean RT same = %f\n',mean([trl_blk_RTs{2}]));
    fprintf('Mean RT minc = %f\n',mean([trl_blk_RTs{3}]));
    [~,pval] = ttest2([trl_blk_RTs{1}],[trl_blk_RTs{3}]);%,'Alpha',alpha);
    legend(block_legend);
    title(strcat('"',cond_lab(t_ix), '" RTs across block conditions: p=',num2str(pval)));
    clear pval
end

fig_filename = strcat(fig_dir,fig_name,'.',fig_type);
if save_fig ==1
    fprintf('Saving %s\n',fig_filename);
    eval(['export_fig ' fig_filename]);
end


%% Gratton Effects
% Trial Type

% Within Trial Type Across Blocks

% Error Likelihood by Proportion congruency

