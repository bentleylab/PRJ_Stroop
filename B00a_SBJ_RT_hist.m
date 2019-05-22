function B00a_SBJ_RT_hist(SBJ,plt_id,fig_vis,save_fig,fig_type)
%% Plot SBJ RTs for give conditions:
%       CNI: single histogram with RTs per trial type overlapping
%       pcon: Plot 1- CNI hist with subplot per block
%             Plot 2- pcon hists with subplot per trial type
% INPUTS:
%   SBJ [str] - 
%   plt_id [str] - defines histogram and line parameters
%   fig_vis ['on'/'off']
%   save_fig [0/1]
%   fig_type [str] - file extension for saving figure

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
if strcmp(fig_type,'eps')
    [~, app_dir] = fn_get_root_dir();
    addpath(genpath([app_dir 'export_fig-master/']));
end

%% Load data
% Load processing variables
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/',SBJ,'_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/',plt_id,'_vars.m']);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

late_RT_cut = 0.6;      %ms window before next stim to eliminate
%!!! implement!!! stat      = 'mean';         % do stats on 'mean' or 'median'

[cni_lab, cni_colors, ~] = fn_condition_label_styles('CNI');
[pcon_lab, pcon_colors, ~] = fn_condition_label_styles('pcon');

cni_idx = zeros(size(trial_info.trial_n));
for cond_ix = 1:numel(cni_lab)
    cni_idx(logical(fn_condition_index(cni_lab{cond_ix},trial_info.condition_n))) = cond_ix;
end
pcon_idx = zeros(size(trial_info.trial_n));
for cond_ix = 1:numel(pcon_lab)
    pcon_idx(logical(fn_condition_index(pcon_lab{cond_ix},trial_info.condition_n))) = cond_ix;
end
% n_bins        = 50;
% line_w        = 2;
% trial_lab     = {'con', 'neu', 'inc'};
% block_lab     = {'mcon', 'same', 'minc'};
% trial_colors   = {'b', 'k', 'r'};    % colors for cond_lab plotting
% block_colors  = {repmat(0.8,3,1), repmat(0.5,3,1), repmat(0.2,3,1)};    % colors for [mcon, same, minc]
% prop_con_lab  = {'con_mcon', 'con_same', 'con_minc', 'neu_mcon', 'neu_same',...
%     'neu_minc', 'inc_mcon', 'inc_same', 'inc_minc'};

% Create figure directory
fig_dir  = strcat([root_dir 'PRJ_Stroop/results/RTs/' SBJ '/']);
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Process data
% con = 1-3, neu = 4-6, inc = 7-9
% within those: same, mic, mcon
% trial_type = NaN(size(trial_info.condition_n));
rts     = trial_info.response_time; % converts sec to ms
edges = linspace(min(rts),max(rts),n_bins);

% Check for RTs overlapping with stim onset or baseline
late_RT_idx = zeros(size(trial_info.trial_n));
for t_ix = 2:numel(trial_info.resp_onset)
    if trial_info.word_onset(t_ix)-late_RT_cut*trial_info.sample_rate <= trial_info.resp_onset(t_ix-1)
        late_RT_idx(t_ix) = 1;
    end
end
fprintf('%i late trials detected\n',sum(late_RT_idx));

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

%% CNI Histogram
fig_name = [SBJ '_RT_hist_CNI'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on;

% Plot Histograms
cni_legend = cell(size(cni_lab));
for cond_ix = [2 1 3]   % plot N first so C is on top if they overlap
    histogram(rts(cni_idx==cond_ix),edges,'FaceColor',cni_colors{cond_ix},'FaceAlpha',hist_alpha);
    cni_legend{cond_ix} = [cni_lab{cond_ix} ' (n=' num2str(sum(cni_idx==cond_ix)) ')'];
end
% Plot Means
for cond_ix = [2 1 3]   % plot N first so C is on top if they overlap
    line([mean(rts(cni_idx==cond_ix)) mean(rts(cni_idx==cond_ix))], ylim,...
        'Color', cni_colors{cond_ix}, 'LineWidth', line_width, 'LineStyle', line_style);
end
% [~,pval] = ttest2([trial_RTs{1}],[trial_RTs{3}]);%,'Alpha',alpha);
title('CNI RT Histogram');% : p=',num2str(pval)));
xlabel('Time (s)');
legend(cni_legend);

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_type];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_type,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Proportion Congruent Bar Violins
fig_name = [SBJ '_RT_violin_pcon'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.8 0.5],'Visible',fig_vis);
hold on;

% Separate data by pcon and CNI
rt_grps = cell([numel(pcon_lab)*numel(cni_lab)+numel(pcon_lab)-1 1]);
grp_lab = cell([numel(pcon_lab)*numel(cni_lab)+numel(pcon_lab)-1 1]);
grp_color = cell([numel(pcon_lab)*numel(cni_lab)+numel(pcon_lab)-1 1]);
grp_ix = 0;
for pcon_ix = 1:numel(pcon_lab)
    for cni_ix = 1:numel(cni_lab)
        grp_ix = grp_ix + 1;
        rt_grps{grp_ix} = rts(pcon_idx==pcon_ix & cni_idx==cni_ix);
        grp_lab{grp_ix} = [pcon_lab{pcon_ix} '-' cni_lab{cni_ix}];
        grp_color{grp_ix} = cni_colors{cni_ix};
    end
    % Insert a blank set of data to create gap between block types for plotting
    if pcon_ix~=numel(pcon_lab)
        grp_ix = grp_ix + 1;
        rt_grps{grp_ix} = mean(rt_grps{grp_ix-1}); % non-empty column vector
        grp_lab{grp_ix} = '';
        grp_color{grp_ix} = [1 1 1];
    end
end

% ANOVA Stats


% Plot Violins
violins = violinplot(padcat(rt_grps{:}), grp_lab, 'ViolinAlpha', violin_alpha);

% Adjust plot propeties
for grp_ix = 1:numel(rt_grps)
    % Change scatter colors to mark condition
    violins(grp_ix).ViolinColor = grp_color{grp_ix};
    violins(grp_ix).BoxPlot.FaceColor = grp_color{grp_ix};
    violins(grp_ix).EdgeColor = grp_color{grp_ix};
    
    % Make spacers transparent
    if isempty(grp_lab{grp_ix})
        violins(grp_ix).MedianPlot.MarkerEdgeColor = grp_color{grp_ix};
        violins(grp_ix).MedianColor = grp_color{grp_ix};
    end
end
ylabel('Time (s)');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_type];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_type,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Compare trial types within each block type
% fig_name = [SBJ '_RT_hist_pcon_TwiB'];
% figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.6 1],'Visible',fig_vis);
% % pval = zeros(size(pcon_lab));
% 
% % Plot Histograms
% ylims = zeros([numel(pcon_lab) 2]);
% pcon_legend = cell([numel(pcon_lab) numel(cni_lab)]);
% for cond_ix = 1:length(pcon_lab)
%     subplot(numel(pcon_lab),1,cond_ix);hold on;
%     for cni_ix = [2 1 3]
%         histogram(RTs(all([pcon_idx==cond_ix cni_idx==cni_ix],2)),edges,...
%             'FaceColor',cni_colors{cni_ix},'FaceAlpha',hist_alpha);
%         pcon_legend{cond_ix,cni_ix} = [cni_lab{cni_ix} ' (n=' ...
%                                        num2str(sum(all([pcon_idx==cond_ix cni_idx==cni_ix],2))) ')'];
%     end
%     ylims(cond_ix,:) = get(gca,'YLim');
% %     [~,pval{cond_ix}] = ttest2([block_RTs{cond_ix,1}],[block_RTs{cond_ix,3}]);%,'Alpha',alpha);
% end
% 
% % Plot Condition Means
% max_ylim = max(ylims(:));
% for cond_ix = 1:length(pcon_lab)
%     subplot(3,1,cond_ix);
%     for cni_ix = [2 1 3]
%         line([mean(RTs(all([pcon_idx==cond_ix cni_idx==cni_ix],2)))...
%               mean(RTs(all([pcon_idx==cond_ix cni_idx==cni_ix],2)))],...
%              [0 max_ylim], 'Color', cni_colors{cni_ix}, 'LineWidth', line_width, 'LineStyle', line_style);
%     end
%     legend(pcon_legend{cond_ix,:});
%     xlabel('Time (s)');
%     title(['CNI Histogram for ' pcon_lab{cond_ix} ' Blocks']);%: p=',num2str(pval{cond_ix})));
% end
% 
% fig_fname = [fig_dir fig_name '.' fig_type];
% if save_fig
%     fprintf('Saving %s\n',fig_fname);
%     if strcmp(fig_type,'eps')
%         eval(['export_fig ' fig_fname]);
%     else
%         saveas(gcf,fig_fname);
%     end
% end
% 
% %% Compare block effects within trial type
% fig_name = [SBJ '_RT_hist_pcon_BwiT'];
% figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.6 1],'Visible',fig_vis);
% % pval = zeros(size(pcon_lab));
% 
% % Plot Histograms
% ylims = zeros([numel(cni_lab) 2]);
% pcon_legend = cell([numel(pcon_lab) numel(cni_lab)]);
% for cni_ix = [2 1 3]
%     subplot(numel(cni_lab),1,cni_ix);hold on;
%     for cond_ix = 1:length(pcon_lab)
%         histogram(RTs(all([pcon_idx==cond_ix cni_idx==cni_ix],2)),edges,...
%             'FaceColor',pcon_colors{cond_ix},'FaceAlpha',hist_alpha);
%         pcon_legend{cond_ix,cni_ix} = [pcon_lab{cond_ix} ' (n=' ...
%                                        num2str(sum(all([pcon_idx==cond_ix cni_idx==cni_ix],2))) ')'];
%     end
%     ylims(cni_ix,:) = get(gca,'YLim');
% %     [~,pval{cond_ix}] = ttest2([block_RTs{cond_ix,1}],[block_RTs{cond_ix,3}]);%,'Alpha',alpha);
% end
% 
% % Plot Condition Means
% max_ylim = max(ylims(:));
% for cni_ix = [2 1 3]
%     subplot(numel(cni_lab),1,cni_ix);
%     for cond_ix = 1:length(pcon_lab)
%         line([mean(RTs(all([pcon_idx==cond_ix cni_idx==cni_ix],2)))...
%               mean(RTs(all([pcon_idx==cond_ix cni_idx==cni_ix],2)))],...
%              [0 max_ylim], 'Color', pcon_colors{cond_ix}, 'LineWidth', line_width, 'LineStyle', line_style);
%     end
%     legend(pcon_legend{:,cni_ix});
%     xlabel('Time (s)');
%     title([cni_lab{cni_ix} ' Trials across pcon blocks']);%: p=',num2str(pval{cond_ix})));
% end
% 
% fig_fname = [fig_dir fig_name '.' fig_type];
% if save_fig
%     fprintf('Saving %s\n',fig_fname);
%     if strcmp(fig_type,'eps')
%         eval(['export_fig ' fig_fname]);
%     else
%         saveas(gcf,fig_fname);
%     end
% end

end