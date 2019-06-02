function B00a_SBJ_RT_violins(SBJ,plt_id,fig_vis,save_fig,fig_ftype)
%% Plot SBJ RTs for give conditions:
%       CNI: single histogram with RTs per trial type overlapping
%       PC: Plot 1- CNI hist with subplot per block
%             Plot 2- PC hists with subplot per trial type
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
if strcmp(fig_ftype,'eps')
    [~, app_dir] = fn_get_root_dir();
    addpath(genpath([app_dir 'export_fig-master/']));
end

%% Load data
% Load processing variables
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

late_RT_cut = 0.6;      %ms window before next stim to eliminate
%!!! implement!!! stat      = 'mean';         % do stats on 'mean' or 'median'

[cni_lab, cni_colors, ~] = fn_condition_label_styles('CNI');
[pc_lab, pc_colors, ~] = fn_condition_label_styles('PC');

cni_idx = zeros(size(trial_info.trial_n));
for cond_ix = 1:numel(cni_lab)
    cni_idx(logical(fn_condition_index(cni_lab{cond_ix},trial_info.condition_n))) = cond_ix;
end
pc_idx = zeros(size(trial_info.trial_n));
for cond_ix = 1:numel(pc_lab)
    pc_idx(logical(fn_condition_index(pc_lab{cond_ix},trial_info.condition_n))) = cond_ix;
end

% Create figure directory
fig_dir  = strcat([root_dir 'PRJ_Stroop/results/RTs/' SBJ '/']);
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Process data
% C = 1-3, N = 4-6, I = 7-9
% within those: EQ, MI, MC
rts     = trial_info.response_time; % converts sec to ms

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

%% ANOVA Statistics
[pval, table] = anovan(rts, {cni_idx, pc_idx},...
    'model', 'interaction', ...%, 'sstype', 3% 'continuous', strmatch('RT',w2.cond),
    'varnames', {'CNI','PC'}, 'display', 'off');
effect_lab = {'CNI','PC','CNI*PC'};

% OLD C vs. I t-test:
% [~,pval] = ttest2([trial_RTs{1}],[trial_RTs{3}]);%,'Alpha',alpha);

%% CNI Violins
fig_name = [SBJ '_RT_violins_CNI'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on; ax = gca;

% Separate data by CNI
rt_grps    = cell(size(cni_lab));
cni_legend = cell(size(cni_lab));
for cni_ix = 1:numel(cni_lab)
    rt_grps{cni_ix} = rts(cni_idx==cni_ix);
    cni_legend{cni_ix} = [cni_lab{cni_ix} ' (n=' num2str(sum(cni_idx==cni_ix)) ')'];
end

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), cni_lab, 'ViolinAlpha', violin_alpha);

% Adjust plot propeties
legend_obj = cell(size(cni_lab));
for cni_ix = 1:numel(cni_lab)
    % Change scatter colors to mark condition
    violins(cni_ix).ViolinColor = cni_colors{cni_ix};
    violins(cni_ix).BoxPlot.FaceColor = cni_colors{cni_ix};
    violins(cni_ix).EdgeColor = cni_colors{cni_ix};
    % Grab violin for legend
    legend_obj{cni_ix} = violins(cni_ix).ViolinPlot;
end
ax.FontSize        = tick_sz;
ax.YLabel.String   = 'Time (s)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['CNI RTs: p=' num2str(pval(1),'%.3f')];
ax.Title.FontSize  = title_sz;
legend([legend_obj{:}],cni_legend,'FontSize',leg_sz,'Location','best');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_ftype,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Proportion Congruent Bar Violins
fig_name = [SBJ '_RT_violins_PC'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on; ax = gca;

% Separate data by PC
rt_grps    = cell(size(pc_lab));
pc_legend = cell(size(pc_lab));
for pc_ix = 1:numel(pc_lab)
    rt_grps{pc_ix} = rts(pc_idx==pc_ix);
    pc_legend{pc_ix} = [pc_lab{pc_ix} ' (n=' num2str(sum(pc_idx==pc_ix)) ')'];
end

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), pc_lab, 'ViolinAlpha', violin_alpha);

% Adjust plot propeties
legend_obj = cell(size(pc_lab));
for pc_ix = 1:numel(pc_lab)
    % Change scatter colors to mark condition
    violins(pc_ix).ViolinColor = pc_colors{pc_ix};
    violins(pc_ix).BoxPlot.FaceColor = pc_colors{pc_ix};
    violins(pc_ix).EdgeColor = pc_colors{pc_ix};
    % Grab violin for legend
    legend_obj{pc_ix} = violins(pc_ix).ViolinPlot;
end
ax.FontSize        = tick_sz;
ax.YLabel.String   = 'Time (s)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['PC RTs: p=' num2str(pval(2),'%.3f')];
ax.Title.FontSize  = title_sz;
legend([legend_obj{:}],pc_legend,'FontSize',leg_sz,'Location','best');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_ftype,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Interaction Violins
fig_name = [SBJ '_RT_violins_CNI-PC_interaction'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.8 0.5],'Visible',fig_vis);
hold on; ax = gca;

% Separate data by PC and CNI
rt_grps    = cell([numel(cni_lab)*numel(pc_lab)+numel(cni_lab)-1 1]);
grp_lab    = cell(size(rt_grps));
grp_color  = cell(size(rt_grps));
grp_legend = cell(size(rt_grps));
grp_ix = 0;
for cni_ix = 1:numel(cni_lab)
    for pc_ix = 1:numel(pc_lab)
        grp_ix = grp_ix + 1;
        rt_grps{grp_ix}    = rts(pc_idx==pc_ix & cni_idx==cni_ix);
        grp_lab{grp_ix}    = [cni_lab{cni_ix} '-' pc_lab{pc_ix}];
        grp_color{grp_ix}  = pc_colors{pc_ix};
        grp_legend{grp_ix} = [grp_lab{grp_ix} ' (n=' num2str(numel(rt_grps{grp_ix})) ')'];
    end
    % Insert a blank set of data to create gap between block types for plotting
    if cni_ix~=numel(cni_lab) % only after first 2 conditions
        grp_ix = grp_ix + 1;
        rt_grps{grp_ix} = mean(rt_grps{grp_ix-1}); % non-empty column vector
        grp_lab{grp_ix} = '';
        grp_color{grp_ix} = [1 1 1];
    end
end
real_grp_idx  = ~cellfun(@isempty,grp_lab);

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), grp_lab, 'ViolinAlpha', violin_alpha);

% Adjust plot propeties
grp_object = cell(size(grp_lab));
for grp_ix = 1:numel(rt_grps)
    % Change scatter colors to mark condition
    violins(grp_ix).ViolinColor = grp_color{grp_ix};
    violins(grp_ix).EdgeColor = grp_color{grp_ix};
    
    % Make spacers transparent
    if ~real_grp_idx(grp_ix)
        violins(grp_ix).MedianPlot.MarkerEdgeColor = grp_color{grp_ix};
        violins(grp_ix).MedianColor = grp_color{grp_ix};
    else
        % Color boxplot bar according to PC block
        cni_ix = find(~cellfun(@isempty,strfind(cni_lab,grp_lab{grp_ix}(1))));
        violins(grp_ix).BoxColor  = cni_colors{cni_ix};
        violins(grp_ix).EdgeColor = cni_colors{cni_ix};
    end
    % Grab object for legend
    grp_object{grp_ix} = violins(grp_ix).ViolinPlot;
end
ax.FontSize        = tick_sz;
ax.YLabel.String   = 'Time (s)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['CNI RTs: p=' num2str(pval(1),'%.3f') '; PC: p=' num2str(pval(2),'%.3f')...
                      '; CNI*PC p=' num2str(pval(3),'%.3f')];
ax.Title.FontSize  = title_sz;
% legend([grp_object{real_grp_idx}],grp_legend{real_grp_idx},'FontSize',leg_sz,'Location','best');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_ftype,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Neutral Trials across MC vs. MI Blocks
% n_ix = find(strcmp(cni_lab,'N'));
% [pval, table] = anova1(rts(cni_idx==n_ix), pc_idx(cni_idx==n_ix), 'off');

fig_name = [SBJ '_RT_violins_PC_N'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on; ax = gca;

% Compare only MI and MC
n_ix = find(strcmp(cni_lab,'N'));
mc_ix = find(strcmp(pc_lab,'MC'));
mi_ix = find(strcmp(pc_lab,'MI'));
n_pc_lab = pc_lab([mc_ix mi_ix]);
n_pc_colors = pc_colors([mc_ix mi_ix]);

% Separate data by PC
rt_grps     = cell(size(n_pc_lab));
n_pc_legend = cell(size(n_pc_lab));
for pc_ix = 1:numel(n_pc_lab)
    rt_grps{pc_ix} = rts(pc_idx==find(strcmp(pc_lab,n_pc_lab{pc_ix})) & cni_idx==n_ix);
    n_pc_legend{pc_ix} = [n_pc_lab{pc_ix} ' (n=' num2str(sum(pc_idx==find(strcmp(pc_lab,n_pc_lab{pc_ix})) & cni_idx==n_ix)) ')'];
end

% Two-sample t-test
[h, pval, ci, stats] = ttest2(rt_grps{strcmp(n_pc_lab,'MC')},...
                              rt_grps{strcmp(n_pc_lab,'MI')});
% % Two-sample t-test, MI RTs < MC RTs
% [h, pval, ci, stats] = ttest2(rt_grps{strcmp(n_pc_lab,'MC')},...
%                               rt_grps{strcmp(n_pc_lab,'MI')},'Tail','right');

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), n_pc_lab, 'ViolinAlpha', violin_alpha);

% Adjust plot properties
legend_obj = cell(size(n_pc_lab));
for pc_ix = 1:numel(n_pc_lab)
    % Change scatter colors to mark condition
    violins(pc_ix).ViolinColor = n_pc_colors{pc_ix};
    violins(pc_ix).BoxPlot.FaceColor = n_pc_colors{pc_ix};
    violins(pc_ix).EdgeColor = n_pc_colors{pc_ix};
    % Grab violin for legend
    legend_obj{pc_ix} = violins(pc_ix).ViolinPlot;
end
ax.FontSize        = tick_sz;
ax.YLabel.String   = 'Time (s)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['N PC RTs: p=' num2str(pval,'%.3f')];
ax.Title.FontSize  = title_sz;
legend([legend_obj{:}],n_pc_legend,'FontSize',leg_sz,'Location','best');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_ftype,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

end