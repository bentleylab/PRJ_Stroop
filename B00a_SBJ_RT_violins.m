function B00a_SBJ_RT_violins(SBJ,plt_id,fig_vis,save_fig,fig_type)
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

%% ANOVA Statistics
[pval, table] = anovan(rts, {cni_idx, pcon_idx},...
    'model', 'interaction', ...%, 'sstype', 3% 'continuous', strmatch('RT',w2.cond),
    'varnames', {'CNI','pcon'}, 'display', 'off');
effect_lab = {'CNI','pcon','CNI*pcon'};

% OLD C vs. I t-test:
% [~,pval] = ttest2([trial_RTs{1}],[trial_RTs{3}]);%,'Alpha',alpha);

%% CNI Violins
fig_name = [SBJ '_RT_violins_CNI'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on;

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
ylabel('Time (s)');
title(['CNI RTs: p=' num2str(pval(1),'%.3f')]);
legend([legend_obj{:}],cni_legend,'Location','best');

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
fig_name = [SBJ '_RT_violins_pcon'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on;

% Separate data by CNI
rt_grps    = cell(size(pcon_lab));
pcon_legend = cell(size(pcon_lab));
for pcon_ix = 1:numel(pcon_lab)
    rt_grps{pcon_ix} = rts(pcon_idx==pcon_ix);
    pcon_legend{pcon_ix} = [pcon_lab{pcon_ix} ' (n=' num2str(sum(pcon_idx==pcon_ix)) ')'];
end

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), pcon_lab, 'ViolinAlpha', violin_alpha);

% Adjust plot propeties
legend_obj = cell(size(pcon_lab));
for pcon_ix = 1:numel(pcon_lab)
    % Change scatter colors to mark condition
    violins(pcon_ix).ViolinColor = pcon_colors{pcon_ix};
    violins(pcon_ix).BoxPlot.FaceColor = pcon_colors{pcon_ix};
    violins(pcon_ix).EdgeColor = pcon_colors{pcon_ix};
    % Grab violin for legend
    legend_obj{pcon_ix} = violins(pcon_ix).ViolinPlot;
end
ylabel('Time (s)');
title(['pcon RTs: p=' num2str(pval(2),'%.3f')]);
legend([legend_obj{:}],pcon_legend,'Location','best');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_type];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_type,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Interaction Violins
fig_name = [SBJ '_RT_violins_CNI-pcon_interaction'];
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.8 0.5],'Visible',fig_vis);
hold on;

% Separate data by pcon and CNI
rt_grps    = cell([numel(pcon_lab)*numel(cni_lab)+numel(pcon_lab)-1 1]);
grp_lab    = cell(size(rt_grps));
grp_color  = cell(size(rt_grps));
grp_legend = cell(size(rt_grps));
grp_ix = 0;
for pcon_ix = 1:numel(pcon_lab)
    for cni_ix = 1:numel(cni_lab)
        grp_ix = grp_ix + 1;
        rt_grps{grp_ix}    = rts(pcon_idx==pcon_ix & cni_idx==cni_ix);
        grp_lab{grp_ix}    = [pcon_lab{pcon_ix} '-' cni_lab{cni_ix}];
        grp_color{grp_ix}  = cni_colors{cni_ix};
        grp_legend{grp_ix} = [grp_lab{grp_ix} ' (n=' num2str(numel(rt_grps{grp_ix})) ')'];
    end
    % Insert a blank set of data to create gap between block types for plotting
    if pcon_ix~=numel(pcon_lab) % only after first 2 conditions
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
        % Color boxplot bar according to pcon block
        pcon_ix = find(~cellfun(@isempty,strfind(pcon_lab,grp_lab{grp_ix}(1:4))));
        violins(grp_ix).BoxColor = pcon_colors{pcon_ix};
        violins(grp_ix).EdgeColor = pcon_colors{pcon_ix};
    end
    % Grab object for legend
    grp_object{grp_ix} = violins(grp_ix).ViolinPlot;
end
ylabel('Time (s)');
title(['CNI RTs: p=' num2str(pval(1),'%.3f') '; pcon: p=' num2str(pval(2),'%.3f') '; CNI*pcon p=' num2str(pval(3),'%.3f')]);
legend([grp_object{real_grp_idx}],grp_legend{real_grp_idx},'Location','best');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_type];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_type,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

end