function B00b_GRP_RT_violins_norm(SBJs,plt_id,save_fig,fig_vis,fig_ftype)
%% RT histogram: overlapping condition hist after normalizing group RTs
%% Set Up Directories
% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
%addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));

%% Analysis Parameters
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);
late_RT_cut = 0.6;      %ms window before next stim to eliminate
%!!! implement!!! stat      = 'mean';         % do stats on 'mean' or 'median'
%!!! add outlier threshold (z-score = 3)

[cni_lab, cni_colors, ~] = fn_condition_label_styles('CNI');
[pc_lab, pc_colors, ~]   = fn_condition_label_styles('PC');

% Process parameters
fig_dir  = strcat([root_dir 'PRJ_Stroop/results/RTs/GRP/']);
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

%% Load data
cni_idx  = cell(size(SBJs));
pc_idx   = cell(size(SBJs));
rts      = cell(size(SBJs));
late_idx = cell(size(SBJs));
sbj_idx  = cell(size(SBJs));
% RT_means  = NaN([numel(SBJs) numel(cond_lab)]);
% RT_vars   = NaN([numel(SBJs) numel(cond_lab)]);
for s_ix = 1:numel(SBJs)
    SBJ = SBJs{s_ix};
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    
    % Normalize RTs within SBJ
    rts{s_ix} = zscore(trial_info.response_time);
    
    % Condition index within SBJ
    cni_idx{s_ix} = zeros(size(trial_info.trial_n));
    for cond_ix = 1:numel(cni_lab)
        cni_idx{s_ix}(logical(fn_condition_index(cni_lab{cond_ix},...
                        trial_info.condition_n))) = cond_ix;
    end
    pc_idx{s_ix} = zeros(size(trial_info.trial_n));
    for cond_ix = 1:numel(pc_lab)
        pc_idx{s_ix}(logical(fn_condition_index(pc_lab{cond_ix},...
                        trial_info.condition_n))) = cond_ix;
    end
    
    % Create SBJ factor for ANOVA
    sbj_idx{s_ix} = ones(size(rts{s_ix}))*s_ix;
    
    % Check for RTs overlapping with stim onset or baseline
    late_idx{s_ix} = false(size(rts{s_ix}));
    for t_ix = 2:length(trial_info.resp_onset)
        if trial_info.word_onset(t_ix)-trial_info.resp_onset(t_ix-1) ...
                    <= late_RT_cut*trial_info.sample_rate 
            late_idx{s_ix}(t_ix) = true;
        end
    end
    if any(late_idx{s_ix})
        fprintf('%s: %i late trials detected\n',SBJ,sum(late_idx{s_ix}));
    end
    rts{s_ix}(late_idx{s_ix}) = [];
    cni_idx{s_ix}(late_idx{s_ix}) = [];
    pc_idx{s_ix}(late_idx{s_ix}) = [];
    clear trial_info
end

% %% Normalize RTs
% RTs_all  = vertcat(rts{:});
% RTs_norm = cell([numel(SBJs) numel(cond_lab)]);
% for s_ix = 1:numel(SBJs)
%     for cond_ix = 1:numel(cond_lab)
%         RTs_norm{s_ix, cond_ix} = (rts{s_ix, cond_ix}-mean(RTs_all))/std(RTs_all);
%     end
% end

%% ANOVA Statistics
factor_lab = {'CNI','PC','SBJ'};
% model_terms = [1 0 0; 0 1 0; 1 1 0];
[pval, table, stats, terms] = anovan(vertcat(rts{:}),...
    {vertcat(cni_idx{:}), vertcat(pc_idx{:}), vertcat(sbj_idx{:})},...
    'model', 'interaction', 'random', 3,...%, 'sstype', 3% 'continuous', strmatch('RT',w2.cond),
    'varnames', factor_lab, 'display', 'off');
% Output pvals = [CNI, PC, SBJ, CNI*PC, CNI*SBJ, PC*SBJ]

%% CNI Violins
% Trial Type
fig_name = 'GRP_RT_violins_norm_CNI';
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on; ax = gca;

% Separate data by CNI
rt_grps    = cell(size(cni_lab));
cni_legend = cell(size(cni_lab));
for cni_ix = 1:numel(cni_lab)
    for s_ix = 1:numel(SBJs)
        rt_grps{cni_ix} = [rt_grps{cni_ix}; rts{s_ix}(cni_idx{s_ix}==cni_ix)];
    end
    cni_legend{cni_ix} = [cni_lab{cni_ix} ' (n=' num2str(numel(rt_grps{cni_ix})) ')'];
end

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), cni_lab,...
                     'ViolinAlpha', violin_alpha, 'ShowData', false);

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
ax.YLabel.String   = 'RT (z-score)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['Group (n=' num2str(numel(SBJs)) ') RTs: CNI p=' num2str(pval(1),'%.3f')];
ax.Title.FontSize  = title_sz;
legend([legend_obj{:}],cni_legend,'FontSize',leg_sz,'Location','best');

% Save Figure
fig_fname = [fig_dir,fig_name,'.',fig_ftype];
if save_fig
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_ftype,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Proportion Congruent Bar Violins
fig_name = 'GRP_RT_violins_PC';
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on; ax = gca;

% Separate data by PC
rt_grps    = cell(size(pc_lab));
pc_legend = cell(size(pc_lab));
for pc_ix = 1:numel(pc_lab)
    for s_ix = 1:numel(SBJs)
        rt_grps{pc_ix} = [rt_grps{pc_ix}; rts{s_ix}(pc_idx{s_ix}==pc_ix)];
    end
    pc_legend{pc_ix} = [pc_lab{pc_ix} ' (n=' num2str(numel(rt_grps{pc_ix})) ')'];
end

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), pc_lab,...
                     'ViolinAlpha', violin_alpha, 'ShowData', false);

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
ax.YLabel.String   = 'RT (z-score)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['Group (n=' num2str(numel(SBJs)) ') RTs: PC p=' num2str(pval(2),'%.3f')];
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
fig_name = 'GRP_RT_violins_CNI-PC_interaction';
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
        for s_ix = 1:numel(SBJs)
            rt_grps{grp_ix} = [rt_grps{grp_ix};...
                               rts{s_ix}(pc_idx{s_ix}==pc_ix & cni_idx{s_ix}==cni_ix)];
        end
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
violins = violinplot(padcat(rt_grps{:}), grp_lab,...
                     'ViolinAlpha', violin_alpha, 'ShowData', false);

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
ax.YLabel.String   = 'RT (z-score)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['Group (n=' num2str(numel(SBJs)) ') RTs: CNI p=' num2str(pval(1),'%.3f')...
                      '; PC p=' num2str(pval(2),'%.3f') '; CNI*PC p='  num2str(pval(4),'%.3f')];
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

%% Compute Stroop Effect Across Blocks
fig_name = 'GRP_RT_violins_PC_stroop';
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on; ax = gca;

stroop = zeros([numel(SBJs) numel(pc_lab)]);
design = zeros([numel(SBJs) numel(pc_lab)]);
con_ix = find(strcmp(cni_lab,'C'));
inc_ix = find(strcmp(cni_lab,'I'));
for pc_ix = 1:numel(pc_lab)
    for s_ix = 1:numel(SBJs)
        con_rts = rts{s_ix}(pc_idx{s_ix}==pc_ix & cni_idx{s_ix}==con_ix);
        inc_rts = rts{s_ix}(pc_idx{s_ix}==pc_ix & cni_idx{s_ix}==inc_ix);
        stroop(s_ix,pc_ix) = mean(inc_rts)-mean(con_rts);
    end
    design(:,pc_ix) = pc_ix;
end

% ANOVA Statistics
[pval, table] = anovan(stroop(:), design(:),...
    'model', 1, 'varnames', 'PC', 'display', 'off');

% Plot Violins
violins = violinplot(stroop, pc_lab,...
                     'ViolinAlpha', violin_alpha, 'ShowData', true);

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
ax.YLabel.String   = 'I-C RT (z-score)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['Group (n=' num2str(numel(SBJs)) ') Stroop RT Effect by PC: p=' num2str(pval,'%.3f')];
ax.Title.FontSize  = title_sz;
legend([legend_obj{:}],pc_lab,'FontSize',leg_sz,'Location','best');

if save_fig
    fig_fname = [fig_dir fig_name '.' fig_ftype];
    fprintf('Saving %s\n',fig_fname);
    if strcmp(fig_ftype,'eps')
        eval(['export_fig ' fig_fname]);
    else
        saveas(gcf,fig_fname);
    end
end

%% Neutral Trials across PC Blocks
fig_name = 'GRP_RT_violins_PC_N';
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.4 0.5],'Visible',fig_vis);
hold on; ax = gca;

% % Separate data by PC
% n_ix = find(strcmp(cni_lab,'N'));
% rt_grps   = cell(size(pc_lab));
% sbj_n_idx = cell(size(pc_lab));
% pc_n_idx  = cell(size(pc_lab));
% pc_legend = cell(size(pc_lab));
% for pc_ix = 1:numel(pc_lab)
%     for s_ix = 1:numel(SBJs)
%         rt_grps{pc_ix}   = [rt_grps{pc_ix}; rts{s_ix}(pc_idx{s_ix}==pc_ix & cni_idx{s_ix}==n_ix)];
%         sbj_n_idx{pc_ix} = [sbj_n_idx{pc_ix}; sbj_idx{s_ix}(pc_idx{s_ix}==pc_ix & cni_idx{s_ix}==n_ix)];
%         pc_n_idx{pc_ix}  = [pc_n_idx{pc_ix}; pc_idx{s_ix}(pc_idx{s_ix}==pc_ix & cni_idx{s_ix}==n_ix)];
%     end
%     pc_legend{pc_ix} = [pc_lab{pc_ix} ' (n=' num2str(numel(rt_grps{pc_ix})) ')'];
% end
% 
% % Run ANOVA on Neutral Trials
% factor_lab = {'PC','SBJ'};
% [pval, table] = anovan(vertcat(rt_grps{:}),...
%     {vertcat(pc_n_idx{:}), vertcat(sbj_n_idx{:})},...
%     'model', 'interaction', 'random', 2,...%, 'sstype', 3% 'continuous', strmatch('RT',w2.cond),
%     'varnames', factor_lab, 'display', 'off');

% % Create PC groups combining EQ + MI
% n_ix = find(strcmp(cni_lab,'N'));
% eq_ix = find(strcmp(pc_lab,'EQ'));
% mc_ix = find(strcmp(pc_lab,'MC'));
% mi_ix = find(strcmp(pc_lab,'MI'));
% n_pc_lab = [pc_lab(mc_ix) strjoin(pc_lab([eq_ix, mi_ix]),'+')];
% n_pc_colors = [pc_colors(mc_ix) {mean([pc_colors{eq_ix}; pc_colors{mi_ix}],1)}];
% 
% % Separate data by PC, but combine EQ and MI
% rt_grps   = cell(size(n_pc_lab));
% sbj_n_idx = cell(size(n_pc_lab));
% pc_n_idx  = cell(size(n_pc_lab));
% pc_legend = cell(size(n_pc_lab));
% for pc_ix = 1:numel(n_pc_lab)
%     for s_ix = 1:numel(SBJs)
%         if pc_ix==1
%             tmp_pc_idx = pc_idx{s_ix}==mc_ix;
%         else
%             tmp_pc_idx = any([pc_idx{s_ix}==eq_ix pc_idx{s_ix}==mi_ix],2);
%         end
%         rt_grps{pc_ix}   = [rt_grps{pc_ix}; rts{s_ix}(tmp_pc_idx & cni_idx{s_ix}==n_ix)];
%         sbj_n_idx{pc_ix} = [sbj_n_idx{pc_ix}; sbj_idx{s_ix}(tmp_pc_idx & cni_idx{s_ix}==n_ix)];
%         pc_n_idx{pc_ix}  = [pc_n_idx{pc_ix}; pc_idx{s_ix}(tmp_pc_idx & cni_idx{s_ix}==n_ix)];
%     end
%     pc_legend{pc_ix} = [n_pc_lab{pc_ix} ' (n=' num2str(numel(rt_grps{pc_ix})) ')'];
% end

% Compare only MI and MC
n_ix = find(strcmp(cni_lab,'N'));
mc_ix = find(strcmp(pc_lab,'MC'));
mi_ix = find(strcmp(pc_lab,'MI'));
n_pc_lab = pc_lab([mc_ix mi_ix]);
n_pc_colors = pc_colors([mc_ix mi_ix]);

% Separate data by PC, but combine EQ and MI
rt_grps   = cell(size(n_pc_lab));
sbj_n_idx = cell(size(n_pc_lab));
pc_n_idx  = cell(size(n_pc_lab));
pc_legend = cell(size(n_pc_lab));
for pc_ix = 1:numel(n_pc_lab)
    for s_ix = 1:numel(SBJs)
        tmp_pc_idx = pc_idx{s_ix}==find(strcmp(pc_lab,n_pc_lab{pc_ix}));
        rt_grps{pc_ix}   = [rt_grps{pc_ix}; rts{s_ix}(tmp_pc_idx & cni_idx{s_ix}==n_ix)];
        sbj_n_idx{pc_ix} = [sbj_n_idx{pc_ix}; sbj_idx{s_ix}(tmp_pc_idx & cni_idx{s_ix}==n_ix)];
        pc_n_idx{pc_ix}  = [pc_n_idx{pc_ix}; pc_idx{s_ix}(tmp_pc_idx & cni_idx{s_ix}==n_ix)];
    end
    pc_legend{pc_ix} = [n_pc_lab{pc_ix} ' (n=' num2str(numel(rt_grps{pc_ix})) ')'];
end

% Run ANOVA on Neutral Trials
factor_lab = {'PC','SBJ'};
[pval, table] = anovan(vertcat(rt_grps{:}),...
    {vertcat(pc_n_idx{:}), vertcat(sbj_n_idx{:})},...
    'model', 'interaction', 'random', 2,...%, 'sstype', 3% 'continuous', strmatch('RT',w2.cond),
    'varnames', factor_lab, 'display', 'off');

% Plot Violins
violins = violinplot(padcat(rt_grps{:}), n_pc_lab,...
                     'ViolinAlpha', violin_alpha, 'ShowData', false);

% Adjust plot propeties
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
ax.YLabel.String   = 'RT (z-score)';
ax.YLabel.FontSize = axis_sz;
ax.Title.String    = ['Group (n=' num2str(numel(SBJs)) ') RTs: N PC p=' num2str(pval(1),'%.3f')];
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

end
