function B00_RT_GRP_hist_norm(SBJs,conditions,save_fig,fig_vis,fig_type)
%% RT Behavioral analysis- draw box plots per SBJ and for Group
% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
%addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));

% Analysis Parameters
% conditions = 'CNI';
% SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};

late_RT_cut   = 0.6;      % window in sec before next stim to eliminate
% Plotting parameters
hist_alpha    = 0.5;
bin_size      = 0.1;      % in normalized z-scores
line_w        = 2;
[cond_lab, cond_colors, cond_styles] = fn_condition_label_styles(conditions);

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
    late_ix = [];
    for t_ix = 2:length(trial_info.resp_onset)
        if trial_info.word_onset(t_ix)-late_RT_cut*trial_info.sample_rate <= trial_info.resp_onset(t_ix-1)
            late_RT_count = late_RT_count + 1;
            late_ix = [late_ix trial_info.trial_n];
        end
    end
    if late_RT_count>0
        fprintf('%s: %i late trials detected- %i\n',SBJ,late_RT_count,late_ix);
    end
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

%% Normalize RTs
RTs_all  = vertcat(RTs{:});
RTs_norm = cell([numel(SBJs) numel(cond_lab)]);
for sbj_ix = 1:numel(SBJs)
    for cond_ix = 1:numel(cond_lab)
        RTs_norm{sbj_ix, cond_ix} = (RTs{sbj_ix, cond_ix}-mean(RTs_all))/std(RTs_all);
    end
end

%% Histograms per condition
% Trial Type
fig_name = strcat('GRP_RT_hist_norm_',conditions);
figure('Name',fig_name,'units','normalized','outerposition',[0 0 0.3 0.4],'Visible',fig_vis);
hold on;
ax = gca;

% Plot Normalized Histograms
bins = [floor(min(vertcat(RTs_norm{:}))):bin_size:ceil(max(vertcat(RTs_norm{:})))];
cond_alphas = repmat(hist_alpha,size(cond_lab));
if strmatch('neu',cond_lab)
    cond_alphas(strmatch('neu',cond_lab)) = 0.3;
end
for cond_ix = 1:numel(cond_lab)
    histogram(vertcat(RTs_norm{:,cond_ix}),bins,'FaceColor',cond_colors(cond_ix),...
        'FaceAlpha',cond_alphas(cond_ix),'Normalization','probability');
end
% Plot Means
for cond_ix = numel(cond_lab):-1:1 % Reverese order because neutral and con overlap, and con is more important
    line([median(RT_means(:,cond_ix)) median(RT_means(:,cond_ix))], ax.YLim,...
        'Color', cond_colors{cond_ix}, 'LineWidth', line_w);
end
% [~,pval] = ttest2([trial_RTs{1}],[trial_RTs{3}]);%,'Alpha',alpha
ax.XLabel.String   = 'RT (z-score)';
ax.XLabel.FontSize = 14;
ax.YLabel.String   = 'Proportion of RTs';
ax.YLabel.FontSize = 14;
ax.Title.String    = 'Normalized Group Reaction Times';
ax.Title.FontSize  = 16;
% trial_legend = {};
% for cond_ix = 1:length(trial_lab)
%     trial_legend{cond_ix} = [trial_lab{cond_ix} '-' num2str(length(trial_RTs{cond_ix}))];
% end
legend(cond_lab);

%% Save Figure
fig_filename = [fig_dir,fig_name,'.',fig_type];
if save_fig ==1
    fprintf('Saving %s\n',fig_filename);
    if strcmp(fig_type,'eps')
        eval(['export_fig ' fig_filename]);
    else
        saveas(gcf,fig_filename);
    end
end

end
