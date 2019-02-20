% function [ output_args ] = untitled3(SBJ,elec_lab,conditions,an_id,stat_id,plt_id,save_fig,fig_vis)
%Plot HFA time series averaged by RT splits
%   Detailed explanation goes here

SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};
conditions = 'CI';
% an_id = 'HGm_S_zbtS_trl2to251_sm0_wn100_stat0';
% elec_lab = {'LAC2-3'};
pipeline_id = 'main_ft';
save_fig = 0;
fig_vis = 'on';
% plt_id = 'ts_errbr_evnt';

RT_splits = [2 3 4];    % 2 = median, 3 = thirds, 4 = quartiles
split_mrk = {'x','o','p','d'};

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(conditions);
SBJ_colors = distinguishable_colors(numel(SBJs));

%% Compute splits
split_cnt = cell([numel(RT_splits) numel(SBJs)]);
split_avg = cell([numel(RT_splits) numel(SBJs)]);
split_cut = cell([numel(RT_splits) numel(SBJs)]);

for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Load RTs
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    
    % Compute mean RT per condition
    RTs = round(1000*trial_info.response_time); % converts sec to ms
    cond_idx = zeros(size(RTs));
    rt_cond = cell(size(cond_lab));
    for cond_ix = 1:numel(cond_lab)
        cond_idx(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1) = cond_ix;
        rt_cond{cond_ix} = RTs(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1);
    end
    
    for split_ix = 1:numel(RT_splits)
        %% Split trials by RT
        split_cnt{split_ix,sbj_ix} = zeros([numel(cond_lab) RT_splits(split_ix)+1]);
        split_avg{split_ix,sbj_ix} = zeros([numel(cond_lab) RT_splits(split_ix)]);
        split_idx = cell(size(cond_lab));
        
        split_cut{split_ix,sbj_ix} = quantile(RTs,0:1/RT_splits(split_ix):1);
        split_cut{split_ix,sbj_ix}(end) = split_cut{split_ix,sbj_ix}(end)+1;  % make sure longest RT gets included in last bin (histc will put an equal value in singleton bin)
        for cond_ix = 1:numel(cond_lab)
            [split_cnt{split_ix,sbj_ix}(cond_ix,:), split_idx{cond_ix}] = histc(rt_cond{cond_ix},split_cut{split_ix,sbj_ix});
            if split_cnt{split_ix,sbj_ix}(cond_ix,end)~=0
                error('histc messed up and put something below the lowest RT!');
            end
            for split_n = 1:RT_splits(split_ix)
                split_avg{split_ix}(cond_ix,split_n) = mean(rt_cond{cond_ix}(split_idx{cond_ix}==split_n));
            end
        end
        split_cnt{split_ix,sbj_ix}(:,end) = [];
        
    end
    clear SBJ SBJ_vars SBJ_vars_cmd split_idx RTs trial_info cond_idx rt_cond
end
%% Plot Results
% fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/' conditions '/RT_split' num2str(RT_split) '/'];
% if ~exist(fig_dir,'dir')
%     mkdir(fig_dir);
% end

figure;
rt_lims = [100 85 70];
for split_ix = 1:numel(RT_splits)
    subplot(1,numel(RT_splits),split_ix); hold on;
%     cond_counts = [];
%     for cond_ix = 1:numel(cond_lab)
%         cond_counts = strcat(cont_counts, num2str(split_cnt(cond_ix,split_ix)));
%     end
%     title(['rts < ' num2str(split_cut{split_ix,sbj_ix}(split_ix)) ', ' conditions ' = ' cond_counts]);
    for sbj_ix = 1:numel(SBJs)
        for split_n = 1:RT_splits(split_ix)
            scatter(split_cnt{split_ix,sbj_ix}(1,split_n),split_cnt{split_ix,sbj_ix}(2,split_n),...
                split_mrk{split_n},'MarkerEdgeColor',SBJ_colors(sbj_ix,:));
        end
    end
    line([0 100],[0 100],'Color','k','LineStyle',':');
    xlabel([cond_lab{1} ' RT count']);
    ylabel([cond_lab{2} ' RT count']);
    xlim([0 100]);%rt_lims(split_ix)]);
    ylim([0 100]);%rt_lims(split_ix)]);
    title(RT_splits(split_ix));
end

%% Line plot
figure;
rt_lims = [100 85 70];
for split_ix = 1:numel(RT_splits)
    subplot(1,numel(RT_splits),split_ix); hold on;
%     cond_counts = [];
%     for cond_ix = 1:numel(cond_lab)
%         cond_counts = strcat(cont_counts, num2str(split_cnt(cond_ix,split_ix)));
%     end
%     title(['rts < ' num2str(split_cut{split_ix,sbj_ix}(split_ix)) ', ' conditions ' = ' cond_counts]);
    for sbj_ix = 1:numel(SBJs)
%         for split_n = 1:RT_splits(split_ix)
            line(split_cnt{split_ix,sbj_ix}(1,:),split_cnt{split_ix,sbj_ix}(2,:),...
                'Color',SBJ_colors(sbj_ix,:));
%         end
    end
    legend(SBJs)
    line([0 100],[0 100],'Color','k','LineStyle',':');
    line([30 100],[30 30],'Color','r','LineWidth',2,'LineStyle','--');
    line([30 30],[30 100],'Color','r','LineWidth',2,'LineStyle','--');
    xlabel([cond_lab{1} ' RT count']);
    ylabel([cond_lab{2} ' RT count']);
    xlim([0 100]);%rt_lims(split_ix)]);
    ylim([0 100]);%rt_lims(split_ix)]);
    title(RT_splits(split_ix));
end

%% Save figure
% if save_fig
%     fig_filename = [fig_dir fig_name '.' fig_filetype];
%     fprintf('Saving %s\n',fig_filename,'svg');
%     saveas(gcf,fig_filename);
%     %eval(['export_fig ' fig_filename]);
% end

