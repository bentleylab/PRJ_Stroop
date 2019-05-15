function SBJ12b_clust_ANOVA_ts_GRP_peaks(SBJs,pipeline_id,stat_id,clust_id,an_id,...
    roi_id,atlas_id,plt_id,fig_vis,save_fig,plot_out)
% Sort the kmeans cluster centroids by their peak time, then combine across
%   patients by peak time calssification, then plot group recon with elecs
%   colored by peak time
fig_filetype = 'svg';

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Get parameters
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
% if rt_correlation
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');
% end
cond_lab = [grp_lab rt_lab];

% Get event timing
if strcmp(an_id(1:5),'HGm_S')
    event_lab = 'stim';
    peak_bins = [0.3 0.6 1 2];
elseif strcmp(an_id(1:5),'HGm_R')
    event_lab = 'resp';
    peak_bins = [7 14 21 40; 7 14 21 40; 400 800 1200 2000];
end

% ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if strcmp(atlas_id,'Yeo7') || strcmp(atlas_id,'Yeo17')
    elec_space = 'mni_v';
else
    elec_space = 'pat';
end

% Cluster info
clust_names  = cell([1 numel(roi_list)]);
for clust_ix = 1:numel(roi_list)
    clust_names{clust_ix} = ['C' num2str(clust_ix)];
end
clust_colors = [0 0 1; 0 1 1; 1 0 1; 1 0 0];  %garish blue, cyan, magenta, red
% clust_colors = [161, 218, 180; 65, 182, 196; 44 127 184; 37 52 148]./255;
% clust_fname = [root_dir 'PRJ_Stroop/data/' SBJs{1} '/04_proc/' SBJs{1} '_' clust_id...
%     '_' stat_id '_' an_id '_' atlas_id '_' roi_id '.mat'];
% tmp_cent = load(clust_fname,'centroids');
% clust_len = size(tmp_cent.centroids{1},2);

%% Load Results
elec       = cell([numel(SBJs) 1]);
clusters   = cell([numel(SBJs) 1]);
clust_data = cell([numel(SBJs) numel(cond_lab)]);
centroids  = cell([numel(SBJs) numel(cond_lab)]);
cent_bins  = cell([numel(SBJs) numel(cond_lab)]);
peak_time  = cell([numel(SBJs) numel(cond_lab)]);
clust_bin = cell([numel(SBJs) numel(cond_lab)]);
% mean_RTs  = zeros(size(SBJs));
for sbj_ix = 1:numel(SBJs)
    %% Load data
    % SBJ processing
    SBJ = SBJs{sbj_ix};
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % one stat to get timing
    if sbj_ix==1
        stat_fname = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
        load(stat_fname,'stat');
        sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));
        % Get Sliding Window Parameters
        win_lim    = fn_sliding_window_lim(stat.time,win_len*sample_rate,win_step*sample_rate);
        win_center = round(mean(win_lim,2));
        
        % Trim data to plotting epoch
        %   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
        %   w2 should fit within that since it's averaging into a smaller window
        cfg_trim = [];
        if strcmp(event_lab,'stim')
            cfg_trim.latency = plt_vars.plt_lim_S;
        else
            cfg_trim.latency = plt_vars.plt_lim_R;
        end
        stat = ft_selectdata(cfg_trim,stat);
    end
    
    % Compute mean RT
    if strcmp(event_lab,'stim')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
        % Compute mean RT
        mean_RTs(sbj_ix) = mean(trial_info.response_time);
    end
    
    % Load clustering data
    % contains: 'clust_data','clusters','clust_colors','clust_names','centroids','dist_sums','distances','elec'
    %   and maybe 'eva' if nCH criterion
    clust_fname = [SBJ_vars.dirs.proc SBJ '_' clust_id '_' stat_id '_' an_id '_' atlas_id '_' roi_id '.mat'];
    tmp = load(clust_fname);
    elec{sbj_ix} = tmp.elec;
    
    for cond_ix = 1:numel(cond_lab)
        clust_data{sbj_ix,cond_ix} = tmp.clust_data{cond_ix};
        clusters{sbj_ix} = tmp.clusters;
        centroids{sbj_ix,cond_ix} = tmp.centroids{cond_ix};
        
        %% Sort clusters by peak time
        [~,peak_time] = max(centroids{sbj_ix,cond_ix},[],2);
        [~,clust_bin{sbj_ix,cond_ix}] = histc(peak_time,peak_bins(cond_ix,:));
%         [~,peak_order{sbj_ix,cond_ix}] = sort(peak_time{sbj_ix,cond_ix});
    end
    clear SBJ SBJ_vars_cmd SBJ_vars tmp
end

%% Plot centroids across SBJ
% all_times{cond_ix} = vertcat(peak_time{:,1});
% Assign centroids to bins
for cond_ix = 1%:numel(cond_lab)
    % Plot all ANOVA ts
    figure('Name',cond_lab{cond_ix},'Visible',fig_vis);
    for sbj_ix = 1:numel(SBJs)
        %         elec{sbj_ix}.clusters = clusters{sbj_ix};
        %         elec{sbj_ix}.clust_bin = [clust_bin{sbj_ix,:}];
        %         elec{sbj_ix}.cond_lab = cond_lab;
        for ch_ix = 1:size(clusters{sbj_ix},1)
            clust_n = clusters{sbj_ix}(ch_ix,cond_ix);
            time_bin = clust_bin{sbj_ix,cond_ix}(clust_n)+1;
            subplot(4,1,time_bin); hold on;
            plot(clust_data{sbj_ix,cond_ix}(ch_ix,:),'Color',elec{sbj_ix}.roi_color(ch_ix,:));%clust_colors(time_bin,:));
        end
        %         [bin_assign,~] = hist(peak_time
        %         all_cent{cond_ix} = vertcat(centroids{:,cond_ix});
    end
    
    % Plot all centroids
    ylims = [-0.4 0.6]; %just visually took this from CNI plots
    yticks = -0.4:0.2:0.6;
    fig_name = ['GRP_' cond_lab{cond_ix} '_' atlas_id '_' roi_id '_centroids'];
    figure('Name',fig_name,'Visible',fig_vis);
    for sbj_ix = 1:numel(SBJs)
        for clust_ix = 1:size(centroids{sbj_ix,cond_ix},1)
            time_bin = clust_bin{sbj_ix,cond_ix}(clust_ix)+1;
            subplot(4,1,time_bin); hold on;
            plot(win_center,centroids{sbj_ix,cond_ix}(clust_ix,:),'Color',clust_colors(time_bin,:));
        end
        %         [bin_assign,~] = hist(peak_time
        %         all_cent{cond_ix} = vertcat(centroids{:,cond_ix});
    end
    
    for clust_ix = 1:numel(roi_list)
        subplot(4,1,clust_ix);
        % Plot event
        if strcmp(event_lab,'stim')
            x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            mean_RT_pre = find(stat.time<=mean_RT);
            event_line = line([mean_RT_pre(end) mean_RT_pre(end)],ylims,...
                'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
        else
            x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            event_line = line([find(stat.time==0) find(stat.time==0)],ylims,...
                'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
        end
        
        % Plotting parameters
        %             ax(bin_ix).Box           = 'off';
        set(gca,'ylim',ylims);
        set(gca,'YTick',yticks);
        %             set(gca,'YLabel',ylab);
        set(gca,'XLim',[0,size(stat.time,2)]);
        set(gca,'XTick',0:plt_vars.x_step_sz*sample_rate:size(stat.time,2));
        set(gca,'XTickLabel',x_tick_lab);
        %             set(gca,'XLabel','Time (s)');
    end
    
    if save_fig
        fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_clust/' clust_id '/' stat_id '/' an_id '/'];
        fig_fname = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
    end
    
end

%% Save time bins as new cluster identity


end
