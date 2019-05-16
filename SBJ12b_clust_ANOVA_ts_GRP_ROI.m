function SBJ12b_clust_ANOVA_ts_GRP_ROI(SBJs,pipeline_id,stat_id,clust_id,an_id,roi_id,atlas_id,plt_id,fig_vis,save_fig)
% Build connectivity matrix based on HFA correlations
%   non-parametric stats via circular shift of trial time series
fig_filetype = 'png';
error('finish this!');
%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/clust_vars/' clust_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
% if rt_correlation
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');
% end
cond_lab = [grp_lab rt_lab];
event_lab = {'stim', 'resp'};

% Get event timing
if strcmp(an_id(1:5),'HGm_S')
    event = 'stim';
elseif strcmp(an_id(1:5),'HGm_R')
    event = 'resp';
end

% Compute mean RT
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
mean_RT = mean(trial_info.response_time);

% Load data
f_name = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id '.mat'];
load(f_name);
% tmp = load(f_name,'hfa'); hfa{1} = tmp.hfa;
% actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
% load(actv_filename,'actv_ch','actv_ch_epochs');
% tmp = load(actv_filename,'hfa'); hfa_actv = tmp.hfa;

sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));

%% Load ROI and GM/WM info
% [roi_list, roi_colors, ~] = fn_roi_label_styles(roi_id);
if strcmp(atlas_id,'Yeo7') || strcmp(atlas_id,'Yeo17')
    elec_space = 'mni_v';
else
    elec_space = 'pat';
end
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_' elec_space '_' atlas_id '.mat'];
load(elec_fname);

% Sort elecs by stat labels
cfgs = []; cfgs.channel = w2.label;
elec = fn_select_elec(cfgs,elec);
if ~all(strcmp(elec.label,w2.label))
    error('need to reorder elec!');
end
elec.roi = fn_atlas2roi_labels(elec.atlas_label,atlas_id,roi_id);
elec.roi_id = roi_id;
elec.roi_color = fn_roi2color(elec.roi);

%% Prep Data
% FDR correct pvalues for ANOVA
for ch_ix = 1:numel(stat.label)
    pvals = squeeze(w2.pval(:,ch_ix,:));
    [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(pvals);%,0.05,'pdep','yes');
end

% Get Sliding Window Parameters
win_lim    = fn_sliding_window_lim(stat.time,round(win_len*sample_rate),round(win_step*sample_rate));
win_center = round(mean(win_lim,2));

% Convert % explained variance to 0-100 scale
w2.trial = w2.trial*100;

% Trim data to plotting epoch
%   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
%   w2 should fit within that since it's averaging into a smaller window
cfg_trim = [];
if strcmp(event,'stim')
    cfg_trim.latency = plt_vars.plt_lim_S;
else
    cfg_trim.latency = plt_vars.plt_lim_R;
end
% hfa{1}  = ft_selectdata(cfg_trim,hfa{1});
stat = ft_selectdata(cfg_trim,stat);

%% Cluster Data
% T = clusterdata(X, CUTOFF);
%   X: MxN matrix, M observations with N variables
%       w2.trial = [factors,chan,time]
%   CUTOFF: int > 2 is max_n_clusters
%   T: Mx1 integer assignment to cluster #
% runs pdist (distance metric), linkage (), cluster ()

roi_list = unique(elec.roi);
if isnumeric(clust_vars.k_method)
    n_clust = clust_vars.k_method;
elseif strcmp(clust_vars.k_method,'roi_match')
    n_clust  = numel(unique(elec.roi));
end
clust_colors = distinguishable_colors(n_clust);
clust_names  = cell([1 n_clust]);
for clust_ix = 1:n_clust
    clust_names{clust_ix} = ['Clust' num2str(clust_ix)];
end
clust_data = cell(size(cond_lab));
clusters = zeros([numel(elec.label) numel(cond_lab)]);
centroids = cell([numel(cond_lab) 1]);
for cond_ix = 1:numel(cond_lab)
    if cond_ix <= numel(grp_lab)
        clust_data{cond_ix} = squeeze(w2.trial(cond_ix,:,:));
    else
        clust_data{cond_ix} = squeeze(stat.rho(:,1,:));
    end
    if strcmp(clust_vars.clust_method,'hier')
        dists = pdist(clust_data{cond_ix}, clust_vars.dist_metric);
        links = linkage(dists, clust_vars.link_method);
        clusters(:,cond_ix) = cluster(links, 'maxclust', n_clust);
    elseif strcmp(clust_vars.clust_method,'kmeans')
        if strcmp(clust_vars.k_method,'CH')
            eva = evalclusters(clust_data{cond_ix},'kmeans','CalinskiHarabasz','KList',[1:numel(unique(elec.roi))]);
            n_clust = eva.OptimalK;
        end
        dist_sums = zeros([numel(cond_lab) n_clust]);
        distances = zeros([numel(cond_lab) numel(elec.label) n_clust]);
        [clusters(:,cond_ix), centroids{cond_ix}, dist_sums(cond_ix,:), distances(cond_ix,:,:)] = ...
            kmeans(clust_data{cond_ix}, n_clust, 'distance', clust_vars.dist_metric, 'replicates', clust_vars.n_iter);
    end
end

out_fname = [SBJ_vars.dirs.proc SBJ '_' clust_id '_' stat_id '_' an_id '.mat'];
if strcmp(clust_vars.k_method,'CH')
    save(out_fname,'-v7.3','eva','clust_data','clusters','clust_colors','clust_names',...
        'centroids','dist_sums','distances','elec');
else
    save(out_fname,'-v7.3','clust_data','clusters','clust_colors','clust_names',...
        'centroids','dist_sums','distances','elec');
end

%% Plot Quality assessment of Clustering
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/clust/' clust_id '/' an_id '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

if strcmp(clust_vars.clust_method,'hier')
% dendrogram plot the tree
elseif strcmp(clust_vars.clust_method,'kmeans')
%     figure;for s=1:6;subplot(2,3,s);histogram(D(:,s));title(['clust ' num2str(s) ': sumD=' num2str(sumD(s))]);end
% figure;for s=1:6;subplot(2,3,s);plot(stat.time,centroids(s,:));title(['clust ' num2str(s) ': sumD=' num2str(sumD(s))]);end
end

% Silhouette plot
fig_name = [SBJ '_ANOVA_clust_' stat_id '_SR_' cond_lab{cond_ix} '_' roi_id '_' atlas_id '_silhouette'];
f = figure('Name',fig_name,'units','normalized',...
    'outerposition',[0 0 1 1],'Visible',fig_vis);
for cond_ix = 1:numel(cond_lab)
    ax = subplot(1,numel(cond_lab),cond_ix);
    [sil,~] = silhouette(clust_data{cond_ix},clusters(:,cond_ix));
    ax.Title.String  = [cond_lab{cond_ix} ' silhouette = ' num2str(mean(sil))];
end
% Save figure
if save_fig
    fig_filename = [fig_dir fig_name '.' fig_filetype];
    fprintf('Saving %s\n',fig_filename);
    saveas(gcf,fig_filename);
    %eval(['export_fig ' fig_filename]);
end

%% Plot ANOVA Time Series by Cluster (colored by ROI)
% Find plot limits
max_w2 = max(max(max(w2.trial)));
min_w2 = min(min(min(w2.trial)));
ylim1_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims1  = [min_w2-ylim1_fudge max_w2+ylim1_fudge];
yticks1 = 0:1:ylims1(2);

max_rho = max(max(squeeze(stat.rho)));
min_rho = min(min(squeeze(stat.rho)));
ylims2  = [round(min_rho*10)/10-0.1 round(max_rho*10)/10+0.1]; % extra on top and bottom for StdErr
yticks2 = ylims2(1):0.1:ylims2(2);
y_sig = zeros([1 numel(grp_lab)+1]);
y_sig(1) = mean([min_w2,max_w2]);
for grp_ix = 2:numel(grp_lab)+2
    y_sig(grp_ix) = y_sig(grp_ix-1)+ylim1_fudge;
end

%!!! SR plotting or separate? maybe need to combine them on one plot axis...
for cond_ix = 1:numel(cond_lab)
    fig_name = [SBJ '_clust_' stat_id '_' an_id '_' cond_lab{cond_ix} '_' roi_id '_' atlas_id];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 1],'Visible',fig_vis);
    [subplot_lay, ~] = fn_num_subplots(n_clust);
    for clust_ix = 1:n_clust
        ax = subplot(subplot_lay(1), subplot_lay(2), clust_ix);
        hold on;
        
        clust_elecs = find(clusters(:,cond_ix)==clust_ix);
        for ch = 1:numel(clust_elecs)
            ch_ix = clust_elecs(ch);
            %         roi_ix = find(strcmp(roi_list,elec.roi{ch_ix}));
            if cond_ix <= numel(grp_lab)
                ylims = ylims1;
                yticks = yticks1;
                ylab = '% Variance Explained';
                plot(win_center,squeeze(w2.trial(cond_ix,ch_ix,:))',...
                    'Color',elec.roi_color(ch_ix,:),'LineStyle',grp_style{cond_ix});
                % Find significant periods
                if strcmp(plt_vars.sig_type,'bold')
                    sig_chunks = fn_find_chunks(squeeze(qvals(cond_ix,ch_ix,:))<0.05);
                    sig_chunks(squeeze(qvals(cond_ix,ch_ix,sig_chunks(:,1)))>0.05,:) = [];
                    for sig_ix = 1:size(sig_chunks,1)
                        line(win_center(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                            squeeze(w2.trial(cond_ix,ch_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                            'Color',elec.roi_color(ch_ix,:),'LineStyle',plt_vars.sig_style,...
                            'LineWidth',plt_vars.sig_width);
                    end
                else
                    error('Plot the bold version of significance!');
                end
%                 % Plot centroid
%                 plot(win_center,centroids{cond_ix}(clust_ix,:),...
%                     'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
            else
                % RT correlation significant time periods
                ylims = ylims2;
                yticks = yticks2;
                ylab = 'Correlation with RT';
                plot(1:numel(stat.time),squeeze(stat.rho(ch_ix,:,:))',...
                    'Color',elec.roi_color(ch_ix,:),'LineStyle','-');
                
                sig_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,:,:)));
                sig_chunks(squeeze(stat.mask(ch_ix,:,sig_chunks(:,1)))==0,:) = [];
                if strcmp(plt_vars.sig_type,'bold')
                    for sig_ix = 1:size(sig_chunks,1)
                        sig_times = sig_chunks(sig_ix,1):sig_chunks(sig_ix,2);
                        line(sig_times,squeeze(stat.rho(ch_ix,:,sig_times)),...
                            'Color',elec.roi_color(ch_ix,:),'LineStyle',plt_vars.sig_style,...
                            'LineWidth',plt_vars.sig_width);
                    end
                else
                    error('plot the bold line version of sig!');
                end
%                 % Plot centroid
%                 plot(1:numel(stat.time),centroids{cond_ix}(clust_ix,:),...
%                     'Color',[0.4 0.4 0.4],'LineStyle','--','LineWidth',2);
            end
        end
                
        % Plot event
        if strcmp(event_lab,'stim')
            x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            mean_RT_pre = find(stat.time<=mean_RT);
            event_line = line([mean_RT_pre(end) mean_RT_pre(end)],ylims,...
                'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
        else
            x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Plot Response Marker
            event_line = line([find(stat.time==0) find(stat.time==0)],ylims,...
                'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
            %                 main_lines = [main_lines event_line];
            %                 lgd_lab = {lgd_lab{:} 'RT'};
        end
        
        % Plotting parameters
        ax.Title.String  = ['Clust ' num2str(clust_ix) ': ' event_lab];
        ax.Box           = 'off';
        ax.YLim          = ylims;
        ax.YTick         = yticks;
        ax.YLabel.String = ylab;
        ax.XLim          = [0,size(stat.time,2)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(stat.time,2);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
    end
    % Save figure
    if save_fig
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
    
    %% Plot cluster centroids
    if strcmp(clust_vars.clust_method,'kmeans')
        fig_name = [SBJ '_clust_' stat_id '_' an_id '_' cond_lab{cond_ix} '_' roi_id '_' atlas_id '_centroids'];
        f = figure('Name',fig_name,'units','normalized',...
            'outerposition',[0 0 0.5 0.5],'Visible',fig_vis);
        
        % Plot centroids
        hold on;
        for clust_ix = 1:n_clust
            if cond_ix<=numel(grp_lab)
                plot(win_center,centroids{cond_ix}(clust_ix,:),...
                    'Color',clust_colors(clust_ix,:),'LineStyle','-','LineWidth',2);
            else
                plot(1:numel(stat.time),centroids{cond_ix}(clust_ix,:),...
                    'Color',clust_colors(clust_ix,:),'LineStyle','-','LineWidth',2);
            end
        end
        legend(clust_names);
        
        % Plot event
        if strcmp(event_lab,'stim')
            x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
            mean_RT_pre = find(stat.time<=mean_RT);
            event_line = line([mean_RT_pre(end) mean_RT_pre(end)],ylim,...
                'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
        else
            x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Plot Response Marker
            event_line = line([find(stat.time==0) find(stat.time==0)],ylim,...
                'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
            %                 main_lines = [main_lines event_line];
            %                 lgd_lab = {lgd_lab{:} 'RT'};
        end
        
        % Plotting parameters
        ax.Title.String  = ['Cluster Centroids: ' event_lab];
        ax.Box           = 'off';
%         ax.YLim          = ylims;
%         ax.YTick         = yticks;
%         ax.YLabel.String = ylab;
        ax.XLim          = [0,size(stat.time,2)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(stat.time,2);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
        
        % Save figure
        if save_fig
            fig_filename = [fig_dir fig_name '.' fig_filetype];
            fprintf('Saving %s\n',fig_filename);
            saveas(gcf,fig_filename);
            %eval(['export_fig ' fig_filename]);
        end
    end
end

%% Print clustering composition report
report_fname = [SBJ_vars.dirs.proc SBJ '_clust_' stat_id '_' an_id '_' cond_lab{cond_ix} '_' roi_id '_' atlas_id '_composition.txt'];
report = fopen(report_fname,'w');
for cond_ix = 1:numel(cond_lab)
    fprintf(report,'=============================================================\n');
    fprintf(report,'\t%s\n',cond_lab{cond_ix});
    %     comp_str = ['Cluster %i:' repmat('\t%.02f %s,',[1 numel(clust_n)]) '\n'];
    for clust_ix = 1:numel(n_clust)
        cfgs = [];
        cfgs.channel = elec.label(clusters(:,cond_ix)==clust_ix);
        elec_clust = fn_select_elec(cfgs, elec);
        
        fprintf(report,'Cluster #%i:',clust_ix);
        for roi_ix = 1:numel(roi_list)
            fprintf(report,'\t%.02f %s,',...
                sum(strcmp(elec_clust.roi,roi_list{roi_ix}))/numel(elec_clust.label),...
                roi_list{roi_ix});
        end
        fprintf(report,'\n');
    end
end
fclose(report);


end



%% ========================================================================
%% Plot cluster means vs. centroids
%%=========================================================================
% for cond_ix = 1:numel(cond_lab)
%     fig_name = [SBJ '_ANOVA_clust_' stat_id '_SR_' cond_lab{cond_ix} '_' roi_id '_' atlas_id];
%     f = figure('Name',fig_name,'units','normalized');
%     ax = subplot(2, 1, 1);
%     hold on;
%     for clust_ix = 1:n_clust
%         clust_elecs = find(clusters(:,cond_ix)==clust_ix);
%         if cond_ix <= numel(grp_lab)
%             ylab = '% Variance Explained';
%             data = squeeze(mean(w2.trial(cond_ix,clust_elecs,:),2));
%             plot(win_center,(data-mean(data))/std(data),'Color',roi_colors(clust_ix,:));
%         else
%             % RT correlation significant time periods
%             data = squeeze(mean(stat.rho(clust_elecs,:,:),1));
%             ylab = 'Correlation with RT';
%             plot(1:numel(stat.time),(data-mean(data))/std(data),'Color',roi_colors(clust_ix,:));
%         end
%     end
%     
%     % Plot event
%     if strcmp(event_lab,'stim')
%         x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
%         mean_RT_pre = find(stat.time<=mean_RT);
%         event_line = line([mean_RT_pre(end) mean_RT_pre(end)],ylims,...
%             'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%             'LineStyle',plt_vars.evnt_style);
%     else
%         x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
%         % Plot Response Marker
%         event_line = line([find(stat.time==0) find(stat.time==0)],ylims,...
%             'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%             'LineStyle',plt_vars.evnt_style);
%     end
%     
%     % Plotting parameters
%     ax.Title.String  = ['z-scored Cluster Means: ' event_lab];
%     ax.Box           = 'off';
%     ax.YLabel.String = ylab;
%     ax.XLim          = [0,size(stat.time,2)];
%     ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(stat.time,2);
%     ax.XTickLabel    = x_tick_lab;
%     ax.XLabel.String = 'Time (s)';
%     
%     %% Centroids
%     ax = subplot(2, 1, 2);
%     hold on;
%     for clust_ix = 1:n_clust
%         if cond_ix <= numel(grp_lab)
%             ylab = '% Variance Explained';
%             plot(win_center,centroids{cond_ix}(clust_ix,:),'Color',roi_colors(clust_ix,:));
%         else
%             % RT correlation significant time periods
%             ylab = 'Correlation with RT';
%             plot(1:numel(stat.time),centroids{cond_ix}(clust_ix,:),'Color',roi_colors(clust_ix,:));
%         end
%     end
%     
%     % Plot event
%     if strcmp(event_lab,'stim')
%         x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
%         mean_RT_pre = find(stat.time<=mean_RT);
%         event_line = line([mean_RT_pre(end) mean_RT_pre(end)],ylim,...
%             'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%             'LineStyle',plt_vars.evnt_style);
%     else
%         x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
%         % Plot Response Marker
%         event_line = line([find(stat.time==0) find(stat.time==0)],ylim,...
%             'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
%             'LineStyle',plt_vars.evnt_style);
%     end
%     
%     % Plotting parameters
%     ax.Title.String  = ['Cluster Centroids: ' event_lab];
%     ax.Box           = 'off';
%     ax.YLabel.String = ylab;
%     ax.XLim          = [0,size(stat.time,2)];
%     ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(stat.time,2);
%     ax.XTickLabel    = x_tick_lab;
%     ax.XLabel.String = 'Time (s)';
% end
