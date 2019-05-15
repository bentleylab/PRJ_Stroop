function SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA_clustBin(SBJs,clust_id,stat_id,pipeline_id,an_id,roi_id,...
                                                    atlas_id,gm_thresh,plt_id,save_fig,fig_vis,fig_filetype)
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot a summary of those time period per electrode
% Normalize all onset times by mean(RT)
% INPUTS:
%   plt_vars.grp_metric [str] - {'avg','mdn','all'}
%       mean/median will compute that metric within each SBJ (variance is across SBJs)
%       all- all electrode onsets are aggregated as if from the same SBJ
% clear all; %close all;
% fig_filetype = 'png';
label_spacer = 0;
groi_label_spacer = '      ';
if ischar(save_fig); save_fig = str2num(save_fig); end
% if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);

%% Prep variables
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/clust_vars/' clust_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
cond_lab = [grp_lab, {'corr(RT)'}];

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(atlas_id,{'DK','Dx'}))
    view_space = 'pat';
elseif any(strcmp(atlas_id,{'Yeo7','Yeo17'}))
    view_space = 'mni_v';
else
    error(['Unknown atlas_id: ' atlas_id]);
end
if strcmp(plt_vars.violin_scat_colors,'SBJ')
    SBJ_colors = distinguishable_colors(numel(SBJs));
end

% Get event timing
mean_RTs = zeros(size(SBJs));
if strcmp(an_id(1:5),'HGm_S')
    error('tbin for stim locked needs more thinking!');
    event_lab = 'stim';
%     peak_bins = [0.3 0.6 1 2];
elseif strcmp(an_id(1:5),'HGm_R')
    event_lab = 'resp';
%     peak_bins = [7 14 21 40; 7 14 21 40; 400 800 1200 2000];
end
if strcmp(clust_vars.k_method,'roi_match')
    n_tbins = numel(roi_list);
else
    error('n_tbins undefined for non roi_match');
end

% Cluster info
clust_names  = cell([1 numel(roi_list)]);
for clust_ix = 1:numel(roi_list)
    clust_names{clust_ix} = ['C' num2str(clust_ix)];
end
bin_colors = {[0 0 1], [0 1 1], [1 0 1], [1 0 0]};  %garish blue, cyan, magenta, red
% bin_colors = {[27 158 119]./255, [117 112 179]./255, [217 95 2]./255, [231 41 138]./255}; %qualitative max diff from gROI, cool to hot

%% Load Results
% Set up onset counts
all_onsets  = cell([numel(SBJs) n_tbins numel(cond_lab)]);
sig_ch = cell([numel(SBJs) numel(cond_lab)]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    if strcmp(event_lab,'stim')
        load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
        % Compute mean RT
        mean_RTs(sbj_ix) = mean(trial_info.response_time);
    end
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_ANOVA_ROI_',stat_id,'_',an_id,'.mat'));
    
    %% Process parameters
    sample_rate = (numel(stat.time)-1)/(stat.time(end)-stat.time(1));
    
    % FDR correct pvalues for ANOVA
    qvals = NaN(size(w2.pval));
    for ch_ix = 1:numel(stat.label)
        [~, ~, ~, qvals(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));%,0.05,'pdep','yes');
    end
    
    % Get Time Bin and Sliding Window Parameters
    win_lim    = fn_sliding_window_lim(stat.time,win_len*sample_rate,win_step*sample_rate);
    win_center = round(mean(win_lim,2));
    % Assume k_method = roi_match for now...
%     if strcmp(tbin_id,'eqROI')  %!!! check for 1st 2 letters = 'eq'
    % 4 ROIs = R time bins: -0.5, -0.1, 0.25, 0.6, 1.0
    %   peak_bins = [7 14 21 40; 7 14 21 40; 400 800 1200 2000];
    % 4 ROIs = R time bins: [0.3 0.6 1 2]
    peak_bins = zeros([numel(cond_lab) n_tbins]);
    for cond_ix = 1:numel(cond_lab)
        if strcmp(cond_lab{cond_ix},'corr(RT)')
            n_time = numel(stat.time);
        else
            n_time = numel(win_center);
        end
%             bin_step = ceil(n_time/n_tbins);% bin_step:bin_step:numel(win_center)
        edges = linspace(0,n_time,n_tbins+1);   % n_tbins+1 to drop 0
        peak_bins(cond_ix,:) = edges(2:end);
        peak_bins(cond_ix,end) = peak_bins(cond_ix,end)+1;  % histc will assign data on the edge to another bin
    end
%     else
%         error(['Unknown tbin_id: ' tbin_id]);
%     end
    
    %% Load clusters
    clust_fname = [SBJ_vars.dirs.proc SBJ '_' clust_id '_' stat_id '_' an_id '_' atlas_id '_' roi_id '.mat'];
    tmp = load(clust_fname);
    centroids = tmp.centroids; clusters = tmp.clusters;
    clust_elec = tmp.elec;
    
    % Sort clusters by peak time
    clust_bin = cell([numel(cond_lab) 1]);
    for cond_ix = 1:numel(cond_lab)
        [~,peak_time] = max(centroids{cond_ix},[],2);
        [~,clust_bin{cond_ix}] = histc(peak_time,peak_bins(cond_ix,:));
    end
    
    %% Load ROI and GM/WM info
    if any(strcmp(atlas_id,{'DK','Dx'})) % these have tissue probabilities
%         elec_tis_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_' view_space '_' atlas_id '_tis.mat'];
        elec_tis_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_' view_space '_' atlas_id '.mat'];
        load(elec_tis_fname);
    else
        % Load tissue prob
        elec_tis_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_pat_Dx_tis.mat'];
        load(elec_tis_fname);
        tiss_prob = elec.tissue_prob;
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_' view_space '_' atlas_id '.mat'];
        load(elec_fname);
        elec.tissue_prob = tiss_prob;
    end
    
    % Sort elecs by stat labels
    cfgs = []; cfgs.channel = stat.label;
    elec    = fn_select_elec(cfgs,elec);
    elec.roi = fn_atlas2roi_labels(elec.atlas_label,atlas_id,roi_id);
    
%     % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
%     gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
        
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        % If elec matches roi_list and is in GM, get stats
        if any(strcmp(elec.roi{ch_ix},roi_list))% && gm_bin(ch_ix)
            clust_elec_ix = strcmp(elec.label{ch_ix},clust_elec.label);
            if ~any(clust_elec_ix)
                error('trying to process elec that isnt in clusters');
            end
            
%             roi_ix = find(strcmp(elec.roi{ch_ix},roi_list));
            % Get ANOVA group onsets
            for cond_ix = 1:numel(cond_lab)
                % Assign cluster bin
                clust_n = clusters(clust_elec_ix,cond_ix);
                time_bin = clust_bin{cond_ix}(clust_n)+1;

                % ANOVA conditions
                if any(strcmp(cond_lab{cond_ix},grp_lab)) && any(squeeze(qvals(cond_ix,ch_ix,:))<0.05)
                    sig_ch{sbj_ix,cond_ix} = [sig_ch{sbj_ix,cond_ix} stat.label{ch_ix}];
                    sig_onsets = stat.time(win_lim(squeeze(qvals(cond_ix,ch_ix,:))<0.05,1));
                    if strcmp(event_lab,'resp')
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} sig_onsets(1)];
                    elseif strcmp(event_lab,'stim') && (sig_onsets(1)<mean_RTs(sbj_ix))
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} sig_onsets(1)];
                    end
                    
                % Get RT correlation onset
                elseif strcmp(cond_lab{cond_ix},'corr(RT)') && sum(squeeze(stat.mask(ch_ix,1,:)))>0
                    sig_ch{sbj_ix,cond_ix} = [sig_ch{sbj_ix,cond_ix} stat.label{ch_ix}];
                    mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
                    mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
                    % Convert the first onset of significance to time
                    onset_time = stat.time(mask_chunks(1,1));
                    % Exclude differences after the mean RT for this SBJ
                    if strcmp(event_lab,'resp')
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} onset_time];
                    elseif strcmp(event_lab,'stim') && (onset_time<mean_RTs(sbj_ix))
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} onset_time];
                    end
                end
            end
        end
    end
    
    % Normalize all onset times by mean reaction time
    if strcmp(event_lab,'stim')
        for time_bin = 1:n_tbins
            for cond_ix = 1:numel(cond_lab)
                all_onsets{sbj_ix,time_bin,cond_ix} = all_onsets{sbj_ix,time_bin,cond_ix}./mean_RTs(sbj_ix);
            end
        end
    end
    clear SBJ SBJ_vars hfa stat einfo w2
end

%% Aggregate/Process onsets per gROI
% Format as struct to fit violinplot
plot_onsets    = cell([numel(cond_lab) 1]);
plot_onset_sbj = cell([numel(cond_lab) 1]);
good_tbin_map   = cell([numel(cond_lab) 1]);
for cond_ix = 1:numel(cond_lab)
    plot_onsets{cond_ix}  = {};
    plot_onset_sbj{cond_ix}  = {};
    good_tbin_map{cond_ix} = zeros([1 n_tbins]);
    for time_bin = 1:n_tbins
        if strcmp(plt_vars.grp_metric,'all')
            plot_onsets{cond_ix}{time_bin} = [all_onsets{:,time_bin,cond_ix}]';
            plot_onset_sbj{cond_ix}{time_bin} = [];
        else
            plot_onsets{cond_ix}{time_bin} = [];
        end
        for sbj_ix = 1:numel(SBJs)
            % Aggregate onsets per ROI within each SBJ
            if strcmp(plt_vars.grp_metric,'all')
                plot_onset_sbj{cond_ix}{time_bin} = [plot_onset_sbj{cond_ix}{time_bin};...
                                                    repmat(sbj_ix,size(all_onsets{sbj_ix,time_bin,cond_ix}))'];
            elseif strcmp(plt_vars.grp_metric,'mdn')
                plot_onsets{cond_ix}{time_bin} = [plot_onsets{cond_ix}{time_bin}; nanmedian(all_onsets{sbj_ix,time_bin,cond_ix})];
            elseif strcmp(plt_vars.grp_metric,'avg')
                plot_onsets{cond_ix}{time_bin} = [plot_onsets{cond_ix}{time_bin}; nanmean(all_onsets{sbj_ix,time_bin,cond_ix})];
            else
                error(['Unknown plt_vars.grp_metric: ' plt_vars.grp_metric]);
            end
            % Report results in text
            %         fprintf('%s , %s: %f (N=%i)\n',SBJs{sbj_ix},groi_list{groi_ix},...
            %             median_onsets(sbj_ix,groi_ix),numel(cond_g_onsets{sbj_ix,groi_ix}));
            %         disp(cond_g_onsets{sbj_ix,groi_ix});
            %         fprintf('\n');
        end
    end
    plot_onsets{cond_ix} = padcat(plot_onsets{cond_ix}{:});
    [~,good_tbin_map{cond_ix}] = find(~all(isnan(plot_onsets{cond_ix}),1)==1);
    if strcmp(plt_vars.grp_metric,'all')
        plot_onset_sbj{cond_ix} = padcat(plot_onset_sbj{cond_ix}{:});
    end
end
save([root_dir 'PRJ_Stroop/data/tmp_onsets_sig_ch_' clust_id '.mat'],'-v7.3','sig_ch','plot_onsets','all_onsets');

%% Plot GROI Results
for cond_ix = 1:numel(cond_lab)
    % Create and format the plot
    fig_name = ['GRP' plt_vars.grp_metric '_HFA_onsets_' cond_lab{cond_ix} '_' clust_id '_' roi_id '_' event_lab...
        '_GM' num2str(gm_thresh) '_normRTout'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
    fprintf('Printing %i onsets across %i groups in %s\n',sum(sum(~isnan(plot_onsets{cond_ix}(:,good_tbin_map{cond_ix})))),...
                    numel(good_tbin_map{cond_ix}),fig_name);
                
    violins = violinplot(plot_onsets{cond_ix}(:,good_tbin_map{cond_ix}),clust_names(good_tbin_map{cond_ix}),...
                            'ViolinAlpha',0.3);
                        
    % Adjust plot propeties
    for time_bin = 1:numel(good_tbin_map{cond_ix})
        if strcmp(plt_vars.violin_scat_colors,'SBJ')
            % Change scatter colors to mark SBJ
            violins(time_bin).ViolinColor = [0.8 0.8 0.8];
            violins(time_bin).BoxPlot.FaceColor = bin_colors{good_tbin_map{cond_ix}(time_bin)};
            violins(time_bin).EdgeColor = bin_colors{good_tbin_map{cond_ix}(time_bin)};
            if sum(~isnan(plot_onsets{cond_ix}(:,good_tbin_map{cond_ix}(time_bin))))>1
                scat_colors = zeros([numel(violins(time_bin).ScatterPlot.XData) 3]);
                for sbj_ix = 1:numel(SBJs)
                    scat_colors(plot_onset_sbj{cond_ix}(:,good_tbin_map{cond_ix}(time_bin))==sbj_ix,:) = repmat(SBJ_colors(sbj_ix,:),...
                        sum(plot_onset_sbj{cond_ix}(:,good_tbin_map{cond_ix}(time_bin))==sbj_ix),1);
                end
                violins(time_bin).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
                violins(time_bin).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
                violins(time_bin).ScatterPlot.CData = scat_colors;
            else   % violin.MedianColor is plotted over only ScatterPlot point
                violins(time_bin).MedianColor = ...
                    SBJ_colors(plot_onset_sbj{cond_ix}(1,good_tbin_map{cond_ix}(time_bin)),:);
            end
        else
            % Change violin color to match ROI
            violins(time_bin).ViolinColor = bin_colors{good_tbin_map{cond_ix}(time_bin)};
        end
    end
    
    %% Save figure
    if save_fig
        fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_onsets_ROI/'...
            stat_id '/' roi_id '/' an_id '/' plt_id '/'];
        if ~exist(fig_dir,'dir')
            [~] = mkdir(fig_dir);
        end
        
        fig_filename = [fig_dir fig_name '.' fig_filetype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

end
