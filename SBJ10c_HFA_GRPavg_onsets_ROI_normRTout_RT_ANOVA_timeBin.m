function SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA_timeBin(SBJs,tbin_id,stat_id,proc_id,an_id,roi_id,...
                                                    atlas_id,gm_thresh,plt_id,save_fig,fig_vis,fig_ftype)
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot a summary of those time period per electrode
% Normalize all onset times by mean(RT)
% INPUTS:
%   plt_vars.grp_metric [str] - {'avg','mdn','all'}
%       mean/median will compute that metric within each SBJ (variance is across SBJs)
%       all- all electrode onsets are aggregated as if from the same SBJ

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
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

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
if strcmp(event_type,'stim')
    error('tbin for stim locked needs more thinking!');
end
if strcmp(tbin_id,'eqROI')
    n_tbins = numel(roi_list);
end

% Cluster info
clust_names  = cell([1 numel(roi_list)]);
for clust_ix = 1:numel(roi_list)
    clust_names{clust_ix} = ['C' num2str(clust_ix)];
end
if n_tbins==4
    bin_colors = [0 0 1; 0 1 1; 1 0 1; 1 0 0];  %garish blue, cyan, magenta, red
else
    bin_colors = parula(n_tbins);
end
% bin_colors = {[27 158 119]./255, [117 112 179]./255, [217 95 2]./255, [231 41 138]./255}; %qualitative max diff from gROI, cool to hot

%% Load Results
% Set up onset counts
all_onsets  = cell([numel(SBJs) n_tbins numel(cond_lab)]);
sig_ch      = cell([numel(SBJs) numel(cond_lab)]);
sig_ch_tbin = cell([numel(SBJs) numel(cond_lab)]);
sig_ch_onset = cell([numel(SBJs) numel(cond_lab)]);
if strcmp(plt_vars.violin_scat_colors,'ROI')
    all_onset_rois  = cell([numel(SBJs) numel(roi_list) numel(cond_lab)]);
end
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    
    % Compute mean RT
    if strcmp(event_type,'stim')
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
    win_lim    = fn_sliding_window_lim(stat.time,round(st.win_len*sample_rate),round(st.win_step*sample_rate));
    win_center = round(mean(win_lim,2));
    if strcmp(tbin_id,'eqROI')  %!!! check for 1st 2 letters = 'eq'
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
        end
    else
        error(['Unknown tbin_id: ' tbin_id]);
    end
    
    %% Load ROI and GM/WM info
    if any(strcmp(atlas_id,{'DK','Dx'})) % these have tissue probabilities
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' atlas_id '_full.mat'];
        load(elec_fname);
    else
        error('use Dx_full!');
%         % Load tissue prob
%         elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_Dx_tis.mat'];
%         load(elec_fname);
%         tiss_prob = elec.tissue_prob;
%         elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' atlas_id '.mat'];
%         load(elec_fname);
%         elec.tissue_prob = tiss_prob;
    end
    
    % Sort elecs by stat labels
    cfgs = []; cfgs.channel = stat.label;
    elec = fn_select_elec(cfgs,elec);
    elec.roi = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
    
%     % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
%     gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
        
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        % If elec matches roi_list and is in GM, get stats
        if any(strcmp(elec.roi{ch_ix},roi_list))% && gm_bin(ch_ix)
            roi_ix = find(strcmp(elec.roi{ch_ix},roi_list));
            % Get ANOVA group onsets
            for cond_ix = 1:numel(cond_lab)               
                % ANOVA conditions
                if any(strcmp(cond_lab{cond_ix},grp_lab)) && any(squeeze(qvals(cond_ix,ch_ix,:))<0.05)
                    sig_ch{sbj_ix,cond_ix} = [sig_ch{sbj_ix,cond_ix} stat.label(ch_ix)];
                    sig_onsets = stat.time(win_lim(squeeze(qvals(cond_ix,ch_ix,:))<0.05,1));
                    % Assign time bin by first onset
                    time_bin = find(find(squeeze(qvals(cond_ix,ch_ix,:))<0.05,1)<=peak_bins(cond_ix,:),1);
                    sig_ch_tbin{sbj_ix,cond_ix} = [sig_ch_tbin{sbj_ix,cond_ix} time_bin];
                    sig_ch_onset{sbj_ix,cond_ix} = [sig_ch_onset{sbj_ix,cond_ix} sig_onsets(1)];
                    if strcmp(event_type,'resp')
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} sig_onsets(1)];
                    elseif strcmp(event_type,'stim') && (sig_onsets(1)<mean_RTs(sbj_ix))
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} sig_onsets(1)];
                    end
                    % Grab ROI to plot individual scatter colors
                    if strcmp(plt_vars.violin_scat_colors,'ROI')
                        all_onset_rois{sbj_ix,time_bin,cond_ix}  = ...
                            [all_onset_rois{sbj_ix,time_bin,cond_ix} roi_ix];
                    end
                    
                % Get RT correlation onset
                elseif strcmp(cond_lab{cond_ix},'corr(RT)') && sum(squeeze(stat.mask(ch_ix,1,:)))>0
                    sig_ch{sbj_ix,cond_ix} = [sig_ch{sbj_ix,cond_ix} stat.label(ch_ix)];
                    mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
                    mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
                    % Convert the first onset of significance to time
                    onset_time = stat.time(mask_chunks(1,1));
                    % Assign time bin by first onset
                    time_bin = find(mask_chunks(1,1)<=peak_bins(cond_ix,:),1);
                    sig_ch_tbin{sbj_ix,cond_ix} = [sig_ch_tbin{sbj_ix,cond_ix} time_bin];
                    sig_ch_onset{sbj_ix,cond_ix} = [sig_ch_onset{sbj_ix,cond_ix} onset_time];
                    % Exclude differences after the mean RT for this SBJ
                    if strcmp(event_type,'resp')
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} onset_time];
                    elseif strcmp(event_type,'stim') && (onset_time<mean_RTs(sbj_ix))
                        all_onsets{sbj_ix,time_bin,cond_ix} = ...
                            [all_onsets{sbj_ix,time_bin,cond_ix} onset_time];
                    end
                    % Grab ROI to plot individual scatter colors
                    if strcmp(plt_vars.violin_scat_colors,'ROI')
                        all_onset_rois{sbj_ix,time_bin,cond_ix}  = ...
                            [all_onset_rois{sbj_ix,time_bin,cond_ix} roi_ix];
                    end
                end
            end
        end
    end
    
    % Normalize all onset times by mean reaction time
    if strcmp(event_type,'stim')
        for time_bin = 1:n_tbins
            for cond_ix = 1:numel(cond_lab)
                all_onsets{sbj_ix,time_bin,cond_ix} = all_onsets{sbj_ix,time_bin,cond_ix}./mean_RTs(sbj_ix);
            end
        end
    end
    clear SBJ SBJ_vars hfa stat elec w2
end

%% Aggregate/Process onsets per gROI
% Format as struct to fit violinplot
plot_onsets     = cell([numel(cond_lab) 1]);
plot_onset_scat = cell([numel(cond_lab) 1]);
good_tbin_map   = cell([numel(cond_lab) 1]);
for cond_ix = 1:numel(cond_lab)
    plot_onsets{cond_ix}  = {};
    plot_onset_scat{cond_ix}  = {};
    good_tbin_map{cond_ix} = zeros([1 n_tbins]);
    for time_bin = 1:n_tbins
        if strcmp(plt_vars.grp_metric,'all')
            plot_onsets{cond_ix}{time_bin} = [all_onsets{:,time_bin,cond_ix}]';
            plot_onset_scat{cond_ix}{time_bin} = [];
        else
            plot_onsets{cond_ix}{time_bin} = [];
        end
        for sbj_ix = 1:numel(SBJs)
            % Aggregate onsets per ROI within each SBJ
            if strcmp(plt_vars.grp_metric,'all')
                if strcmp(plt_vars.violin_scat_colors,'SBJ')
                    plot_onset_scat{cond_ix}{time_bin} = [plot_onset_scat{cond_ix}{time_bin};...
                                                    repmat(sbj_ix,size(all_onsets{sbj_ix,time_bin,cond_ix}))'];
                elseif strcmp(plt_vars.violin_scat_colors,'ROI')
                    plot_onset_scat{cond_ix}{time_bin} = [plot_onset_scat{cond_ix}{time_bin};...
                                                    all_onset_rois{sbj_ix,time_bin,cond_ix}'];
                end
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
        plot_onset_scat{cond_ix} = padcat(plot_onset_scat{cond_ix}{:});
    end
end
save_dir = [root_dir 'PRJ_Stroop/data/HFA_onsets/' stat_id '/' an_id '/'];
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
save([save_dir 'onsets_' tbin_id '.mat'],'-v7.3','SBJs','sig_ch','sig_ch_tbin','sig_ch_onset','plot_onsets','all_onsets');

%% Plot GROI Results
for cond_ix = 1:numel(cond_lab)
    % Create and format the plot
    fig_name = ['GRP' plt_vars.grp_metric '_HFA_onsets_' cond_lab{cond_ix} '_' tbin_id '_' roi_id '_' event_type...
        '_GM' num2str(gm_thresh) '_normRTout'];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.5 0.6],'Visible',fig_vis);
    fprintf('Printing %i onsets across %i groups in %s\n',sum(sum(~isnan(plot_onsets{cond_ix}(:,good_tbin_map{cond_ix})))),...
                    numel(good_tbin_map{cond_ix}),fig_name);
                
    violins = violinplot(plot_onsets{cond_ix}(:,good_tbin_map{cond_ix}),clust_names(good_tbin_map{cond_ix}),...
                            'ViolinAlpha',0.3);
                        
    % Adjust plot propeties
    for bin_ix = 1:numel(good_tbin_map{cond_ix})
        if strcmp(plt_vars.violin_scat_colors,'')
            % No individual scatter colors, so change violin color to match time bin
            violins(bin_ix).ViolinColor = bin_colors(good_tbin_map{cond_ix}(bin_ix),:);
        else
            % Adjust properties to see individual scatter colors
            violins(bin_ix).ViolinColor = [0.8 0.8 0.8];  % shaded region
            violins(bin_ix).BoxPlot.FaceColor = bin_colors(good_tbin_map{cond_ix}(bin_ix),:);
            violins(bin_ix).EdgeColor = bin_colors(good_tbin_map{cond_ix}(bin_ix),:);
            % Plot individual scatter colors
            if sum(~isnan(plot_onsets{cond_ix}(:,good_tbin_map{cond_ix}(bin_ix))))>1
                scat_colors = zeros([numel(violins(bin_ix).ScatterPlot.XData) 3]);
                if strcmp(plt_vars.violin_scat_colors,'SBJ')
                    for sbj_ix = 1:numel(SBJs)
                        scat_colors(plot_onset_scat{cond_ix}(:,good_tbin_map{cond_ix}(bin_ix))==sbj_ix,:) = repmat(SBJ_colors(sbj_ix,:),...
                            sum(plot_onset_scat{cond_ix}(:,good_tbin_map{cond_ix}(bin_ix))==sbj_ix),1);
                    end
                elseif strcmp(plt_vars.violin_scat_colors,'ROI')
                    for roi_ix = 1:numel(roi_list)
                        scat_colors(plot_onset_scat{cond_ix}(:,good_tbin_map{cond_ix}(bin_ix))==roi_ix,:) = repmat(roi_colors{roi_ix},...
                            sum(plot_onset_scat{cond_ix}(:,good_tbin_map{cond_ix}(bin_ix))==roi_ix),1);
                    end
                end
                violins(bin_ix).ScatterPlot.MarkerFaceColor = 'flat';   % Necessary for CData to work
                violins(bin_ix).ScatterPlot.MarkerEdgeColor = 'flat';   % Necessary for CData to work
                violins(bin_ix).ScatterPlot.CData = scat_colors;
            else
                % only 1 point, so violin.MedianColor is plotted over only ScatterPlot point
                if strcmp(plt_vars.violin_scat_colors,'SBJ')
                    violins(bin_ix).MedianColor = ...
                        SBJ_colors(plot_onset_scat{cond_ix}(1,good_tbin_map{cond_ix}(bin_ix)),:);
                elseif strcmp(plt_vars.violin_scat_colors,'ROI')
                    violins(bin_ix).MedianColor = ...
                        roi_colors{plot_onset_scat{cond_ix}(1,good_tbin_map{cond_ix}(bin_ix))};                    
                end
            end
        end
    end
    
    %% Save figure
    if save_fig
        fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_onsets_ROI/'...
            stat_id '/' roi_id '/' an_id '/' plt_id '/'];
        if ~exist(fig_dir,'dir')
            [~] = mkdir(fig_dir);
        end
        
        fig_filename = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_filename);
        saveas(gcf,fig_filename);
        %eval(['export_fig ' fig_filename]);
    end
end

end
