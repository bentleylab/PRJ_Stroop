function SBJ10c_HFA_GRP_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,proc_id,an_id,median_yn,roi_id,...
                                                    atlas_id,gm_thresh,plt_id,save_fig,fig_vis,fig_ftype)
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot music staff style summary of those time period per electrode
%   also has a histogram across patients summary at the bottom
% Normalize all onset times by mean(RT)
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
addpath(ft_dir);
ft_defaults

%% Prep variables
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Get condition info
[grp_lab, grp_colors, grp_style] = fn_group_label_styles(st.model_lab);
cond_lab = [grp_lab, {'corr(RT)'}];

% Get event timing
mean_RTs = zeros(size(SBJs));
if strcmp(an_id(1:5),'HGm_S')
    evnt_lab = 'S';
elseif strcmp(an_id(1:5),'HGm_R')
    evnt_lab = 'R';
end

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(atlas_id,{'DK','Dx'}))
    view_space = 'pat';
elseif any(strcmp(atlas_id,{'Yeo7','Yeo17'}))
    view_space = 'mni_v';
else
    error(['Unknown atlas_id: ' atlas_id]);
end

% Set up onset counts
all_onsets  = cell([numel(SBJs) numel(roi_list) numel(grp_lab)+1]);

%% Load Results
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    srate = trial_info.sample_rate;
    
    % Compute mean RT
    if strcmp(evnt_lab,'S')
        % Compute mean RT
        mean_RTs(sbj_ix) = mean(trial_info.response_time);
    end
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_mANOVA_ROI_',stat_id,'_',an_id,'.mat'));
    
    % Get Sliding Window Parameters
    win_lim    = fn_sliding_window_lim(stat.time,round(st.win_len*srate),round(st.win_step*srate));
    win_center = round(mean(win_lim,2));
    
    %% Load ROI and GM/WM info
    if any(strcmp(atlas_id,{'DK','Dx'})) % these have tissue probabilities
        elec_tis_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' atlas_id '_full.mat'];
        load(elec_tis_fname);
    else
        % Load tissue prob
        elec_tis_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_Dx_tis.mat'];
        load(elec_tis_fname);
        tiss_prob = elec.tissue_prob;
        elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' atlas_id '.mat'];
        load(elec_fname);
        elec.tissue_prob = tiss_prob;
    end
    
    % Sort elecs by stat labels
    cfgs = []; cfgs.channel = stat.label;
    elec    = fn_select_elec(cfgs,elec);
    roi_lab = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,roi_id);
    
    % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
    gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
        
    %% Aggregate results per ROI
    for ch_ix = 1:numel(stat.label)
        % If elec matches roi_list and is in GM, get stats
        if any(strcmp(roi_lab{ch_ix},roi_list)) && gm_bin(ch_ix)
            roi_ix = find(strcmp(roi_lab{ch_ix},roi_list));
            % Get ANOVA group onsets
            for grp_ix = 1:numel(grp_lab)
                if any(squeeze(w2.qval(grp_ix,ch_ix,:))<st.alpha)
                    sig_onsets = stat.time(win_lim(squeeze(w2.qval(grp_ix,ch_ix,:))<st.alpha,1));
                    if strcmp(evnt_lab,'R')
                        all_onsets{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onsets{sbj_ix,roi_ix,grp_ix} sig_onsets(1)];
                    elseif strcmp(evnt_lab,'S') && (sig_onsets(1)<mean_RTs(sbj_ix))
                        all_onsets{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onsets{sbj_ix,roi_ix,grp_ix} sig_onsets(1)];
                    end
                end
            end
            
            % Get RT correlation onset
            if sum(squeeze(stat.mask(ch_ix,1,:)))>0
                mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
                mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
                % Convert the first onset of significance to time
                onset_time = stat.time(mask_chunks(1,1));
                % Exclude differences after the mean RT for this SBJ
                if strcmp(evnt_lab,'R')
                    all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
                        [all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} onset_time];
                elseif strcmp(evnt_lab,'S') && (onset_time<mean_RTs(sbj_ix))
                    all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
                        [all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} onset_time];
                end
            end
        end
    end
    
    % Normalize all onset times by mean reaction time
    if strcmp(evnt_lab,'S')
        for roi_ix = 1:size(all_onsets,2)
            for cond_ix = 1:size(all_onsets,3)
                all_onsets{sbj_ix,roi_ix,cond_ix} = all_onsets{sbj_ix,roi_ix,cond_ix}./mean_RTs(sbj_ix);
            end
        end
    end
    clear SBJ SBJ_vars hfa stat elec w2 st trial_info
end

%% Compute medians per gROI
if median_yn
    median_onsets = NaN([numel(SBJs) numel(roi_list) numel(grp_lab)+1]);
    for sbj_ix = 1:numel(SBJs)
        for roi_ix = 1:numel(roi_list)
            for grp_ix = 1:numel(grp_lab)
                median_onsets(sbj_ix,roi_ix,grp_ix) = nanmedian(all_onsets{sbj_ix,roi_ix,grp_ix});
            end
            median_onsets(sbj_ix,roi_ix,numel(grp_lab)+1) = nanmedian(all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1});
            %         fprintf('%s , %s: %f (N=%i)\n',SBJs{sbj_ix},groi_list{groi_ix},...
            %             median_onsets(sbj_ix,groi_ix),numel(cond_g_onsets{sbj_ix,groi_ix}));
            %         disp(cond_g_onsets{sbj_ix,groi_ix});
            %         fprintf('\n');
        end
    end
end

%% Plot GROI Results
for cond_ix = 1:numel(cond_lab)
    % Create and format the plot
    fig_name = ['GRP_HFA_onsets_' cond_lab{cond_ix} '_' roi_id '_' evnt_lab...
        '_GM' num2str(gm_thresh) '_normRTout_med' num2str(median_yn)];
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 0.9 1],'Visible',fig_vis);
    
    % Plot Condition Difference Onsets per SBJ
    ax = subplot(2,1,1);
    hold on;
    median_line_height = 0.3;
    SBJ_y_spacer       = 0.35;
    SBJ_ys             = SBJ_y_spacer*[1:numel(SBJs)]-median_line_height/3;
    rt_marker_width    = 3;
    rt_marker_size     = 250;
    onset_marker_size  = 150;
    marker_offset      = 0.08;
    roi_y_offset       = linspace(-marker_offset,marker_offset,numel(roi_list));
    scat     = cell(size(roi_list));
    lgd_lab  = cell(size(roi_list));
    roi_flag = ones(size(roi_list));
    for sbj_ix = 1:numel(SBJs)
        for roi_ix = 1:numel(roi_list)
            if ~isempty(all_onsets{sbj_ix,roi_ix,cond_ix})
                tmp_s = scatter(all_onsets{sbj_ix,roi_ix,cond_ix},...
                    repmat(SBJ_ys(sbj_ix)+roi_y_offset(roi_ix),...
                    [1 numel(all_onsets{sbj_ix,roi_ix,cond_ix})]),...
                    onset_marker_size,'o','filled');
                set(tmp_s,'MarkerFaceColor',roi_colors{roi_ix},'MarkerEdgeColor','k');
                
                scat{roi_ix} = tmp_s(1);
                lgd_lab{roi_ix} = roi_list{roi_ix};
            end
        end
        if strcmp(evnt_lab,'S')
            r = scatter(1-plt_vars.plt_lim(1),SBJ_ys(sbj_ix),rt_marker_size,'+');
            if sbj_ix==1
                scat     = {scat{:} r};
                roi_flag = [roi_flag 1];
                lgd_lab  = {lgd_lab{:} 'Mean RT'};
            end
        else
%             error('plotting RT markers for resp-locked looks wrong...');
            r = scatter(mean_RTs(sbj_ix),SBJ_ys(sbj_ix),rt_marker_size,'+');
        end
        set(r,'MarkerEdgeColor','k','LineWidth',rt_marker_width);
        if median_yn
            for roi_ix = 1:numel(roi_list)
                line([median_onsets(sbj_ix,roi_ix,cond_ix) median_onsets(sbj_ix,roi_ix,cond_ix)],...
                    [SBJ_ys(sbj_ix)-median_line_height/2 SBJ_ys(sbj_ix)+median_line_height/2],...
                    'Color',roi_colors{roi_ix},'LineStyle','-','LineWidth',3);
            end
        end
        line([plt_vars.plt_lim(1) plt_vars.plt_lim(2)],[SBJ_ys(sbj_ix) SBJ_ys(sbj_ix)],...
            'Color',[0.2 0.2 0.2],'LineStyle',':');
    end
    for roi_ix = 1:numel(roi_list)
        if isempty(lgd_lab{roi_ix})
            lgd_lab{roi_ix} = [];
            roi_flag(roi_ix) = 0;
        end
    end
    legend([scat{logical(roi_flag)}],lgd_lab{logical(roi_flag)},'Location','northeast');
    
    % ax = gca;
    % Plot labels
%     ax.XLabel.String   = '% Reaction Time';
%     ax.XLabel.FontSize = 14;
    ax.XLim    = plt_vars.plt_lim;
    ax.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    ax.XColor  = 'k';
    
    ax.YLabel.String   = 'Subject';
    ax.YLabel.FontSize = 14;
    ax.YLim            = [0 max(SBJ_ys)+median_line_height/2];
    ax.YTick           = SBJ_ys;        % Ticks anywhere non-zero in plot_idx
    ax.YTickLabel      = SBJs;
    % ax.YTickLabelRotation = 45;
    ax.YColor  = 'k';
    
    ax.Title.String = [cond_lab{cond_ix} ' Onsets by ' roi_id ' in Individuals'];
    ax.Title.FontSize = 16;
    
    % Plot histogram of the median onsets by ROI across subjects
    ax2 = subplot(2,1,2);
    hold on;
    % ax2 = gca;
    b = [];
    l = [];
    % hist_alpha = 0.6;
    hist_data = zeros([numel(roi_list) numel(plt_vars.x_data)-1]);
    for roi_ix = 1:numel(roi_list)
        if median_yn
            hist_data(roi_ix,:)  = histcounts(median_onsets(:,roi_ix,cond_ix),plt_vars.x_data);
        else
            hist_data(roi_ix,:)  = histcounts([all_onsets{:,roi_ix,cond_ix}],plt_vars.x_data);
        end
    end
    b = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),hist_data',1,'stacked');
    for roi_ix = numel(roi_list):-1:1 %backwards so OFC doesn't overwrite LPFC mean onset
        set(b(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
        if median_yn
            l(roi_ix) = line([nanmean(median_onsets(:,roi_ix,cond_ix))...
                nanmean(median_onsets(:,roi_ix,cond_ix))],...
                ax2.YLim,'Color',roi_colors{roi_ix},'LineStyle','-','LineWidth',3);
            fprintf('%s mean(RT) for %s:\t%f\n',cond_lab{cond_ix},roi_list{roi_ix},...
                nanmean(median_onsets(:,roi_ix,cond_ix)));
        end
    end
    legend(b,roi_list{:},'Location','northeast');
    % Plot labels
    ax2.XLabel.String   = '% Reaction Time';
    ax2.XLabel.FontSize = 14;
    ax2.XLim    = plt_vars.plt_lim;
    ax2.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
    ax2.XColor  = 'k';
    
    if median_yn
        ax2.YLabel.String   = '# Median Onsets';
    else
        ax2.YLabel.String   = '# Onsets';
    end
    ax2.YLabel.FontSize = 14;
    % ax2.YLim            = [0 numel(SBJs)+1];
    ax2.YTick           = 1:1:ax2.YLim(2);        % Ticks anywhere non-zero in plot_idx
    % ax2.YTickLabel      = SBJs;
    % ax2.YTickLabelRotation = 45;
    ax2.YColor  = 'k';
    
    if median_yn
        ax2.Title.String = ['Group-Level Median ' cond_lab{cond_ix} ' Onsets by ' roi_id];
    else
        ax2.Title.String = ['Group-Level ' cond_lab{cond_ix} ' Onsets by ' roi_id];
    end
    ax2.Title.FontSize = 16;
    
    % Reposition axes
    set(ax, 'Position',[0.05 0.33 0.9 0.63]);    %[left bottom width height]
    set(ax2,'Position',[0.05 0.05 0.9 0.21]);    %[left bottom width height]
    
    %% Save figure
    if save_fig
        fig_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_onsets_ROI/'...
            stat_id '/' roi_id '/' an_id '/'];
        if ~exist(fig_dir,'dir')
            [~] = mkdir(fig_dir);
        end
        
        fig_fname = [fig_dir fig_name '.' fig_ftype];
        fprintf('Saving %s\n',fig_fname);
        saveas(gcf,fig_fname);
        %eval(['export_fig ' fig_fname]);
    end
end

end
