function SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,proc_id,an_id,roi_id,...
                                                    atlas_id)%,gm_thresh,plt_id,save_fig,fig_vis) %,fig_ftype)
% Load HFA analysis results for active and condition-differentiating
%   epochs, plot a summary of those time period per electrode
% clear all; %close all;
% fig_ftype = 'png';
% label_spacer = 0;
% groi_label_spacer = '      ';
% if ischar(save_fig); save_fig = str2num(save_fig); end
% if isnumeric(actv_win); actv_win = num2str(actv_win); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Prep variables
% an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
% eval(an_vars_cmd);
% plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
% eval(plt_vars_cmd);

% Get condition info
load([root_dir,'PRJ_Stroop/data/',SBJs{1},'/04_proc/',SBJs{1},'_smANOVA_ROI_',stat_id,'_',an_id,'.mat'],'st');
cond_lab = st.groups;

% Get event timing
mean_RTs = zeros(size(SBJs));

% Load all ROI info
[roi_list, roi_colors] = fn_roi_label_styles(roi_id);
if any(strcmp(roi_id,{'mgROI','gROI','main3','lat','deep','gPFC'}))
    roi_field = 'gROI';
else
    roi_field = 'ROI';
end

%% Load Results
% Set up onset counts
all_onsets          = cell([numel(SBJs) numel(roi_list) numel(st.groups)]);
all_onset_elec_lab  = cell([numel(SBJs) numel(roi_list) numel(st.groups)]);
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    % Load variables
    SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
    eval(SBJ_vars_cmd);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    srate = trial_info.sample_rate;
    
    % Load data
    load(strcat(SBJ_vars.dirs.proc,SBJ,'_smANOVA_ROI_',stat_id,'_',an_id,'.mat'));
    % Compute mean RT
    if strcmp(st.evnt_lab,'S')
        % Compute mean RT
        mean_RTs(sbj_ix) = mean(trial_info.response_time);
    end
    
    %% Load ROI and GM/WM info
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_pat_' atlas_id '_final.mat'];
    load(elec_fname);
    
    % Sort elecs by stat labels
    if numel(elec.label)~=numel(w2.label) || ~all(strcmp(elec.label,w2.label))
        error('elec w2 label mismatch');
        cfgs = []; cfgs.channel = stat.label;
        elec = fn_select_elec(cfgs,elec);
    end
    
%     % Get GM probability from tissue labels {'GM','WM','CSF','OUT'}
%     gm_bin  = elec.tissue_prob(:,1)>gm_thresh;
        
    %% Aggregate results per ROI
    for ch_ix = 1:numel(w2.label)
        % If elec matches roi_list and is in GM, get stats
        if any(strcmp(elec.(roi_field){ch_ix},roi_list)) %&& w2.max_hfa_z(ch_ix)>=z_thresh% && gm_bin(ch_ix)
            roi_ix = find(strcmp(elec.(roi_field){ch_ix},roi_list));
            % Get ANOVA group onsets
            for grp_ix = 1:numel(st.groups)
                if any(squeeze(w2.qval(grp_ix,ch_ix,:))<=st.alpha)
                    sig_onsets = w2.win_lim_s(squeeze(w2.qval(grp_ix,ch_ix,:))<=st.alpha,1);
                    if strcmp(st.evnt_lab,'R')
                        all_onsets{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onsets{sbj_ix,roi_ix,grp_ix} sig_onsets(1)];
                        all_onset_elec_lab{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onset_elec_lab{sbj_ix,roi_ix,grp_ix} w2.label(ch_ix)];
                    elseif ~isempty(strfind(st.ep_lab,'S')) && (sig_onsets(1)<mean_RTs(sbj_ix))
                        all_onsets{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onsets{sbj_ix,roi_ix,grp_ix} sig_onsets(1)];
                        all_onset_elec_lab{sbj_ix,roi_ix,grp_ix} = ...
                            [all_onset_elec_lab{sbj_ix,roi_ix,grp_ix} w2.label(ch_ix)];
                    else
                        error('st.ep_lab didnt match R or S*');
                    end
                end
            end
            
%             % Get RT correlation onset
%             if sum(squeeze(stat.mask(ch_ix,1,:)))>0
%                 mask_chunks = fn_find_chunks(squeeze(stat.mask(ch_ix,1,:)));
%                 mask_chunks(squeeze(stat.mask(ch_ix,1,mask_chunks(:,1)))==0,:) = [];
%                 % Convert the first onset of significance to time
%                 onset_time = stat.time(mask_chunks(1,1));
%                 % Exclude differences after the mean RT for this SBJ
%                 if strcmp(event_type,'resp') && (onset_time<0)
%                     all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
%                         [all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} onset_time];
%                     all_onset_elec_lab{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
%                         [all_onset_elec_lab{sbj_ix,roi_ix,numel(grp_lab)+1} stat.label(ch_ix)];
%                 elseif strcmp(event_type,'stim') && (onset_time<mean_RTs(sbj_ix))
%                     all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
%                         [all_onsets{sbj_ix,roi_ix,numel(grp_lab)+1} onset_time];
%                     all_onset_elec_lab{sbj_ix,roi_ix,numel(grp_lab)+1} = ...
%                         [all_onset_elec_lab{sbj_ix,roi_ix,numel(grp_lab)+1} stat.label(ch_ix)];
%                 end
%             end
        end
    end
    
    clear SBJ SBJ_vars hfa rt elec w2 st trial_info
end

%% Compute median onsets and onset differences per ROI pair
roi_pairs = nchoosek(1:numel(roi_list),2);
median_onsets = NaN([numel(SBJs) numel(roi_list) numel(cond_lab)]);
onset_diffs   = cell([numel(SBJs) size(roi_pairs,1) numel(cond_lab)]);
for sbj_ix = 1:numel(SBJs)
    for pair_ix = 1:size(roi_pairs,1)
        for grp_ix = 1:numel(cond_lab)
            median_onsets(sbj_ix,roi_ix,grp_ix) = nanmedian(all_onsets{sbj_ix,roi_ix,grp_ix});
            
            % If elecs in both ROIs
            if ~isempty(all_onsets{sbj_ix,roi_pairs(pair_ix,1),grp_ix}) && ...
               ~isempty(all_onsets{sbj_ix,roi_pairs(pair_ix,2),grp_ix})
                roi1_elecs = all_onset_elec_lab{sbj_ix,roi_pairs(pair_ix,1),grp_ix};
                roi2_elecs = all_onset_elec_lab{sbj_ix,roi_pairs(pair_ix,2),grp_ix};
                if numel(roi1_elecs)~=numel(unique(roi1_elecs))
                    error(['repeated onsets in roi1_elecs for ' SBJs{sbj_ix} ':' [roi1_elecs{:}]]);
                end
                if numel(roi2_elecs)~=numel(unique(roi2_elecs))
                    error(['repeated onsets in roi2_elecs for ' SBJs{sbj_ix} ': ' [roi2_elecs{:}]]);
                end
                % Compute difference for all pairs
                for roi1_eix = 1:numel(roi1_elecs)
                    for roi2_eix = 1:numel(roi2_elecs)
                        onset_diffs{sbj_ix,pair_ix,grp_ix} = [onset_diffs{sbj_ix,pair_ix,grp_ix}...
                                                            all_onsets{sbj_ix,roi_pairs(pair_ix,1),grp_ix}(roi1_eix) - ...
                                                            all_onsets{sbj_ix,roi_pairs(pair_ix,2),grp_ix}(roi2_eix)];
                    end
                end
            end
        end
%         fprintf('%s , %s: %f (N=%i)\n',SBJs{sbj_ix},groi_list{groi_ix},...
%             median_onsets(sbj_ix,groi_ix),numel(cond_g_onsets{sbj_ix,groi_ix}));
%         disp(cond_g_onsets{sbj_ix,groi_ix});
%         fprintf('\n');
    end
end

%% Print onset differences
hist_bins = linspace(-0.3,0.3,10);
for cond_ix = 1:numel(cond_lab)
    fprintf('========================== %s ==========================\n',cond_lab{cond_ix});
%     fig_name = ['GRP_HFA_onset_diffs_' cond_lab{cond_ix} '_' event_lab ...
%         '_GM' num2str(gm_thresh)];
%     figure('Name',fig_name,'units','normalized',...
%         'outerposition',[0 0 0.9 1],'Visible',fig_vis);
    for pair_ix = 1:size(roi_pairs,1)
%         ax = subplot(1,size(roi_pairs,1),pair_ix);
        sbj_cnt = 0;
        fprintf('-------------------------- %s-%s --------------------------\n',...
            roi_list{roi_pairs(pair_ix,1)},roi_list{roi_pairs(pair_ix,2)});
        for sbj_ix = 1:numel(SBJs)
            if ~isempty([onset_diffs{sbj_ix,pair_ix,cond_ix}])
                fprintf('%s (n_pairs = %i) mean = %.3f; sd = %.3f\n', ...
                    SBJs{sbj_ix}, numel(onset_diffs{sbj_ix,pair_ix,cond_ix}), ...
                    nanmean(onset_diffs{sbj_ix,pair_ix,cond_ix}), nanstd([onset_diffs{sbj_ix,pair_ix,cond_ix}]));
                sbj_cnt= sbj_cnt+1;
            end
        end
        % aggregate across subjects
        if ~isempty([onset_diffs{:,pair_ix,cond_ix}])
%             histogram([onset_diffs{:,pair_ix,cond_ix}],hist_bins);
%             line([nanmean([onset_diffs{:,pair_ix,cond_ix}]) ...
%                 nanmean([onset_diffs{:,pair_ix,cond_ix}])], ylim,...
%                 'Color','k','LineWidth',2);
            fprintf('GRP (n_pairs = %i, n_sbj = %i) mean = %.3f; sd = %.3f\n', ...
                numel([onset_diffs{:,pair_ix,cond_ix}]), sbj_cnt,...
                nanmean([onset_diffs{:,pair_ix,cond_ix}]), nanstd([onset_diffs{:,pair_ix,cond_ix}]));
        end
        ax.Title.String = [roi_list{roi_pairs(pair_ix,1)} '-' roi_list{roi_pairs(pair_ix,2)}];
    end
    fprintf('\n\n');
end

% %% Plot ROI Results
% error('copy over new version of plotting that uses elec and gm_thresh!');
% for cond_ix = 1:numel(cond_lab)
%     % Create and format the plot
%     fig_name = ['GRP_HFA_onsets_' cond_lab{cond_ix} '_' roi_id '_' event_lab '_meanRTout'];
%     figure('Name',fig_name,'units','normalized',...
%         'outerposition',[0 0 0.9 1],'Visible',fig_vis);
%     
%     % Plot Condition Difference Onsets per SBJ
%     ax = subplot(2,1,1);
%     hold on;
%     median_line_height = 0.3;
%     SBJ_y_spacer       = 0.35;
%     SBJ_ys             = SBJ_y_spacer*[1:numel(SBJs)]-median_line_height/3;
%     rt_marker_width    = 3;
%     rt_marker_size     = 250;
%     onset_marker_size  = 150;
%     marker_offset      = 0.08;
%     roi_y_offset       = linspace(-marker_offset,marker_offset,numel(roi_list));
%     s        = cell(size(roi_list));
%     lgd_lab  = cell(size(roi_list));
%     roi_flag = ones(size(roi_list));
%     for sbj_ix = 1:numel(SBJs)
%         for roi_ix = 1:numel(roi_list)
%             if ~isempty(all_onsets{sbj_ix,roi_ix,cond_ix})
%                 tmp_s = scatter(all_onsets{sbj_ix,roi_ix,cond_ix},...
%                     repmat(SBJ_ys(sbj_ix)+roi_y_offset(roi_ix),[1 numel(all_onsets{sbj_ix,roi_ix,cond_ix})]),...
%                     onset_marker_size,'o','filled');
%                 set(tmp_s,'MarkerFaceColor',roi_colors{roi_ix},'MarkerEdgeColor','k');
%                 
%                 s{roi_ix} = tmp_s(1);
%                 lgd_lab{roi_ix} = roi_list{roi_ix};
%             end
%         end
%         if strcmp(event_lab,'stim')
%             r = scatter(mean_RTs(sbj_ix)-plt_vars.plt_lim(1),SBJ_ys(sbj_ix),rt_marker_size,'+');
%             if sbj_ix==1
%                 s        = {s{:} r};
%                 roi_flag = [roi_flag 1];
%                 lgd_lab  = {lgd_lab{:} 'Mean RT'};
%             end
%         else
% %             error('plotting RT markers for resp-locked looks wrong...');
%             r = scatter(mean_RTs(sbj_ix),SBJ_ys(sbj_ix),rt_marker_size,'+');
%         end
%         set(r,'MarkerEdgeColor','k','LineWidth',rt_marker_width);
%         for roi_ix = 1:numel(roi_list)
%             line([median_onsets(sbj_ix,roi_ix,cond_ix) median_onsets(sbj_ix,roi_ix,cond_ix)],...
%                 [SBJ_ys(sbj_ix)-median_line_height/2 SBJ_ys(sbj_ix)+median_line_height/2],...
%                 'Color',roi_colors{roi_ix},'LineStyle','-','LineWidth',3);
%         end
%         line([plt_vars.plt_lim(1) plt_vars.plt_lim(2)],[SBJ_ys(sbj_ix) SBJ_ys(sbj_ix)],...
%             'Color',[0.2 0.2 0.2],'LineStyle',':');
%     end
%     for roi_ix = 1:numel(roi_list)
%         if isempty(lgd_lab{roi_ix})
%             lgd_lab{roi_ix} = [];
%             roi_flag(roi_ix) = 0;
%         end
%     end
%     legend([s{logical(roi_flag)}],lgd_lab{logical(roi_flag)},'Location','northeast');
%     
%     % ax = gca;
%     % Plot labels
%     % ax.XLabel.String   = 'Time (s)';
%     % ax.XLabel.FontSize = 14;
%     ax.XLim    = plt_vars.plt_lim;
%     ax.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
%     ax.XColor  = 'k';
%     
%     ax.YLabel.String   = 'Subject';
%     ax.YLabel.FontSize = 14;
%     ax.YLim            = [0 max(SBJ_ys)+median_line_height/2];
%     ax.YTick           = SBJ_ys;        % Ticks anywhere non-zero in plot_idx
%     ax.YTickLabel      = SBJs;
%     % ax.YTickLabelRotation = 45;
%     ax.YColor  = 'k';
%     
%     ax.Title.String = [cond_lab{cond_ix} ' Onsets by ' roi_id ' in Individuals'];
%     ax.Title.FontSize = 16;
%     
%     % Plot histogram of the median onsets by ROI across subjects
%     ax2 = subplot(2,1,2);
%     hold on;
%     % ax2 = gca;
%     b = [];
%     l = [];
%     % hist_alpha = 0.6;
%     hist_data = zeros([numel(roi_list) numel(plt_vars.x_data)-1]);
%     for roi_ix = 1:numel(roi_list)
%         hist_data(roi_ix,:)  = histcounts(median_onsets(:,roi_ix,cond_ix),plt_vars.x_data);
%     end
%     b = bar(plt_vars.x_data(1:end-1)+(diff(plt_vars.x_data(1:2))/2),hist_data',1,'stacked');
%     for roi_ix = numel(roi_list):-1:1 %backwards so OFC doesn't overwrite LPFC mean onset
%         set(b(roi_ix),'FaceColor',roi_colors{roi_ix},'EdgeColor','k');
%         l(roi_ix) = line([nanmean(median_onsets(:,roi_ix,cond_ix)) nanmean(median_onsets(:,roi_ix,cond_ix))],...
%             ax2.YLim,'Color',roi_colors{roi_ix},'LineStyle','-','LineWidth',3);
%         fprintf('%s mean(RT) for %s:\t%f\n',cond_lab{cond_ix},roi_list{roi_ix},...
%             nanmean(median_onsets(:,roi_ix,cond_ix)));
%     end
%     legend(b,roi_list{:},'Location','northeast');
%     % Plot labels
%     ax2.XLabel.String   = 'Time (s)';
%     ax2.XLabel.FontSize = 14;
%     ax2.XLim    = plt_vars.plt_lim;
%     ax2.XTick   = plt_vars.plt_lim(1):plt_vars.x_step_sz:plt_vars.plt_lim(2);
%     ax2.XColor  = 'k';
%     
%     ax2.YLabel.String   = '# Median Onsets';
%     ax2.YLabel.FontSize = 14;
%     % ax2.YLim            = [0 numel(SBJs)+1];
%     ax2.YTick           = 1:1:ax2.YLim(2);        % Ticks anywhere non-zero in plot_idx
%     % ax2.YTickLabel      = SBJs;
%     % ax2.YTickLabelRotation = 45;
%     ax2.YColor  = 'k';
%     
%     ax2.Title.String = ['Group-Level Median ' cond_lab{cond_ix} ' Onsets by ' roi_id];
%     ax2.Title.FontSize = 16;
%     
%     % Reposition axes
%     set(ax, 'Position',[0.05 0.33 0.9 0.63]);    %[left bottom width height]
%     set(ax2,'Position',[0.05 0.05 0.9 0.21]);    %[left bottom width height]
%     
%     %% Save figure
%     if save_fig
%         fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/GRP_onsets_ROI/'...
%             stat_id '_' roi_id '_' an_id '/'];
%         if ~exist(fig_dir,'dir')
%             mkdir(fig_dir);
%         end
%         
%         fig_filename = [fig_dir fig_name '.' fig_ftype];
%         fprintf('Saving %s\n',fig_filename);
%         saveas(gcf,fig_filename);
%         %eval(['export_fig ' fig_filename]);
%     end
% end

end
