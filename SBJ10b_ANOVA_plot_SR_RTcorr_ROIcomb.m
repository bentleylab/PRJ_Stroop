function SBJ10b_ANOVA_plot_SR_RTcorr_ROIcomb(SBJ,pipeline_id,stat_id,an_id_s,an_id_r,...
                                        atlas_id,roi_id,plt_id,save_fig,fig_vis,fig_filetype)
% Plots ANOVA results
% clear all; %close all;
% fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);plt_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
eval(plt_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');
cond_lab = [grp_lab rt_lab];
event_lab = {'stim', 'resp'};

[roi_list, roi_colors, ~] = fn_roi_label_styles(roi_id);

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

f_name_s = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id_s '.mat'];
f_name_r = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id_r '.mat'];
tmp = load(f_name_s,'w2'); w2{1} = tmp.w2;
tmp = load(f_name_r,'w2'); w2{2} = tmp.w2;
% tmp = load(f_name_s,'hfa'); hfa{1} = tmp.hfa;
% tmp = load(f_name_r,'hfa'); hfa{2} = tmp.hfa;
tmp = load(f_name_s,'stat'); stat{1} = tmp.stat;
tmp = load(f_name_r,'stat'); stat{2} = tmp.stat;
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(stat{1}.time)-1)/(stat{1}.time(end)-stat{1}.time(1));

% Load ROI and GM/WM info
elec_tis_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_pat_' atlas_id '_tis.mat'];
load(elec_tis_fname);

% Sort elecs by stat labels
cfgs = []; cfgs.channel = stat{1}.label;
elec    = fn_select_elec(cfgs,elec);
roi_lab = fn_atlas2roi_labels(elec.atlas_label,atlas_id,roi_id);

%% Prep Data
% FDR correct pvalues for ANOVA
win_lim = {}; win_center = {};
qvals = {NaN(size(w2{1}.pval)) NaN(size(w2{2}.pval))};
for sr_ix = 1:2
    for ch_ix = 1:numel(stat{1}.label)
        pvals = squeeze(w2{sr_ix}.pval(:,ch_ix,:));
        [~, ~, ~, qvals{sr_ix}(:,ch_ix,:)] = fdr_bh(pvals);%,0.05,'pdep','yes');
    end
    
    % Get Sliding Window Parameters
    win_lim{sr_ix}    = fn_sliding_window_lim(stat{sr_ix}.time,round(win_len*sample_rate),round(win_step*sample_rate));
    win_center{sr_ix} = round(mean(win_lim{sr_ix},2));
    
    % Convert % explained variance to 0-100 scale
    w2{sr_ix}.trial = w2{sr_ix}.trial*100;
end

% Trim data to plotting epoch
%   NOTE: stat should be on stat_lim(1):stat_lim(2)+0.001 time axis
%   w2 should fit within that since it's averaging into a smaller window
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
% hfa{1}  = ft_selectdata(cfg_trim,hfa{1});
stat{1} = ft_selectdata(cfg_trim,stat{1});
cfg_trim.latency = plt_vars.plt_lim_R;
% hfa{2}  = ft_selectdata(cfg_trim,hfa{2});
stat{2} = ft_selectdata(cfg_trim,stat{2});

% % Compute mean RT per condition
% RTs = round(1000*trial_info.response_time); % converts sec to ms
% for cond_ix = 1:numel(grp_lab{congr_ix})
%     RT_means{cond_ix} = mean(RTs(fn_condition_index([grp_lab{congr_ix}{cond_ix}], trial_info.condition_n)==1));
%     % Add in the baseline offset to plot correctly
%     RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*1000;
% end

%% Plot Results
fig_dir = [root_dir 'PRJ_Stroop/results/HFA/' SBJ '/' stat_id '/SR/' an_id_s '-' an_id_r '/ROIcomb/'];
if ~exist(fig_dir,'dir')
    [~] = mkdir(fig_dir);
end

% Find plot limits
max_w2 = max([max(max(max(w2{1}.trial))) max(max(max(w2{2}.trial)))]);
min_w2 = min([min(min(min(w2{1}.trial))) min(min(min(w2{2}.trial)))]);
ylim1_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims1  = [min_w2-ylim1_fudge max_w2+ylim1_fudge];
yticks1 = 0:1:ylims1(2);

max_rho = max([max(max(squeeze(stat{1}.rho))) max(max(squeeze(stat{2}.rho)))]);
min_rho = min([min(min(squeeze(stat{1}.rho))) min(min(squeeze(stat{2}.rho)))]);
ylims2  = [round(min_rho*10)/10-0.1 round(max_rho*10)/10+0.1]; % extra on top and bottom for StdErr
yticks2 = ylims2(1):0.1:ylims2(2);
y_sig = zeros([1 numel(grp_lab)+1]);
y_sig(1) = mean([min_w2,max_w2]);
for grp_ix = 2:numel(grp_lab)+2
    y_sig(grp_ix) = y_sig(grp_ix-1)+ylim1_fudge;
end

% Create a figure for each statistic and each ROI
for cond_ix = 1:numel(cond_lab)
    %     sig_ch = cell([1+numel(grp_lab) 2]);
    fig_name = [SBJ '_ANOVA_' stat_id '_SR_' cond_lab{cond_ix} '_' roi_id];
    f = figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.5],'Visible',fig_vis); %twice as wide for the double plot
    any_plot = 0;
    for sr_ix = 1:2
        ax = subplot(1,2,sr_ix);
        hold on;
        
        for ch_ix = 1:numel(stat{sr_ix}.label)
            roi_ix = find(strcmp(roi_list,roi_lab{ch_ix}));
            %                 % Plot var_exp
            %                 [axs,grp_lines,rt_line] = plotyy(...
            %                     repmat(win_center{sr_ix},size(grp_lab)),squeeze(w2{sr_ix}.trial(:,ch_ix,:))',...
            %                     1:numel(stat{sr_ix}.time),squeeze(stat{sr_ix}.rho(ch_ix,:,:))');
            %                 for grp_ix = 1:numel(grp_lab)
            %                     set(grp_lines(grp_ix),'Color',[grp_colors{grp_ix}],'LineStyle',grp_style{grp_ix});
            %                 end
            %                 set(rt_line,'Color',rt_color{1},'LineStyle',rt_style{1});
            %                 main_lines = [main_lines grp_lines' rt_line];
            %                 lgd_lab = {grp_lab{:} 'corr(RT)'};
            
            % Plot significant time periods
            if cond_ix <= numel(grp_lab)
                ylims = ylims1;
                yticks = yticks1;
                ylab = '% Variance Explained';
                if any(squeeze(qvals{sr_ix}(cond_ix,ch_ix,:))<0.05) && ~isempty(roi_ix)
                    any_plot = 1;
                    plot(win_center{sr_ix},squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,:))',...
                        'Color',[roi_colors{roi_ix}],'LineStyle',grp_style{cond_ix});
                    %                         % Track channels with significant var_exp
                    %                         if ~iscell(sig_ch{cond_ix,sr_ix})
                    %                             sig_ch{cond_ix,sr_ix} = {};
                    %                         end
                    %                         sig_ch{cond_ix,sr_ix} = {sig_ch{cond_ix,sr_ix}{:} stat{sr_ix}.label{ch_ix}};
                    
                    % Find significant periods
                    if strcmp(plt_vars.sig_type,'bold')
                        sig_chunks = fn_find_chunks(squeeze(qvals{sr_ix}(cond_ix,ch_ix,:))<0.05);
                        sig_chunks(squeeze(qvals{sr_ix}(cond_ix,ch_ix,sig_chunks(:,1)))>0.05,:) = [];
                        for sig_ix = 1:size(sig_chunks,1)
                            line(win_center{sr_ix}(sig_chunks(sig_ix,1):sig_chunks(sig_ix,2)),...
                                squeeze(w2{sr_ix}.trial(cond_ix,ch_ix,sig_chunks(sig_ix,1):sig_chunks(sig_ix,2))),...
                                'Color',[roi_colors{roi_ix}],'LineStyle',plt_vars.sig_style,...
                                'LineWidth',plt_vars.sig_width);
                        end
                    elseif strcmp(plt_vars.sig_type,'scatter')
                        sig_times = win_center{sr_ix}(squeeze(qvals{sr_ix}(cond_ix,ch_ix,:))<0.05);
                        scatter(sig_times,repmat(y_sig(cond_ix),size(sig_times)),...
                            plt_vars.sig_scat_size,roi_colors{roi_ix},plt_vars.sig_scat_mrkr);
                    elseif strcmp(plt_vars.sig_type,'patch')
                        sig_chunks = fn_find_chunks(squeeze(qvals{sr_ix}(cond_ix,ch_ix,:))<0.05);
                        sig_chunks(squeeze(qvals{sr_ix}(cond_ix,ch_ix,sig_chunks(:,1)))>0.05,:) = [];
                        error('sig_type = patch needs sig_times variable!');
                        % Plot Significance Shading
                        fprintf('%s %s (%s) -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
                            stat{sr_ix}.label{ch_ix},grp_lab{cond_ix},event_lab{sr_ix},size(sig_chunks,1));
                        for sig_ix = 1:size(sig_chunks,1)
                            sig_times = win_lim{sr_ix}(sig_chunks(sig_ix,:),:);
                            patch([sig_times(1,1) sig_times(1,1) sig_times(2,2) sig_times(2,2)], ...
                                [ylims1(1) ylims1(2) ylims1(2) ylims1(1)],...
                                roi_colors{roi_ix},'FaceAlpha',plt_vars.sig_alpha);
                        end
                    end
                    %                     else
                    %                         fprintf('%s %s (%s) -- NO SIGNIFICANT CLUSTERS FOUND...\n',...
                    %                             stat{sr_ix}.label{ch_ix},grp_lab{cond_ix},event_lab{sr_ix});
                end
            else
                % RT correlation significant time periods
                ylims = ylims2;
                yticks = yticks2;
                ylab = 'Correlation with RT';
                if sum(stat{sr_ix}.mask(ch_ix,:))>0 && ~isempty(roi_ix)
                    any_plot = 1;
                    %                         % Track channels with significant var_exp
                    %                         if ~iscell(sig_ch{cond_ix,sr_ix})
                    %                             sig_ch{cond_ix,sr_ix} = {};
                    %                         end
                    %                         sig_ch{cond_ix,sr_ix} = {sig_ch{cond_ix,sr_ix}{:} stat{sr_ix}.label{ch_ix}};
                    plot(1:numel(stat{sr_ix}.time),squeeze(stat{sr_ix}.rho(ch_ix,:,:))',...
                        'Color',[roi_colors{roi_ix}],'LineStyle','-');
                    
                    sig_chunks = fn_find_chunks(squeeze(stat{sr_ix}.mask(ch_ix,:,:)));
                    sig_chunks(squeeze(stat{sr_ix}.mask(ch_ix,:,sig_chunks(:,1)))==0,:) = [];
                    sig_times = find(squeeze(stat{sr_ix}.mask(ch_ix,:,:)==1));
                    if strcmp(plt_vars.sig_type,'bold')
                        for sig_ix = 1:size(sig_chunks,1)
                            sig_times = sig_chunks(sig_ix,1):sig_chunks(sig_ix,2);
                            line(sig_times,squeeze(stat{sr_ix}.rho(ch_ix,:,sig_times)),...
                                'Color',roi_colors{roi_ix},'LineStyle',plt_vars.sig_style,...
                                'LineWidth',plt_vars.sig_width);
                        end
                    elseif strcmp(plt_vars.sig_type,'scatter')
                        scatter(sig_times,repmat(y_sig(cond_ix),size(sig_times)),...
                            plt_vars.sig_scat_size2,roi_colors{roi_ix},plt_vars.sig_scat_mrkr2);
                    end
                    fprintf('%s RT (%s) -- SIGNIFICANT CLUSTERS FOUND, plotting with significance shading...\n',...
                        stat{sr_ix}.label{ch_ix},event_lab{sr_ix});%,size(sig_chunks,1));
                    %                     else
                    %                         fprintf('%s RT (%s) -- NO SIGNIFICANT CLUSTERS FOUND, plotting without significance shading...\n',...
                    %                             stat{sr_ix}.label{ch_ix},event_lab{sr_ix});
                end
            end
        end
        
        % Plot event
        if strcmp(event_lab{sr_ix},'stim')
            x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
        else
            x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Plot Response Marker
            event_line = line([find(stat{sr_ix}.time==0) find(stat{sr_ix}.time==0)],ylims,...
                'LineWidth',plt_vars.evnt_width, 'Color',plt_vars.evnt_color,...
                'LineStyle',plt_vars.evnt_style);
            %                 main_lines = [main_lines event_line];
            %                 lgd_lab = {lgd_lab{:} 'RT'};
        end
        
        % Plotting parameters
        ax.Title.String  = [roi_id ': ' event_lab{sr_ix}];
        ax.Box           = 'off';
        ax.YLim          = ylims;
        ax.YTick         = yticks;
        ax.YLabel.String = ylab;
        ax.XLim          = [0,size(stat{sr_ix}.time,2)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(stat{sr_ix}.time,2);
        ax.XTickLabel    = x_tick_lab;
        ax.XLabel.String = 'Time (s)';
        %             ax(2).YColor        = 'k';
        %             legend(main_lines,lgd_lab{:},'Location',plt_vars.legend_loc);
    end
    
    % Save figure
    if any_plot
        if save_fig
            fig_filename = [fig_dir fig_name '.' fig_filetype];
            fprintf('Saving %s\n',fig_filename);
            saveas(gcf,fig_filename);
            %eval(['export_fig ' fig_filename]);
        end
    else
        close(f);
    end
end

%% Save out list of channels with significant differences
% sig_report_filename = [fig_dir 'sig_ch_list.txt'];
% sig_report = fopen(sig_report_filename,'a');
% fprintf(sig_report,'%s: %s for %s - %s\n',SBJ,stat_id,an_id_s,an_id_r);
% for grp_ix = 1:numel(grp_lab)+1
%     for sr_ix = 1:numel(event_lab)
%         fprintf(sig_report,'%s, %s :\n',lgd_lab{grp_ix},event_lab{sr_ix});
%         if ~isempty(sig_ch{grp_ix,sr_ix})
%             fprintf(sig_report,'\t%s\n',strjoin(sig_ch{grp_ix,sr_ix},','));
%         else
%             fprintf(sig_report,'\n');
%         end
%     end
% end
% fclose(sig_report);

end
