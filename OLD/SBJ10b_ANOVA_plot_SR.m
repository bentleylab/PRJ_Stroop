function SBJ10b_ANOVA_plot_SR(SBJ,stat_id,an_id_s,an_id_r,plt_id,save_fig,fig_vis)
% Plots ANOVA results
% clear all; %close all;

fig_filetype = 'png';
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);

[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
event_lab = {'stim', 'resp'};

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

f_name_s = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id_s '.mat'];
f_name_r = [SBJ_vars.dirs.proc SBJ '_ANOVA_ROI_' stat_id '_' an_id_r '.mat'];
tmp = load(f_name_s,'w2'); w2{1} = tmp.w2;
tmp = load(f_name_r,'w2'); w2{2} = tmp.w2;
tmp = load(f_name_s,'hfa'); hfa{1} = tmp.hfa;
tmp = load(f_name_r,'hfa'); hfa{2} = tmp.hfa;
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(hfa{1}.time)-1)/(hfa{1}.time(end)-hfa{1}.time(1));

%% Prep Data
% Trim data to plotting epoch
cfg_trim = [];
cfg_trim.latency = plt_vars.plt_lim_S;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
cfg_trim.latency = plt_vars.plt_lim_R;
hfa{2} = ft_selectdata(cfg_trim,hfa{2});

win_lim    = fn_sliding_window_lim(squeeze(hfa{1}.powspctrm(1,1,1,:)),win_len,win_step);
win_center = round(mean(win_lim,2));

% % Compute mean RT per condition
% RTs = round(1000*trial_info.response_time); % converts sec to ms
% for cond_ix = 1:numel(grp_lab{congr_ix})
%     RT_means{cond_ix} = mean(RTs(fn_condition_index([grp_lab{congr_ix}{cond_ix}], trial_info.condition_n)==1));
%     % Add in the baseline offset to plot correctly
%     RT_means{cond_ix} = RT_means{cond_ix}-plt_vars.plt_lim_S(1)*1000;
% end

%% Plot Results
fig_dir = ['/home/knight/hoycw/PRJ_Stroop/results/HFA/' SBJ '/' stat_id '/SR/' an_id_s '-' an_id_r '/'];
if ~exist(fig_dir,'dir')
    mkdir(fig_dir);
end

% Find plot limits
max_w2 = max([max(max(max(w2{1}.trial))) max(max(max(w2{1}.trial)))]);
min_w2 = min([min(min(min(w2{1}.trial))) min(min(min(w2{1}.trial)))]);
ylim_fudge = (max_w2-min_w2)*plt_vars.ylim_fudge;
ylims  = [min_w2-ylim_fudge max_w2+ylim_fudge];
y_sig = zeros(size(grp_lab));
y_sig(1) = mean([min_w2,max_w2]);
for grp_ix = 2:numel(grp_lab)
    y_sig(grp_ix) = y_sig(grp_ix-1)+ylim_fudge;
end

% Create a figure for each channel
sig_ch = cell(size(grp_lab));
for ch_ix = 1:numel(hfa{1}.label)    
    fig_name = [SBJ '_ANOVA_' stat_id '_SR_' hfa{1}.label{ch_ix}];    
    figure('Name',fig_name,'units','normalized',...
        'outerposition',[0 0 1 0.5],'Visible',fig_vis); %twice as wide for the double plot
    
    for set_ix = 1:2
        subplot(1,2,set_ix);
        hold on;
        ax = gca;
        ax.Title.String  = [hfa{set_ix}.label{ch_ix} ':' event_lab{set_ix}];
        lgd_lab          = grp_lab;
        
        % Compute explained variance and standard error
        var_exp = squeeze(w2{set_ix}.trial(:,ch_ix,:));
        error('I dont think computing error bars based on permutation distributions is valid! need repeats across channels or something...');
        var = NaN([numel(grp_lab) size(var_exp,2)]);
        for grp_ix = 1:numel(grp_lab)
            var(grp_ix,:) = squeeze(std(w2{set_ix}.boot(grp_ix,ch_ix,:,:),[],4))./...
                sqrt(size(w2{set_ix}.boot,4))';
        end
        
        % Plot beta from RT
        main_lines = [];
%         if rt_confound
%             beta_line = plot(win_center, squeeze(hfa{set_ix}.beta(:,ch_ix,:,win_center)),...
%                 'Color',[127,201,127]./256);
%             main_lines = [main_lines beta_line];
%             lgd_lab = {'RT' lgd_lab{:}};
%         end
        
        % Plot var_exp
        ebars = {};
        for grp_ix = 1:length(grp_lab)
            ebars{grp_ix} = shadedErrorBar(win_center, var_exp(grp_ix,:), var(grp_ix,:),...
                {'Color',[grp_colors{grp_ix}],...
                'LineStyle',grp_style{grp_ix}},plt_vars.errbar_alpha);
            main_lines = [main_lines ebars{grp_ix}.mainLine];
        end
        
        % Plot event
        if strcmp(event_lab{set_ix},'stim')
            x_tick_lab       = plt_vars.plt_lim_S(1):plt_vars.x_step_sz:plt_vars.plt_lim_S(2);
        else
            x_tick_lab       = plt_vars.plt_lim_R(1):plt_vars.x_step_sz:plt_vars.plt_lim_R(2);
            % Plot Response Marker
            event_line = line([find(hfa{set_ix}.time==0) find(hfa{set_ix}.time==0)],ylims,...
                'LineWidth',plt_vars.evnt_width, 'Color',...
                plt_vars.evnt_color, 'LineStyle',plt_vars.evnt_style);
            main_lines = [main_lines event_line];
            lgd_lab = {lgd_lab{:} 'RT'};
        end
        
        % Plot significant time periods
        for grp_ix = 1:numel(grp_lab)
            if any(w2{set_ix}.pval(grp_ix,ch_ix,:)<0.05)
                % Track channels with significant var_exp
                if ~iscell(sig_ch{grp_ix})
                    sig_ch{grp_ix} = {};
                end
                sig_ch{grp_ix} = {sig_ch{grp_ix}{:} hfa{set_ix}.label{ch_ix}};
                
                % Find significant periods
                if strcmp(plt_vars.sig_type,'scatter')
                    sig_times = win_center(squeeze(w2{set_ix}.pval(grp_ix,ch_ix,:))<0.05);
                    scatter(sig_times,repmat(y_sig(grp_ix),size(sig_times)),...
                        plt_vars.sig_scat_size,grp_colors{grp_ix},plt_vars.sig_scat_mrkr);
                elseif strcmp(plt_vars.sig_type,'patch')
                    sig_chunks = fn_find_chunks(squeeze(w2{set_ix}.pval(grp_ix,ch_ix,:))<0.05);
                    sig_chunks(w2{set_ix}.pval(grp_ix,ch_ix,sig_chunks(:,1))>0.05,:) = [];
                    
                    % Plot Significance Shading
                    fprintf('%s %s (%s) -- %i SIGNIFICANT CLUSTERS FOUND...\n',...
                        hfa{set_ix}.label{ch_ix},grp_lab{grp_ix},event_lab{set_ix},size(sig_chunks,1));
                    for sig_ix = 1:size(sig_chunks,1)
                        sig_times = win_lim(sig_chunks(sig_ix,:),:);
                        patch([sig_times(1,1) sig_times(1,1) sig_times(2,2) sig_times(2,2)], ...
                            [ylims(1) ylims(2) ylims(2) ylims(1)],...
                            grp_colors{grp_ix},'FaceAlpha',plt_vars.sig_alpha);
                    end
                end
            else
                fprintf('%s %s (%s) -- NO SIGNIFICANT CLUSTERS FOUND...\n',...
                    hfa{set_ix}.label{ch_ix},grp_lab{grp_ix},event_lab{set_ix});
            end
        end
        
        % Plotting parameters
        legend(main_lines,lgd_lab{:},'Location',plt_vars.legend_loc);
        ax.YLim          = ylims;
        ax.YLabel.String = '% Variance Explained';
        ax.XLim          = [0,size(hfa{set_ix}.time,2)];
        ax.XTick         = 0:plt_vars.x_step_sz*sample_rate:size(hfa{set_ix}.time,2);
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
end

% % Save out list of channels with significant differences
% sig_report_filename = [fig_dir 'sig_ch_list.txt'];
% sig_report = fopen(sig_report_filename,'a');
% fprintf(sig_report,'%s - %s\n',an_id_s,an_id_r);
% fprintf(sig_report,'%s\n',[sig_ch{:}]);
% fclose(sig_report);

end
