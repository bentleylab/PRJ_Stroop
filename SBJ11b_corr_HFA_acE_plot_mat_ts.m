function SBJ11b_corr_HFA_acE_plot_mat_ts(SBJ,pipeline_id,stat_id,an_id_main,plt_id,...
    plot_ts,plot_scnd_mtx,save_fig,fig_vis,fig_filetype)
% Build connectivity matrix based on HFA correlations, plot as matrix
%   non-parametric stats via circular shift of trial time series

if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

if strcmp(stat_vars.event_lab,'stim')
    [cond_lab, cond_colors, ~] = fn_condition_label_styles(stat_vars.conditions);
    
    % Load RTs
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');
    % Compute mean RT per condition
    RTs = round(1000*trial_info.response_time); % converts sec to ms
    for cond_ix = 1:numel(cond_lab)
        RT_means{cond_ix} = mean(RTs(fn_condition_index([cond_lab{cond_ix}], trial_info.condition_n)==1));
        % Add in the baseline offset to plot correctly
        RT_means{cond_ix} = RT_means{cond_ix}-stat_vars.stat_lim(1)*1000;
    end
end

% Load data
corr_filename = [SBJ_vars.dirs.proc SBJ '_' stat_id '_' an_id_main '.mat'];
load(corr_filename);
sig_pairs = find(qvals<=plt_vars.sig_cut);
corr_cut_pairs = find(abs(corr_vals)>=plt_vars.corr_cut);
plot_pairs = intersect(sig_pairs,corr_cut_pairs);

%% Load ROI and GM/WM info
einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
load(einfo_filename);
% Select active electrodes
if stat_vars.only_actv_ch

actv_ix = zeros(size(einfo,1),1);
for e_ix = 1:numel(hfa.label)
    if any(strcmp(einfo(:,1),hfa.label{e_ix}))
        actv_ix(strcmp(einfo(:,1),hfa.label{e_ix})) = 1;
    end
end
einfo = einfo(logical(actv_ix),:);
einfo = sortrows(einfo,plt_vars.sort_vec);
einfo_fields = {'label','ROI','gROI','ROI2','tissue','GM weight','Out'};
% Electrode Info Table:
%   label- name of electrode
%   ROI- specific region
%   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%   ROI2- specific region of second electrode
%   tissue- primary tissue type
%   GM weight- percentage of electrode pair in GM
%   Out- 0/1 flag for whether this is partially out of the brain

% Major and Minor divisions in einfo for plotting lines to divide matrix
major_lab = unique(einfo(:,plt_vars.sort_vec(1)));
maj_line_pos = NaN(length(major_lab),1);
for lab_ix = 2:length(major_lab)      % Skip #1 because that'll just be line 1 of the first section
    tmp_idx = find(ismember(einfo(:,plt_vars.sort_vec(1)),major_lab(lab_ix)));
    maj_line_pos(lab_ix) = tmp_idx(1);
end
einfo(:,8)={0};
einfo(maj_line_pos(2:end),8)={1};

minor_lab = unique(einfo(:,plt_vars.sort_vec(2)));     % these are labels with "sort1-sort2"
min_line_pos = NaN(length(minor_lab),1);
for lab_ix=2:length(minor_lab)
    tmp_idx = find(ismember(einfo(:,plt_vars.sort_vec(2)),minor_lab(lab_ix)));
    min_line_pos(lab_ix) = tmp_idx(1);
end
einfo(:,9)={0};
einfo(min_line_pos(2:end),9)={1};

%% Plot Correlation Matrix with significance
fprintf('=========================== Plotting Correlation Matrix ===========================\n');
% Thin out the labels
if any(plt_vars.thin_labels)
    [x_lab,y_lab] = fn_thin_labels(einfo,plt_vars.sort_vec,plt_vars.thin_labels);
end

fig_name = strcat(SBJ,'_corr_HFA_acE_mat_',plt_id);
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
% colormap(color_map);
hold on;

% Draw matrix with significant dots
imagesc(corr_mat); axis xy;
scatter(pairs(plot_pairs,1),pairs(plot_pairs,2),...
    plt_vars.sig_dot_size,plt_vars.sig_dot_color,'filled');
if plt_vars.double_sig_dot
    scatter(pairs(plot_pairs,1),pairs(plot_pairs,2),...
        floor(plt_vars.sig_dot_size2),plt_vars.sig_dot_color2,'filled');
end

% Draw Labels
% set(gcf,'OuterPosition',plt_vars.subplot_pos(sr_ix,:),'PlotBoxAspectRatio',[1 1 1]);
ax_tick = 1:numel(hfa.label);
set(gca,'XLim',[0.5 numel(hfa.label)+0.5],'XTick',ax_tick,'XTickLabel',x_lab,'FontSize',plt_vars.font_sz);
set(gca,'YLim',[0.5 numel(hfa.label)+0.5],'YTick',ax_tick,'YTickLabel',y_lab,'FontSize',plt_vars.font_sz);
xlabel(einfo_fields{plt_vars.sort_vec(1)});
ylabel(einfo_fields{plt_vars.sort_vec(2)});

% Set color scale
max_corr = max(abs(corr_vals));
caxis([-max_corr max_corr]);
colorbar;
%                 func_xticklabel_rotate(axTick,rotAng,xLabsubels,'FontSize',fontSize);

% Get divisions between major zones of matrix
if plt_vars.plt_maj_div
    maj_divs = find([einfo{:,8}]==1)-0.5;
    for line_ix = 1:length(maj_divs)
        line(xlim,[maj_divs(line_ix) maj_divs(line_ix)],'Color','k','LineStyle','-','LineWidth',plt_vars.div_line_width);
        line([maj_divs(line_ix) maj_divs(line_ix)],ylim,'Color','k','LineStyle','-','LineWidth',plt_vars.div_line_width);
    end
end
if plt_vars.plt_min_div
    min_divs = find(([einfo{:,9}]==1)&([einfo{:,8}]~=1))-0.5;
    for line_ix = 1:length(min_divs)
        line(xlim,[min_divs(line_ix) min_divs(line_ix)],'Color','k','LineStyle','--','LineWidth',plt_vars.div_line_width);
        line([min_divs(line_ix) min_divs(line_ix)],ylim,'Color','k','LineStyle','--','LineWidth',plt_vars.div_line_width);
    end
end

subplot_name = strcat(stat_vars.event_lab,': ',num2str(numel(plot_pairs)),' sig');
title(subplot_name,'interpreter','none');
colorbar

% Save Fig
if save_fig
    results_dir = ['~/PRJ_Stroop/results/HFA/',SBJ,'/',stat_id,'/',an_id_main,'/'];
    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end
    saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
end

%% Averaged within main ROI Matrix
fig_name = strcat(SBJ,'_corr_HFA_acE_mat_avg_main_',plt_id);
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
hold on;

rois = unique(einfo(:,plt_vars.sort_vec(1)));% numel(unique(einfo(:,plt_vars.sort_vec(2))))];
roi_corr_cnt = zeros(numel(rois));
roi_elec_cnt = zeros(numel(rois));
for pair_ix = 1:size(pairs,1)
    roi_ix = sort([find(strcmp(rois,einfo{pairs(pair_ix,1),plt_vars.sort_vec(1)})),...
                    find(strcmp(rois,einfo{pairs(pair_ix,2),plt_vars.sort_vec(1)}))]);
    roi_elec_cnt(roi_ix(1),roi_ix(2)) = roi_elec_cnt(roi_ix(1),roi_ix(2))+1;
    if qvals(pair_ix)<=plt_vars.sig_cut
        roi_corr_cnt(roi_ix(1),roi_ix(2)) = roi_corr_cnt(roi_ix(1),roi_ix(2))+1;
    end
end
imagesc(roi_corr_cnt./roi_elec_cnt); axis xy;

ax_tick = 1:numel(rois);
set(gca,'XLim',[0.5 numel(rois)+0.5],'XTick',ax_tick,'XTickLabel',rois,'FontSize',plt_vars.font_sz);
set(gca,'YLim',[0.5 numel(rois)+0.5],'YTick',ax_tick,'YTickLabel',rois,'FontSize',plt_vars.font_sz);
xlabel(einfo_fields{plt_vars.sort_vec(1)});
ylabel(einfo_fields{plt_vars.sort_vec(1)});

subplot_name = strcat(stat_vars.event_lab,': Proportion sig by ',einfo_fields{plt_vars.sort_vec(1)});
title(subplot_name,'interpreter','none');
colorbar;
caxis([0 1]);


% % One color bar for whole figure
% c=axes('OuterPosition', [0.5 0 0.5 1], 'Visible', 'off');
% tickDiv=(colors(2)-colors(1))/5;
% colorbar('YTick',[colors(1):tickDiv:colors(2)]);
% caxis(colors);

%         suptitle(strcat(SBJ,',',task,': ',fband_pairs{fband},' (combined) ',sortID,'_cscale.',colorscale));
% Save plot
if save_fig
    saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
end

%% Averaged within secondary ROI Matrix
if plot_scnd_mtx
    fig_name = strcat(SBJ,'_corr_HFA_acE_mat_avg_scnd_',plt_id);
    figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
    hold on;
    
    rois = unique(einfo(:,plt_vars.sort_vec(2)));
    roi_corr_cnt = zeros(numel(rois));
    roi_elec_cnt = zeros(numel(rois));
    for pair_ix = 1:size(pairs,1)
        roi_ix = sort([find(strcmp(rois,einfo{pairs(pair_ix,1),plt_vars.sort_vec(2)})),...
                    find(strcmp(rois,einfo{pairs(pair_ix,2),plt_vars.sort_vec(2)}))]);   % use sort to make sure all on one side of matrix
        roi_elec_cnt(roi_ix(1),roi_ix(2)) = roi_elec_cnt(roi_ix(1),roi_ix(2))+1;
        if qvals(pair_ix)<=plt_vars.sig_cut
            roi_corr_cnt(roi_ix(1),roi_ix(2)) = roi_corr_cnt(roi_ix(1),roi_ix(2))+1;
        end
    end
    imagesc(roi_corr_cnt./roi_elec_cnt); axis xy;
    
    ax_tick = 1:numel(rois);
    set(gca,'XLim',[0.5 numel(rois)+0.5],'XTick',ax_tick,'XTickLabel',rois,'FontSize',plt_vars.font_sz);
    set(gca,'YLim',[0.5 numel(rois)+0.5],'YTick',ax_tick,'YTickLabel',rois,'FontSize',plt_vars.font_sz);
    xlabel(einfo_fields{plt_vars.sort_vec(2)});
    ylabel(einfo_fields{plt_vars.sort_vec(2)});
    
    subplot_name = strcat(stat_vars.event_lab,': Proportion sig by ',einfo_fields{plt_vars.sort_vec(2)});
    title(subplot_name,'interpreter','none');
    colorbar;
    caxis([0 1]);
    
    
    % % One color bar for whole figure
    % c=axes('OuterPosition', [0.5 0 0.5 1], 'Visible', 'off');
    % tickDiv=(colors(2)-colors(1))/5;
    % colorbar('YTick',[colors(1):tickDiv:colors(2)]);
    % caxis(colors);
    
    %         suptitle(strcat(SBJ,',',task,': ',fband_pairs{fband},' (combined) ',sortID,'_cscale.',colorscale));
    % Save plot
    if save_fig
        saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
    end
end

%% Plot Distribution of correlation values
fig_name = strcat(SBJ,'_corr_HFA_acE_hist_',plt_id);
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
hold on;

%     ksdensity(corr_vals);
histogram(corr_vals,plt_vars.hist_bins,'Normalization','probability');
%     y_vals = ylim;
%     heights = rand(sum(qvals<plt_vars.sig_cut),1)*plt_vars.y_jitter*diff(y_vals)+mean(y_vals);
%     scatter(corr_vals(plot_pairs),heights);

%     subplot(2,2,sr_ix+2); hold on;
scatter(corr_vals(setdiff(1:numel(corr_vals),plot_pairs)),qvals(setdiff(1:numel(qvals),plot_pairs)),'k');
scatter(corr_vals(plot_pairs),qvals(plot_pairs),'r');
xlabel('Correlation');
ylabel('qval / normalized count');
title(stat_vars.event_lab);

% Compare Stim and Response locked
% subplot(1,3,3); hold on;
% sig_pairs1 = find(qvals{1}<=plt_vars.sig_cut);
% sig_pairs2 = find(qvals{2}<=plt_vars.sig_cut);
% sig_pairs = union(sig_pairs1,sig_pairs2);
% ns_pairs = setdiff(1:numel(qvals{1}),sig_pairs);
% scatter(corr_vals{1}(ns_pairs),corr_vals{2}(ns_pairs),'k');
% scatter(corr_vals{1}(sig_pairs1),corr_vals{2}(sig_pairs1),'r','Marker','x');
% scatter(corr_vals{1}(sig_pairs2),corr_vals{2}(sig_pairs2),'r','Marker','+');
% scatter(corr_vals{1}(intersect(sig_pairs1,sig_pairs2)),corr_vals{2}(intersect(sig_pairs1,sig_pairs2)),'c','Marker','*');
% xlabel([stat_vars.event_lab{1} ' Correlations']);
% ylabel([stat_vars.event_lab{2} ' Correlations']);
% legend('not sig',[stat_vars.event_lab{1} '-sig'],[stat_vars.event_lab{2} '-sig'],'both sig');
% title([stat_vars.event_lab{1} ' vs. ' stat_vars.event_lab{2}]);

% Save plot
if save_fig
    saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
end

%% Plot Time Series of Correlated Electrode Pairs
if plot_ts
    error('sig_pairs indexing isnt right and the sort must be handled!');
    sig_pairs = find(sort(qvals)<=plt_vars.sig_cut);
    if plot_ts==1
        plot_ts = numel(sig_pairs);
    else
        sig_pairs = sig_pairs(1:plot_ts);
    end
    fprintf('=========================== Plotting Time Series: n = %i/%i ===========================\n',...
        plot_ts,numel(plot_pairs));
    if plot_ts~=numel(sig_pairs)
        sig_pairs = sig_pairs(1:plot_ts);
    end
    for pair_ix = 1:numel(sig_pairs)
        fprintf('%i..',pair_ix);
        if mod(pair_ix,20)==0
            fprintf('\n');
        end
        epair_ix(1) = pairs(sig_pairs(pair_ix),1);
        epair_ix(2) = pairs(sig_pairs(pair_ix),2);
        fig_name = strcat(SBJ,'_corr_HFA_acE_SR_ts_',plt_id,'_',...
            hfa.label{epair_ix(1)},'-',hfa.label{epair_ix(2)});
        figure('Name',fig_name,'Visible','on','units','normalized','outerposition',plt_vars.fig_dim);
        plot_info.fig        = gcf;
        plot_info.x_step     = plt_vars.x_step_sz*proc_vars.resample_freq;
        plot_info.legend_loc = plt_vars.legend_loc;
        plot_info.sig_alpha  = plt_vars.errbar_alpha;
        plot_info.sig_color  = plt_vars.errbar_color;
        % Condition plotting params
        cond_info.name       = {hfa.label{epair_ix(1)} hfa.label{epair_ix(2)}};
        cond_info.style      = plt_vars.main_style;
        cond_info.color      = plt_vars.main_colors;
        cond_info.alpha      = repmat(plt_vars.errbar_alpha,[1 numel(cond_lab)]);
        % colormap(color_map);
        
        % Average HFA time series and compute SEM
        hfa_ts = NaN([2 size(hfa.powspctrm,4)]);
        sem_ts = hfa_ts;
        for e_ix = 1:2
            hfa_ts(e_ix,:) = squeeze(mean(hfa.powspctrm(:,epair_ix(e_ix),1,:),1))';
            sem_ts(e_ix,:) = squeeze(std(hfa.powspctrm(:,epair_ix(e_ix),1,:),[],1)./...
                sqrt(size(hfa.powspctrm,1)))';
        end
        
        % Find plot limits
        max_y = max(hfa_ts(:));
        min_y = min(hfa_ts(:));
        ylim_fudge = (max_y-min_y)*plt_vars.ylim_fudge;
        ylims  = [min_y-ylim_fudge max_y+ylim_fudge];
        yticks = 0:plt_vars.y_tick_int:ylims(2);
        
        %             subplot(1,2,sr_ix);
        hold on;
        
        plot_info.ax     = gca;
        plot_info.title  = [hfa.label{epair_ix(1)} '+' hfa.label{epair_ix(2)}...
            ': r=' num2str(corr_vals(sig_pairs(pair_ix))) ',q=' num2str(qvals(sig_pairs(pair_ix)))];
        plot_info.legend = plt_vars.legend;
        if strcmp(stat_vars.event_lab,'stim')
            plot_info.x_lab = stat_vars.stat_lim(1):plt_vars.x_step_sz:stat_vars.stat_lim(2);
            % Stimulus plotting params
            event_info.name  = {stat_vars.event_lab, cond_lab{:}};
            event_info.color = {[0 0 0], cond_colors{:}};
            event_info.width = repmat(plt_vars.evnt_width,[1 numel(event_info.name)]);
            event_info.style = repmat({plt_vars.evnt_style},[1 numel(event_info.name)]);
            event_info.time  = [-stat_vars.stat_lim(1)*proc_vars.resample_freq, RT_means{:}];
        else
            plot_info.x_lab = stat_vars.stat_lim(1):plt_vars.x_step_sz:stat_vars.stat_lim(2);
            % Stimulus plotting params
            event_info.name  = {stat_vars.event_lab};
            event_info.width = plt_vars.evnt_width;
            event_info.color = {plt_vars.evnt_color};
            event_info.style = {plt_vars.evnt_style};
            event_info.time  = -stat_vars.stat_lim(1)*proc_vars.resample_freq;
        end
        
        % Plot time series
        fn_plot_ts_error_bar(plot_info,hfa_ts,sem_ts,event_info,cond_info);
        %         set(gca,'YLim',[min_y max_y]);
        
        % Save Figure
        if save_fig
            saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
        end
        %     close(gcf);
    end
    fprintf('\n');
end
end
