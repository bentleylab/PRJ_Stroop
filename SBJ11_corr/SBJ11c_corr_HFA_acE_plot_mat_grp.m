function SBJ11c_corr_HFA_acE_plot_mat_grp(SBJs,pipeline_id,stat_id,an_id_main,roi_id,atlas_id,plt_id,...
    plot_scnd_mtx,save_fig,fig_vis,fig_filetype)
% Build connectivity matrix based on HFA correlations, plot as matrix
%   non-parametric stats via circular shift of trial time series

if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

if strcmp(stat_vars.event_lab,'stim')
    [cond_lab, cond_colors, ~] = fn_condition_label_styles(stat_vars.conditions);
end

[rois, roi_colors] = fn_roi_label_styles(roi_id);
if plt_vars.exclude_FWM
    rois(strcmp(rois,'FWM')) = [];
end
if plt_vars.exclude_OUT
    rois(strcmp(rois,'OUT')) = [];
end

%% Load and compile data
roi_corr_cnt = zeros(numel(SBJs), numel(rois), numel(rois));
roi_elec_cnt = zeros(numel(SBJs), numel(rois), numel(rois));
for sbj_ix = 1:numel(SBJs)
    SBJ = SBJs{sbj_ix};
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    
    % Load data
    % No need to load elec, it was saved in SBJ11a with the correct sorting
    corr_filename = [SBJ_vars.dirs.proc SBJ '_' stat_id '_' an_id_main '_' roi_id '_' atlas_id '.mat'];
    load(corr_filename);
    
    % Get ROIs of each electrode pair
    pair_rois = zeros(size(pairs));
    for pair_ix = 1:size(pairs,1)
        pair_rois(pair_ix,:) = [find(strcmp(elec.roi{pairs(pair_ix,1)},rois)) ...
            find(strcmp(elec.roi{pairs(pair_ix,2)},rois))];
    end
    
    sig_pairs = find(qvals<=plt_vars.sig_cut);
    corr_cut_pairs = find(abs(corr_vals)>=plt_vars.corr_cut);
    good_pairs = intersect(sig_pairs,corr_cut_pairs);
    
    % Count sig and total pairs
    for pair_ix = 1:numel(good_pairs)
        % Sort (now min/max) roi_ix so that all results are aggregated on one side of
        % diagonal (it doesn't matter which elec was in which ROI)
        roi_ix = sort(pair_rois(good_pairs(pair_ix),:));
        %         roi_ix = sort([find(strcmp(elec.roi{pairs(pair_ix,1)},rois)) ...
        %             find(strcmp(elec.roi{pairs(pair_ix,2)},rois))]);
        roi_elec_cnt(sbj_ix,roi_ix(1),roi_ix(2)) = roi_elec_cnt(sbj_ix,roi_ix(1),roi_ix(2))+1;
        if qvals(pair_ix)<=plt_vars.sig_cut
            roi_corr_cnt(sbj_ix,roi_ix(1),roi_ix(2)) = roi_corr_cnt(sbj_ix,roi_ix(1),roi_ix(2))+1;
        end
    end
    clear SBJ corr_mat corr_vals hfa pairs qvals elec
end

%% Averaged within main ROI Matrix
fig_name = strcat('GRP_corr_HFA_acE_mat_avg_main_',roi_id,'_',atlas_id,plt_id);
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
hold on;

imagesc(squeeze(mean(roi_corr_cnt,1))./squeeze(mean(roi_elec_cnt,1))); axis xy;

ax_tick = 1:numel(rois);
set(gca,'XLim',[0.5 numel(rois)+0.5],'XTick',ax_tick,'XTickLabel',rois,'FontSize',plt_vars.font_sz);
set(gca,'YLim',[0.5 numel(rois)+0.5],'YTick',ax_tick,'YTickLabel',rois,'FontSize',plt_vars.font_sz);
xlabel(roi_id);
ylabel(roi_id);

subplot_name = strcat(stat_vars.event_lab,': Proportion sig by ',roi_id);
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
    results_dir = [root_dir 'PRJ_Stroop/results/HFA/GRP_',stat_id,'/',an_id_main,'/'];
    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end
    saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
end

%% Averaged within secondary ROI Matrix
% Not needed with atlas_id assignments

% if plot_scnd_mtx
%     fig_name = strcat(SBJ,'_corr_HFA_acE_mat_avg_scnd_',plt_id);
%     figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
%     hold on;
%     
%     rois = unique(einfo(:,sort_vec(2)));
%     roi_corr_cnt = zeros(numel(rois));
%     roi_elec_cnt = zeros(numel(rois));
%     for pair_ix = 1:size(pairs,1)
%         roi_ix = sort([find(strcmp(rois,einfo{pairs(pair_ix,1),sort_vec(2)})),...
%                     find(strcmp(rois,einfo{pairs(pair_ix,2),sort_vec(2)}))]);   % use sort to make sure all on one side of matrix
%         roi_elec_cnt(roi_ix(1),roi_ix(2)) = roi_elec_cnt(roi_ix(1),roi_ix(2))+1;
%         if qvals(pair_ix)<=plt_vars.sig_cut
%             roi_corr_cnt(roi_ix(1),roi_ix(2)) = roi_corr_cnt(roi_ix(1),roi_ix(2))+1;
%         end
%     end
%     imagesc(roi_corr_cnt./roi_elec_cnt); axis xy;
%     
%     ax_tick = 1:numel(rois);
%     set(gca,'XLim',[0.5 numel(rois)+0.5],'XTick',ax_tick,'XTickLabel',rois,'FontSize',plt_vars.font_sz);
%     set(gca,'YLim',[0.5 numel(rois)+0.5],'YTick',ax_tick,'YTickLabel',rois,'FontSize',plt_vars.font_sz);
%     xlabel(einfo_fields{sort_vec(2)});
%     ylabel(einfo_fields{sort_vec(2)});
%     
%     subplot_name = strcat(stat_vars.event_lab,': Proportion sig by ',einfo_fields{sort_vec(2)});
%     title(subplot_name,'interpreter','none');
%     colorbar;
%     caxis([0 1]);
%     
%     
%     % % One color bar for whole figure
%     % c=axes('OuterPosition', [0.5 0 0.5 1], 'Visible', 'off');
%     % tickDiv=(colors(2)-colors(1))/5;
%     % colorbar('YTick',[colors(1):tickDiv:colors(2)]);
%     % caxis(colors);
%     
%     %         suptitle(strcat(SBJ,',',task,': ',fband_pairs{fband},' (combined) ',sortID,'_cscale.',colorscale));
%     % Save plot
%     if save_fig
%         saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
%     end
% end


end