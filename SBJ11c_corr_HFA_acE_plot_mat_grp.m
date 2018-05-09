function SBJ11c_corr_HFA_acE_plot_mat_grp(SBJs,pipeline_id,stat_id,an_id_main,plt_id,...
    plot_scnd_mtx,save_fig,fig_vis,fig_filetype)
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
actv_ix = zeros(size(einfo,1),1);
for e_ix = 1:numel(hfa.label)
    if any(strcmp(einfo(:,1),hfa.label{e_ix}))
        actv_ix(strcmp(einfo(:,1),hfa.label{e_ix})) = 1;
    end
end
einfo = einfo(logical(actv_ix),:);
einfo = sortrows(einfo,sort_vec);
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
major_lab = unique(einfo(:,sort_vec(1)));
maj_line_pos = NaN(length(major_lab),1);
for lab_ix = 2:length(major_lab)      % Skip #1 because that'll just be line 1 of the first section
    tmp_idx = find(ismember(einfo(:,sort_vec(1)),major_lab(lab_ix)));
    maj_line_pos(lab_ix) = tmp_idx(1);
end
einfo(:,8)={0};
einfo(maj_line_pos(2:end),8)={1};

minor_lab = unique(einfo(:,sort_vec(2)));     % these are labels with "sort1-sort2"
min_line_pos = NaN(length(minor_lab),1);
for lab_ix=2:length(minor_lab)
    tmp_idx = find(ismember(einfo(:,sort_vec(2)),minor_lab(lab_ix)));
    min_line_pos(lab_ix) = tmp_idx(1);
end
einfo(:,9)={0};
einfo(min_line_pos(2:end),9)={1};

%% Averaged within main ROI Matrix
fig_name = strcat(SBJ,'_corr_HFA_acE_mat_avg_main_',plt_id);
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
hold on;

rois = unique(einfo(:,sort_vec(1)));% numel(unique(einfo(:,sort_vec(2))))];
roi_corr_cnt = zeros(numel(rois));
roi_elec_cnt = zeros(numel(rois));
for pair_ix = 1:size(pairs,1)
    roi_ix = sort([find(strcmp(rois,einfo{pairs(pair_ix,1),sort_vec(1)})),...
                    find(strcmp(rois,einfo{pairs(pair_ix,2),sort_vec(1)}))]);
    roi_elec_cnt(roi_ix(1),roi_ix(2)) = roi_elec_cnt(roi_ix(1),roi_ix(2))+1;
    if qvals(pair_ix)<=plt_vars.sig_cut
        roi_corr_cnt(roi_ix(1),roi_ix(2)) = roi_corr_cnt(roi_ix(1),roi_ix(2))+1;
    end
end
imagesc(roi_corr_cnt./roi_elec_cnt); axis xy;

ax_tick = 1:numel(rois);
set(gca,'XLim',[0.5 numel(rois)+0.5],'XTick',ax_tick,'XTickLabel',rois,'FontSize',plt_vars.font_sz);
set(gca,'YLim',[0.5 numel(rois)+0.5],'YTick',ax_tick,'YTickLabel',rois,'FontSize',plt_vars.font_sz);
xlabel(einfo_fields{sort_vec(1)});
ylabel(einfo_fields{sort_vec(1)});

subplot_name = strcat(stat_vars.event_lab,': Proportion sig by ',einfo_fields{sort_vec(1)});
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
    results_dir = 'grp'['~/PRJ_Stroop/results/HFA/',SBJ,'/',stat_id,'/',an_id_main,'/'];
    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end
    saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
end

%% Averaged within secondary ROI Matrix
if plot_scnd_mtx
    fig_name = strcat(SBJ,'_corr_HFA_acE_mat_avg_scnd_',plt_id);
    figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
    hold on;
    
    rois = unique(einfo(:,sort_vec(2)));
    roi_corr_cnt = zeros(numel(rois));
    roi_elec_cnt = zeros(numel(rois));
    for pair_ix = 1:size(pairs,1)
        roi_ix = sort([find(strcmp(rois,einfo{pairs(pair_ix,1),sort_vec(2)})),...
                    find(strcmp(rois,einfo{pairs(pair_ix,2),sort_vec(2)}))]);   % use sort to make sure all on one side of matrix
        roi_elec_cnt(roi_ix(1),roi_ix(2)) = roi_elec_cnt(roi_ix(1),roi_ix(2))+1;
        if qvals(pair_ix)<=plt_vars.sig_cut
            roi_corr_cnt(roi_ix(1),roi_ix(2)) = roi_corr_cnt(roi_ix(1),roi_ix(2))+1;
        end
    end
    imagesc(roi_corr_cnt./roi_elec_cnt); axis xy;
    
    ax_tick = 1:numel(rois);
    set(gca,'XLim',[0.5 numel(rois)+0.5],'XTick',ax_tick,'XTickLabel',rois,'FontSize',plt_vars.font_sz);
    set(gca,'YLim',[0.5 numel(rois)+0.5],'YTick',ax_tick,'YTickLabel',rois,'FontSize',plt_vars.font_sz);
    xlabel(einfo_fields{sort_vec(2)});
    ylabel(einfo_fields{sort_vec(2)});
    
    subplot_name = strcat(stat_vars.event_lab,': Proportion sig by ',einfo_fields{sort_vec(2)});
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


end