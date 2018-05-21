function SBJ11c_corr_HFA_cond_acE_plot_grp_mat(SBJs,pipeline_id,stat_id,an_id_main,plt_id,...
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
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

[cond_lab, cond_colors, ~] = fn_condition_label_styles(stat_vars.conditions);

einfo_fields = {'label','ROI','gROI','ROI2','tissue','GM weight','Out'};
load('/home/knight/hoycw/PRJ_Stroop/data/full_roi_lists.mat');
eval(['rois = all_' lower(einfo_fields{plt_vars.sort_vec(1)}) 's;']);
if plt_vars.exclude_FWM
    rois(strcmp(rois,'FWM')) = [];
end
if plt_vars.exclude_OUT
    rois(strcmp(rois,'OUT')) = [];
end
roi_corr_cnt = zeros(numel(SBJs), numel(rois), numel(rois));
roi_elec_cnt = zeros(numel(SBJs), numel(rois), numel(rois));

for sbj_ix = 1:numel(SBJs)
    %% Load data
    SBJ = SBJs{sbj_ix};
    eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    
    corr_filename = [SBJ_vars.dirs.proc SBJ '_' stat_id '_' an_id_main '.mat'];
    load(corr_filename,'qvals','diff_vals','pairs');
    sig_pairs = find(qvals<=plt_vars.sig_cut);
    corr_cut_pairs = find(abs(diff_vals)>=plt_vars.corr_cut);
    good_pairs = intersect(sig_pairs,corr_cut_pairs);
    
    %% Load ROI and GM/WM info
    einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
    load(einfo_filename);
    % Select active electrodes
    load(corr_filename,'hfa_cond');
    % Clean up einfo after removing channels
    good_ix = NaN(size(hfa_cond{1}.label));
    for e_ix = 1:numel(hfa_cond{1}.label)
        good_ix(e_ix) = find(strcmp(einfo(:,1),hfa_cond{1}.label{e_ix}));
    end
    einfo = einfo(good_ix,:);
    
    einfo = sortrows(einfo,plt_vars.sort_vec);
    
    %% Compute significant percentages
    for pair_ix = 1:size(pairs,1)
        roi_ix = sort([find(strcmp(rois,einfo{pairs(pair_ix,1),plt_vars.sort_vec(1)})),...
            find(strcmp(rois,einfo{pairs(pair_ix,2),plt_vars.sort_vec(1)}))]);
        roi_elec_cnt(sbj_ix,roi_ix(1),roi_ix(2)) = roi_elec_cnt(sbj_ix,roi_ix(1),roi_ix(2))+1;
        if find(good_pairs==pair_ix)
            roi_corr_cnt(sbj_ix,roi_ix(1),roi_ix(2)) = roi_corr_cnt(sbj_ix,roi_ix(1),roi_ix(2))+1;
        end
    end
end

% Average across subjesyscts
roi_corr_per = roi_corr_cnt./roi_elec_cnt;
grp_roi_corr_per = squeeze(nanmean(roi_corr_per,1));

%% Averaged within main ROI Matrix
fprintf('=========================== GROUP Main ROI Averaged Correlation Matrix ===========================\n');
fig_name = strcat('GRP_corr_HFA_cond_acE_mat_avg_main_',plt_id);
figure('Name',fig_name,'Visible',fig_vis,'units','normalized','outerposition',plt_vars.fig_dim);
hold on;

imagesc(grp_roi_corr_per); axis xy;

ax_tick = 1:numel(rois);
set(gca,'XLim',[0.5 numel(rois)+0.5],'XTick',ax_tick,'XTickLabel',rois,'FontSize',plt_vars.font_sz);
set(gca,'YLim',[0.5 numel(rois)+0.5],'YTick',ax_tick,'YTickLabel',rois,'FontSize',plt_vars.font_sz);
xlabel(einfo_fields{plt_vars.sort_vec(1)});
ylabel(einfo_fields{plt_vars.sort_vec(1)});

subplot_name = strcat(stat_vars.event_lab,': Proportion sig by ',einfo_fields{plt_vars.sort_vec(1)});
% title(subplot_name,'interpreter','none');
cbar = colorbar;
caxis([0 max(grp_roi_corr_per(:))]);
title('% Sig. Different Correlations Across Electrodes');

% % One color bar for whole figure
% c=axes('OuterPosition', [0.5 0 0.5 1], 'Visible', 'off');
% tickDiv=(colors(2)-colors(1))/5;
% colorbar('YTick',[colors(1):tickDiv:colors(2)]);
% caxis(colors);

%         suptitle(strcat(SBJ,',',task,': ',fband_pairs{fband},' (combined) ',sortID,'_cscale.',colorscale));
% Save plot
if save_fig
    results_dir = ['~/PRJ_Stroop/results/HFA/GRP_',stat_id,'/',an_id_main,'/'];
    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end
    saveas(gcf,[results_dir,fig_name,'.',fig_filetype]);
end

% Print results
fprintf('WITHIN ROIs:\n');
for roi_ix = 1:numel(rois)
    fprintf('\t%s: %.03f (%i / %i)\n',rois{roi_ix},grp_roi_corr_per(roi_ix,roi_ix),...
        sum(roi_corr_cnt(:,roi_ix,roi_ix)),sum(roi_elec_cnt(:,roi_ix,roi_ix)));
end
fprintf('\nACROSS ROIs:\n');
for roi_ix1 = 1:numel(rois)
    for roi_ix2 = roi_ix1+1:numel(rois)
        fprintf('\t%s-%s: %.03f (%i / %i)\n',rois{roi_ix1},rois{roi_ix2},grp_roi_corr_per(roi_ix1,roi_ix2),...
            sum(roi_corr_cnt(:,roi_ix1,roi_ix2)),sum(roi_elec_cnt(:,roi_ix1,roi_ix2)));
    end
end
fprintf('\n');


%% Print results
fprintf('=========================== SBJ Main ROI Averaged Correlation Matrix ===========================\n');
% Print WITHIN ROI
fprintf('WITHIN ROIs:\n');
% SBJ Header
fprintf('\t');
for sbj_ix = 1:numel(SBJs)
    fprintf('\t%s\t',SBJs{sbj_ix});
end
fprintf('\n');
for roi_ix = 1:numel(rois)
    fprintf('\t%s:\t',rois{roi_ix});
    for sbj_ix = 1:numel(SBJs)
        fprintf('%.03f (%i/%i)\t',roi_corr_per(sbj_ix,roi_ix,roi_ix),...
            roi_corr_cnt(sbj_ix,roi_ix,roi_ix),roi_elec_cnt(sbj_ix,roi_ix,roi_ix));
    end
    fprintf('\n');
end
fprintf('\n\nACROSS ROIs:\n');
% SBJ Header
fprintf('\t\t');
for sbj_ix = 1:numel(SBJs)
    fprintf('\t%s\t',SBJs{sbj_ix});
end
fprintf('\n');
for roi_ix1 = 1:numel(rois)
    for roi_ix2 = roi_ix1+1:numel(rois)
        fprintf('\t%s-%s:\t',rois{roi_ix1},rois{roi_ix2});
        for sbj_ix = 1:numel(SBJs)
            fprintf('%.03f (%i/%i)\t',roi_corr_per(sbj_ix,roi_ix1,roi_ix2),...
                roi_corr_cnt(sbj_ix,roi_ix1,roi_ix2),roi_elec_cnt(sbj_ix,roi_ix1,roi_ix2));
        end
        fprintf('\n');
    end
end
fprintf('\n');

%% Averaged within secondary ROI Matrix
if plot_scnd_mtx
    error('2nd matrix not set for GRP yet');
    fprintf('=========================== Secondary ROI Averaged Correlation Matrix ===========================\n');
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


end
