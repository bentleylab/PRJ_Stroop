if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
%%
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Set up variables to enter function
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR65','IR68','IR72','IR74'};%ALMOST: 'CP24','IR26',  %NEVER: 'IR27','IR37','IR48',
pipeline_id = 'main_ft';
actv_win    = '100';
save_fig    = 1;
fig_vis     = 'on';
fig_type    = 'svg';

%% Plot onsets of ANOVA+RT
pipeline_id = 'main_ft';
stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
clust_id    = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
tbin_id     = 'eqROI';
an_opts     = {'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};%'HGm_S_zbtS_trl2to151_sm0_wn100_stat15',
gm_thresh   = 0;
median_yn   = 0;
roi_opts    = {'mgROI'};%,'Yeo7'};%{'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
atlas_opts  = {'Dx'};%,'Yeo7'};%'Dx';
% plt_opts    = {{'onsets_trl5to1_evnt_roi'}};%{'onsets_trl0to15_evnt_roi'},
plt_opts    = {{'onsets_trl0to15_violin_allROI','onsets_trl5to1_violin_all'}};%{{'onsets_trl0to15_violin_allSBJ','onsets_trl0to15_violin_allROI','onsets_trl0to15_violin_avgROI'},...
%                {'onsets_trl5to1_violin_allSBJ','onsets_trl5to1_violin_allROI','onsets_trl5to1_violin_avgROI'}};
fig_filetype = 'png';%'svg'
an_ix = 1;
for roi_ix = 1%:numel(roi_opts)
    for plt_ix = 1:numel(plt_opts)
%         for plt_ix = 1:numel(plt_opts{1})
            fprintf('roi: %s; an: %s; plt: %s\n',roi_opts{roi_ix},an_opts{an_ix},plt_opts{plt_ix}{plt_ix});
%             % Violin Plots
            SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA_timeBin(SBJs,tbin_id,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
                                                    atlas_opts{roi_ix},gm_thresh,plt_opts{plt_ix}{plt_ix},1,'on',fig_filetype)
%             SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
%                                                     atlas_opts{roi_ix},gm_thresh,plt_opts{an_ix}{plt_ix},1,'on',fig_filetype)
%             SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA_clustBin(SBJs,clust_id,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
%                                                     atlas_opts{roi_ix},gm_thresh,plt_opts{an_ix}{plt_ix},1,'on',fig_filetype)
%             % Pair differences
%         SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
%                                                     atlas_id,0,plt_opts{1}{2},save_fig,fig_vis)
%         end
    end
end

%% ================================================================================
%   RECONS with EFFECTS
%  =================================================================================
%% Plot group recon with mgROI
fn_view_recon_atlas_grp_sphr(SBJs,pipeline_id,'v',0,'l','Dx','mgROI',0);
fn_view_recon_atlas_grp_sphr(SBJs,pipeline_id,'v',0,'r','Dx','mgROI',0);
fn_view_recon_atlas_grp_sphr(SBJs,pipeline_id,'v',0,'l','Dx','MPFC',0);
fn_view_recon_atlas_grp_sphr(SBJs,pipeline_id,'v',0,'r','Dx','MPFC',0);
% fn_view_recon_atlas_grp_sphr_ROI(SBJs, pipeline_id, 'v', 0, 'l', 'Dx', 'MPFC', '', 0);

%% Plot CNI sig elecs
stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
hemi_opts = {'r','l'};
roi_id   = 'mgROI';
atlas_id = 'Dx';
plot_out = 0;
show_labels = 0;

an_ix = 2;
fn_view_recon_atlas_grp_stat_sphr(SBJs,pipeline_id,stat_id,an_opts{an_ix},'v',0,'r',atlas_id,roi_id,plot_out);
fn_view_recon_atlas_grp_stat_sphr(SBJs,pipeline_id,stat_id,an_opts{an_ix},'v',0,'l',atlas_id,roi_id,plot_out);

%% Plot binned clusters with sig CNI
clust_id = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
tbin_id     = 'eqROI';
stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
hemi_opts = {'r','l'};
roi_id   = 'mgROI';
atlas_id = 'Dx';
% plt_id   = 'ts_S0to15_R5to10_evnt_sigline';
plot_ns  = 0;
plot_out = 0;
plot_clusters = 0;
plot_recon = 1;
reg_type = 'v';
show_labels = 0;
fig_vis  = 'on';
save_fig = 0;
an_ix = 2;

% fn_view_recon_atlas_grp_stat_sphr_clust(SBJs, pipeline_id, clust_id, stat_id, an_opts{an_ix}, ...
%                         reg_type, show_labels, 'r', atlas_id, roi_id, plot_out)
% fn_view_recon_atlas_grp_stat_sphr_clust(SBJs, pipeline_id, clust_id, stat_id, an_opts{an_ix}, ...
%                         reg_type, show_labels, 'l', atlas_id, roi_id, plot_out)
fn_view_recon_atlas_grp_stat_sphr_tbin(SBJs, pipeline_id, tbin_id, stat_id, an_opts{an_ix}, ...
                        reg_type, show_labels, 'r', atlas_id, roi_id, plot_out)
% for sbj_ix = 1:numel(SBJs)
%     fn_view_recon_clust(SBJs{sbj_ix}, clust_id, stat_id, an_opts{an_ix}, atlas_id, roi_id,...
%         'pat', '', show_labels, 'b', plot_out, plot_ns, plot_clusters, plot_recon, fig_vis, save_fig);%hemi_opts{hemi_ix}
%     close all;
% end

%% TFRs
conditions  = 'CI';
an_id       = 'TFR_S_wvlt_f3to40_zsc3to1_trl3to15_stat15';
plt_id      = 'tfr_BSLNabs_CBARmaxabs_tk2';
save_fig    = 1;
fig_vis     = 'off';
for s = 1:numel(SBJs)
    SBJ09b_TFR_plot_stats(SBJs{s},conditions,pipeline_id,an_id,plt_id,save_fig,fig_vis);
    close all;
end

%% ERPs
conditions  = 'CI';
an_id_s     = 'ERP_S_trl1_stat1';
an_id_r     = 'ERP_R_trl5to1_stat5to1';
plt_id      = 'ts_S1R1_errbr_evnt';
save_fig    = 1;
fig_vis     = 'off';
for s = 1:numel(SBJs)
    SBJ07b_ERP_plot_SR_stats_stack(SBJs{s},conditions,an_id_s,an_id_r,plt_id,save_fig,fig_vis)
    close all;
end

% %% ================================================================================
% %   BEHAVIOR
% %  =================================================================================
% %% RT behavior group level
% B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis,fig_type);
% 
% %% ================================================================================
% %   Effects by ROI
% %  =================================================================================
% %% HG active examples
% SBJ         = 'IR74';
% conditions  = 'CNI';
% pipeline_id = 'main_ft';
% actv_win    = '100';
% an_id_s     = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r     = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% plt_id      = 'stack_S2to15_R5to10_evnt_c5';
% save_fig    = 1;
% fig_vis     = 'on';
% fig_filetype= 'svg';
% SBJ08b_HFA_plot_SR_stack_cond_onset_noLab(SBJ,conditions,an_id_s,an_id_r,...
%     pipeline_id,actv_win,plt_id,save_fig,fig_vis,fig_filetype)
% 
% % HG condition examples
% % might be interesting to compare mean traces, but I didn't run it because
% % I haven't wanted to run the SBJ08a_HFA_stats for 'CI'...
% conditions = 'CI';
% plt_id     = 'ts_S15R1_errbr_evnt';
% SBJ08b_HFA_plot_SR_stats_svg('IR57','LAC5-6',conditions,pipeline_id,...
%     an_id_s,an_id_r,plt_id,save_fig,fig_vis);
% % SBJ08b_HFA_plot_SR_stats_svg('IR72','LPC2-3',conditions,pipeline_id,...
% %     an_id_s,an_id_r,plt_id,save_fig,fig_vis);
% % SBJ08b_HFA_plot_SR_stats_svg('IR74','LIN1-2',conditions,pipeline_id,...
% %     an_id_s,an_id_r,plt_id,save_fig,fig_vis);
% 
% %% SU examples
% SBJ = 'IR82';% also ran IR75
% pipeline_id = 'SU_nlx';
% an_id = 'PSTH_R_trl5to1_bn20_sm10';
% plt_id = 'ts_trl5to1_errbr_evnt';
% plot_ISI = 0;
% fig_vis = 'off';
% save_plots = 1;
% close_plots = 0;
% fig_ftype= 'svg';
% SU02_PSTH_ft(SBJ,conditions,pipeline_id,an_id,plt_id,plot_ISI,fig_vis,save_plots,fig_ftype,close_plots)
% close all;
% 
% %% Proportions of significant effects across ROI
% stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
% atlas_id    = 'Dx';
% an_opts     = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
% plt_id      = 'onsets_trl0to15_evnt_roi';
% roi_id      = 'mgROI';%{'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
% gm_thresh   = 0;
% plot_out    = 0;
% fig_vis     = 'on';
% fig_filetype= 'svg';
% for an_ix = 1:2
%     SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_CSE_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
%         an_opts{an_ix},roi_id,atlas_id,gm_thresh,plt_id,plot_out,save_fig,fig_vis,fig_filetype);
% end

% %% ================================================================================
% %   EFFECT ONSETS
% %  =================================================================================
% %% Cluster ANOVA time series
% pipeline_id = 'main_ft';
% clust_id = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
% stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
% an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
% an_id = an_opts{2};
% roi_id   = 'mgROI';
% atlas_id = 'Dx';
% plt_id   = 'ts_S0to15_R5to10_evnt_sigline';
% plot_out = 0;
% fig_vis  = 'off';
% save_fig = 1;
% for sbj_ix = 1:numel(SBJs)
%     for an_ix = 1:2%1:numel(clust_opts)
%         SBJ12a_corr_ANOVA_ts_cluster(SBJs{sbj_ix},pipeline_id,stat_id,clust_id,an_opts{an_ix},...
%             roi_id,atlas_id,plt_id,fig_vis,save_fig,plot_out);
% %         fn_view_recon_clust(SBJ, pipeline_id, clust_id, stat_id, an_id, ...
% %                     view_space, reg_type, show_labels, hemi, plot_out)
%         close all;
%     end
% end
% 
% %% Plot binned cluster ANOVA ts
% SBJ = 'IR74';
% clust_id = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
% stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
% an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
% hemi_opts = {'r','l'};
% roi_id   = 'mgROI';
% atlas_id = 'Dx';
% plt_id   = 'ts_S0to15_R5to10_evnt_sigline';
% plot_ns  = 1;
% plot_out = 0;
% plot_clusters = 1;
% plot_recon = 0;
% show_labels = 1;
% fig_vis  = 'on';
% save_fig = 1;
% an_ix = 2
% for sig_ix = 0:1%:2
%     fn_view_recon_clust(SBJ, clust_id, stat_id, an_opts{an_ix}, atlas_id, roi_id,...
%         'pat', '', show_labels, 'b', plot_out, sig_ix, plot_clusters, plot_recon, fig_vis, save_fig);%hemi_opts{hemi_ix}
% end
% 
% % Plot clust centroid within SBJ
% SBJ12a_corr_ANOVA_ts_cluster_plotOnly(SBJ,pipeline_id,stat_id,clust_id,an_opts{an_ix},...
%                                         roi_id,atlas_id,plt_id,fig_vis,save_fig,plot_out)
% 
% % Plot clust centroids across SBJs
% plt_id   = 'ts_S0to15_R5to10_evnt_sigline';
% SBJ12b_clust_ANOVA_ts_GRP_peaks(SBJs,pipeline_id,stat_id,clust_id,an_opts{an_ix},...
%     roi_id,atlas_id,plt_id,fig_vis,save_fig,plot_out)
% 
