if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
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

%% ================================================================================
%   BEHAVIOR
%  =================================================================================
%% RT behavior group level
% didn't change...
for s = 1:numel(SBJs)
    B00_RT_SBJ_hist(SBJs{s},save_fig,fig_vis,'png');
end
B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis);
B00_RT_GRP_hist(SBJs,'CNI',save_fig,fig_vis);
% IR21 is the only one with a significant CSE behavioral effect (cI > iI)
% B00_RT_SBJ_hist_CSE(SBJs(1),'CSE',0,fig_vis);

%% ================================================================================
%   QUALITY CHECKS
%  =================================================================================
% %% Print Results according to GM %
% stat_id = 'corrRT_CNI_pcon_WL200_WS50';
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% roi_id  = 'gROI';
% atlas_id = 'Dx';
% gm_thresh = 0.01;
% save_out = 0;
% 
% % SBJ10c_HFA_GRP_elec_cnt_tis_RT_ANOVA(SBJs,stat_id,an_id_s,roi_id,atlas_id,gm_thresh,pipeline_id,save_out)
% SBJ10c_HFA_GRP_elec_cnt_tis_RT_ANOVA(SBJs,stat_id,an_id_r,roi_id,atlas_id,gm_thresh,pipeline_id,save_out)

% %% Print Results according to activation (what elecs are non-active but sig?)
% stat_id = 'corrRT_CNI_pcon_WL200_WS50';
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% roi_id  = 'ROI';
% atlas_id = 'Dx';
% gm_thresh = 0;
% 
% % SBJ10c_HFA_GRP_elec_cnt_nact_RT_ANOVA(SBJs,stat_id,an_id_s,actv_win,roi_id,atlas_id,gm_thresh,pipeline_id)
% SBJ10c_HFA_GRP_elec_cnt_nact_RT_ANOVA(SBJs,stat_id,an_id_r,actv_win,roi_id,atlas_id,gm_thresh,pipeline_id)

%% ================================================================================
%   Effects by ROI
%  =================================================================================
%% Proportions of significant effects across ROI
stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
atlas_id    = {'Yeo7','Dx'};
an_opts     = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
plt_id      = 'onsets_trl0to15_evnt_roi';
roi_opts    = {'Yeo7','mgROI'};%{'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
gm_lim      = [0 0.01 0.1];
plot_out    = 0;
fig_vis     = 'on';
fig_filetype= 'png';
% roi_ix = 1;
gm_ix = 1;
for roi_ix = 1:numel(roi_opts)
    for an_ix = 1%1:numel(an_opts)
        SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_CSE_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
            an_opts{an_ix},roi_opts{roi_ix},atlas_id{roi_ix},gm_lim(gm_ix),plt_id,plot_out,save_fig,fig_vis,fig_filetype);
%         SBJ10c_HFA_GRP_summary_bar_perc_actv_GMlim_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
%             an_opts{an_ix},actv_win,roi_opts{roi_ix},atlas_id{roi_ix},gm_lim(gm_ix),plt_id,save_fig,fig_vis,fig_filetype);
    end
end

% Individual SBJs
% roi_ix = 1;
% for an_ix = 1:numel(an_opts)
%     SBJ10c_HFA_SBJ_summary_bar_perc_actv_GMlim_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
%         an_opts{an_ix},actv_win,roi_opts{roi_ix},atlas_id,0,plt_id,save_fig,fig_vis,fig_filetype);
% end

%% ================================================================================
%   EFFECT ONSETS
%  =================================================================================
%% Plot onsets of ANOVA+RT
stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
clust_id    = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
an_opts     = {'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};%'HGm_S_zbtS_trl2to151_sm0_wn100_stat15',
gm_thresh   = 0;
median_yn   = 0;
roi_opts    = {'mgROI'};%,'Yeo7'};%{'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
atlas_opts  = {'Dx'};%,'Yeo7'};%'Dx';
% plt_opts    = {{'onsets_trl5to1_evnt_roi'}};%{'onsets_trl0to15_evnt_roi'},
plt_opts    = {{'onsets_trl0to15_violin_allROI','onsets_trl5to1_violin_avgROI'}};%{{'onsets_trl0to15_violin_allSBJ','onsets_trl0to15_violin_allROI','onsets_trl0to15_violin_avgROI'},...
%                {'onsets_trl5to1_violin_allSBJ','onsets_trl5to1_violin_allROI','onsets_trl5to1_violin_avgROI'}};
fig_filetype = 'png';%'svg'
for roi_ix = 1%:numel(roi_opts)
    for an_ix = 1:numel(an_opts)
        for plt_ix = 1%:numel(plt_opts{1})
            fprintf('roi: %s; an: %s; plt: %s\n',roi_opts{roi_ix},an_opts{an_ix},plt_opts{an_ix}{plt_ix});
%             % Violin Plots
            SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
                                                    atlas_opts{roi_ix},gm_thresh,plt_opts{an_ix}{plt_ix},1,'on',fig_filetype)
            SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA_clustBin(SBJs,clust_id,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
                                                    atlas_opts{roi_ix},gm_thresh,plt_opts{an_ix}{plt_ix},1,'on',fig_filetype)
            % Music plots- Normalized
%             SBJ10c_HFA_GRP_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},median_yn,roi_opts{roi_ix},...
%                                                     atlas_opts{roi_ix},0,plt_opts{an_ix}{plt_ix},1,'off',fig_filetype);
%             close all;
%             % Music plots- subject-specific RTs
%         SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
%                                                     atlas_id,0,plt_opts{an_ix},1,'on',fig_filetype)
%             % Pair differences
%         SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
%                                                     atlas_id,0,plt_opts{an_ix},save_fig,fig_vis)
        end
    end
end

%% ================================================================================
%   EFFECT COMMONALITIES ACROSS ROIS
%  =================================================================================
% %% Plot ANOVA with RT correlation by ROI to see any patterns
% stat_id = 'corrRT_CNI_pcon_WL200_WS50';
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
% roi_opts    = {'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
% 
% roi_ix = 2;
% for sbj_ix = 1:numel(SBJs)
%     fprintf('Plotting for %s\n',SBJs{sbj_ix});
%     SBJ10b_ANOVA_plot_SR_RTcorr_ROIcomb(SBJs{sbj_ix},pipeline_id,stat_id,an_id_s,an_id_r,atlas_id,...
%         roi_opts{roi_ix},plt_id,1,'on','svg');
%     close all;
% end
% 
% %% Plot HFA Correlations
% stat_id    = 'corr_HFA_ts_filt0to6_mS0to25_aS0to1_aR5to1_nb1k';
% an_id_main = 'HGm_S_zbtS_trl2to251_sm0_wn100_stat0';
% roi_id     = 'Yeo7';%'gROI','thryROI','Yeo7'
% atlas_id   = 'Yeo7';%'Dx','Yeo7'
% plt_id     = 'corr_mat_ts_S0to15_R5to10_gR_q01_r1';
% plot_hist  = 0;
% plot_ts    = 0;        % 0=none, 1=all, n=highest n correlations
% plot_scnd_mtx = 0;
% fig_filetype = 'png';
% for sbj_ix = 1:numel(SBJs)
%     fprintf('Running for %s\n',SBJs{sbj_ix});
% %     SBJ11a_corr_HFA_acE_surr(SBJs{sbj_ix},pipeline_id,stat_id,an_id_main,atlas_id,roi_id,plt_id)
%     SBJ11b_corr_HFA_acE_plot_mat_ts(SBJs{sbj_ix},pipeline_id,stat_id,an_id_main,roi_id,atlas_id,plt_id,...
%         plot_ts,plot_scnd_mtx,plot_hist,save_fig,'off',fig_filetype)
% %     pause;
% end
% 
% SBJ11c_corr_HFA_acE_plot_mat_grp(SBJs,pipeline_id,stat_id,an_id_main,roi_id,atlas_id,plt_id,...
%                                                 plot_scnd_mtx,save_fig,fig_vis,fig_filetype);

%% Cluster ANOVA time series
clust_id = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
an_id = an_opts{2};
roi_id   = 'mgROI';
atlas_id = 'Dx';
plt_id   = 'ts_S0to15_R5to10_evnt_sigline';
plot_out = 0;
fig_vis  = 'off';
save_fig = 1;
for sbj_ix = 1:numel(SBJs)
    for an_ix = 1:2%1:numel(clust_opts)
        SBJ12a_corr_ANOVA_ts_cluster(SBJs{sbj_ix},pipeline_id,stat_id,clust_id,an_opts{an_ix},...
            roi_id,atlas_id,plt_id,fig_vis,save_fig,plot_out);
%         fn_view_recon_clust(SBJ, pipeline_id, clust_id, stat_id, an_id, ...
%                     view_space, reg_type, show_labels, hemi, plot_out)
        close all;
    end
end

%% Plot clusters on the brain
clust_id = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
hemi_opts = {'r','l'};
roi_id   = 'mgROI';
atlas_id = 'Dx';
% plt_id   = 'ts_S0to15_R5to10_evnt_sigline';
plot_ns  = 1;
plot_out = 0;
plot_clusters = 1;
plot_recon = 0;
show_labels = 1;
fig_vis  = 'off';
save_fig = 1;
for sbj_ix = 1:numel(SBJs)
    for an_ix = 2%1:2
        for sig_ix = 0:1%:2
            fn_view_recon_clust(SBJs{sbj_ix}, clust_id, stat_id, an_opts{an_ix}, atlas_id, roi_id,...
                    'pat', '', show_labels, 'b', plot_out, sig_ix, plot_clusters, plot_recon, fig_vis, save_fig);%hemi_opts{hemi_ix}
        end
        close all;
    end
end

% %% Plot TFRs
% an_id = 'TFR_R_wvlt_f4to60_zsc3to1_trl5to101_stat5to1';
% conditions = 'CI';
% plt_id = 'tfr_R_BSLNabs_CBARmaxabs_tk2';
% for sbj_ix = 1:numel(SBJs)
%     SBJ09b_TFR_plot_stats(SBJs{sbj_ix},conditions,an_id,plt_id,save_fig,fig_vis)
% end
% 
% % %% Run ANOVA
% an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% % an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% % % submit: SBJ10a_corrRT_regressRT_ANOVA(SBJ,an_id,stat_id)
% % 
%% # elecs and # sig per ROI per SBJ
% WARNING: This reveals a non-trivial amount of change after last ANOVA
% window added! No idea why the RT analysis results are changing t hough!
% not really using these in the poster
% stat_id = 'corrRT_CNI_pcon_WL200_WS50';
% roi_id = 'gROI';
% an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)
% 
% an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)

% %% Plot time series of significant effects by ROI
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% roi_id  = 'gROI';
% fig_vis = 'off';
% plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
% SBJ10b_HFA_plot_SR_ROI_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
%                                                 roi_id,plt_id,save_fig,fig_vis,'png')
% 
% plt_id  = 'ts_S0to15_R5to10_evnt_siglineSM';
% SBJ10c_HFA_GRP_plot_SR_ROIavg_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
%                                                 roi_id,plt_id,save_fig,fig_vis,'png')

%% HG active examples
conditions  = 'CNI';
pipeline_id = 'main_ft';
actv_win    = '100';
an_id_s     = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
an_id_r     = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
plt_id      = 'stack_S2to15_R5to10_evnt_c5';
save_fig    = 1;
fig_vis     = 'off';
fig_filetype= 'png';
for sbj_ix = 1:numel(SBJs)
    fprintf('Plotting for %s\n',SBJs{sbj_ix});
    SBJ08b_HFA_plot_SR_stack_cond_onset_noLab(SBJs{sbj_ix},conditions,an_id_s,an_id_r,...
        pipeline_id,actv_win,plt_id,save_fig,fig_vis,fig_filetype)
    close all;
end

% %% Link all sig results
% SBJ10b_HFA_link_sig_ANOVA(SBJs,stat_id,an_id_s,an_id_r)

%% CSE Effects
% % might be interesting to compare mean traces, but I didn't run it because
% % I haven't wanted to run the SBJ08a_HFA_stats for 'CI'...
conditions  = 'CSE';
an_opts     = {'HGm_S_zbtS_trl2to151_sm10_wn100_stat15','HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1'};
plt_id      = 'ts_S15R1_errbr_evnt';
save_fig    = 1;
fig_vis     = 'off';
view_space  = 'pat';
reg_type    = '';
show_labels = 1;
hemi        = 'b';
for ix = 1:numel(SBJs)
    for an_ix = 1:numel(an_opts)
        fn_view_recon_stat(SBJs{ix}, pipeline_id, conditions, an_opts{an_ix}, view_space, reg_type, show_labels, hemi);
    end
    pause
%     SBJ08b_HFA_plot_SR_stats(SBJs{ix},conditions,an_opts{1},an_opts{2},plt_id,save_fig,fig_vis);
%     close all;
end

% SBJ_elecs = {{'IR35','LOF6-7'},{'IR41','LIN6-7'},{'IR35','LPC5-6'},{'IR35','LAC1-2'},{'IR39','ROF7-8'},...
%              {'IR41','RSM3-4'},{'IR41','RIN3-4'},{'IR35','RIN3-4'},{'IR39','LOF3-4'}};
%          %IR35 RIN3-4 used for SFN poster
% for ix = 1:numel(SBJ_elecs)
%     SBJ08b_HFA_plot_SR_stats_svg(SBJ_elecs{ix}{1},SBJ_elecs{ix}{2},conditions,pipeline_id,...
%         an_id_s,an_id_r,plt_id,save_fig,fig_vis);
% end

%% Plot recons and save nice figs
% an_opts     = {'HGm_S_zbtS_trl2to151_sm10_wn100_stat15','HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1'};
an_opts     = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};

an_ix = 1;
fn_view_recon_atlas_grp_stat_sphr(SBJs,pipeline_id,'corrRT_CNI_pcon_WL200_WS50',an_opts{an_ix},'v',0,'r','Yeo7','Yeo7',1);
fn_view_recon_atlas_grp_stat_sphr(SBJs,pipeline_id,'corrRT_CNI_pcon_WL200_WS50',an_opts{an_ix},'v',0,'l','Yeo7','Yeo7',1);

