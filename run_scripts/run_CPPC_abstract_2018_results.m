%% Set up variables to enter function
SBJs = {'IR21','IR27','IR31','IR35','IR39','IR41','IR52','IR54','IR57','IR61','IR68','IR72','IR74'};%'IR26','IR32','IR65',
pipeline_id = 'main_ft';
actv_win    = '100';
save_fig    = 1;
fig_vis     = 'on';

%% RT behavior group level
% didn't change...
for s = 1:numel(SBJs)
    B00_RT_SBJ_hist(SBJs{s},save_fig,fig_vis,'png');
end
B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis);
B00_RT_GRP_hist(SBJs,'CNI',save_fig,fig_vis);

%% Print Results according to GM %
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
roi_id  = 'gROI';
atlas_id = 'Dx';
gm_thresh = 0.01;
save_out = 0;

% SBJ10c_HFA_GRP_elec_cnt_tis_RT_ANOVA(SBJs,stat_id,an_id_s,roi_id,atlas_id,gm_thresh,pipeline_id,save_out)
SBJ10c_HFA_GRP_elec_cnt_tis_RT_ANOVA(SBJs,stat_id,an_id_r,roi_id,atlas_id,gm_thresh,pipeline_id,save_out)

%% Proportions of significant effects across ROI
stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
atlas_id    = 'Dx';
an_list     = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
plt_id      = 'onsets_trl0to15_evnt_roi';
roi_list    = {'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
gm_lim      = [0 0.01 0.1];
fig_vis     = 'on';
fig_filetype= 'png';
roi_ix = 2;
for gm_ix = 1:numel(gm_lim)
    for an_ix = 1:numel(an_list)
        SBJ10c_HFA_GRP_summary_bar_perc_actv_GMlim_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
            an_list{an_ix},actv_win,roi_list{roi_ix},atlas_id,gm_lim(gm_ix),plt_id,save_fig,fig_vis,fig_filetype);
    end
end

%% Plot HFA Correlations
stat_id = 'corr_HFA_CI_ts_filt0to6_mS0to25_aS0to1_aR5to1_nb1k';
an_id_main = 'HGm_S_zbtS_trl2to251_sm0_wn100_stat0';
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
plt_id  = 'corr_mat_ts_S0to15_R5to10_gR_q01_r0';
plot_ts = 0;        % 0=none, 1=all, n=highest n correlations
plot_scnd_mtx = 0;
fig_filetype = 'png';
for sbj_ix = 1:numel(SBJs)
    fprintf('Running for %s\n',SBJs{sbj_ix});
%     SBJ11a_corr_HFA_acE_surr(SBJs{sbj_ix},pipeline_id,stat_id,an_id_s,plt_id)
%     SBJ11a_corr_HFA_acE_surr(SBJs{sbj_ix},pipeline_id,stat_id,an_id_r,plt_id)
    SBJ11b_corr_HFA_acE_plot_mat_ts(SBJs{sbj_ix},pipeline_id,stat_id,an_id_main,plt_id,...
        plot_ts,plot_scnd_mtx,save_fig,fig_vis,fig_filetype)
%     pause;
end

SBJ11c_corr_HFA_cond_acE_plot_grp_mat(SBJs,pipeline_id,stat_id,an_id_main,plt_id,...
                                                plot_scnd_mtx,save_fig,fig_vis,fig_filetype)
                                            
%% Plot TFRs
an_id = 'TFR_R_wvlt_f4to60_zsc3to1_trl5to101_stat5to1';
conditions = 'CI';
plt_id = 'tfr_R_BSLNabs_CBARmaxabs_tk2';
for sbj_ix = 1:numel(SBJs)
    SBJ09b_TFR_plot_stats(SBJs{sbj_ix},conditions,an_id,plt_id,save_fig,fig_vis)
end

% %% Run ANOVA
% an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% % submit: SBJ10a_corrRT_regressRT_ANOVA(SBJ,an_id,stat_id)
% 
% %% Plot ANOVA with RT correlation
% stat_id = 'corrRT_CNI_pcon_WL200_WS50';
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
% for sbj_ix = 1:numel(SBJs)
%     fprintf('Plotting for %s\n',SBJs{sbj_ix});
%     SBJ10b_ANOVA_plot_SR_RTcorr(SBJs{sbj_ix},pipeline_id,stat_id,an_id_s,an_id_r,plt_id,1,'off','svg');
%     close all;
% end
% 
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

% %% Plot onsets of ANOVA+RT
% stat_id = 'corrRT_CNI_pcon_WL200_WS50';
% roi_id  = 'gROI';
% an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% plt_id      = 'onsets_trl0to15_evnt_roi';
% fig_filetype = 'svg';
% SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on',fig_filetype)
% 
% SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,...
%                                                     plt_id,save_fig,fig_vis,fig_filetype)
%                                                 
% an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% plt_id      = 'onsets_trl5to1_evnt_roi';
% roi_id = 'thryROI';
% SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on',fig_filetype)
% 
% SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,...
%                                                     plt_id,save_fig,fig_vis,fig_filetype)
% 
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

%% HG condition examples
% % might be interesting to compare mean traces, but I didn't run it because
% % I haven't wanted to run the SBJ08a_HFA_stats for 'CI'...
% conditions = 'CI';
% an_id_s = 'HGm_S_zbtS_trl2to15_sm10_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to1_sm10_wn100_stat5to1';
% plt_id  = 'ts_S15R1_errbr_evnt';
% SBJ_elecs = {{'IR35','LOF6-7'},{'IR41','LIN6-7'},{'IR35','LPC5-6'},{'IR35','LAC1-2'},{'IR39','ROF7-8'},...
%              {'IR41','RSM3-4'},{'IR41','RIN3-4'},{'IR35','RIN3-4'},{'IR39','LOF3-4'}};
%          %IR35 RIN3-4 used for SFN poster
% for ix = 1:numel(SBJ_elecs)
%     SBJ08b_HFA_plot_SR_stats_svg(SBJ_elecs{ix}{1},SBJ_elecs{ix}{2},conditions,pipeline_id,...
%         an_id_s,an_id_r,plt_id,save_fig,fig_vis);
% end
