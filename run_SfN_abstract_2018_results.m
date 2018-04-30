%% Set up variables to enter function
SBJs = {'IR21','IR27','IR31','IR32','IR35','IR39','IR41','IR52','IR54','IR57','IR61'};

stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
% conditions  = 'CI';
pipeline_id = 'main_ft';
actv_win    = '100';
save_fig    = 1;
fig_vis     = 'on';

%% RT behavior group level
% didn't change...
B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis);
B00_RT_GRP_hist(SBJs,'CNI',save_fig,fig_vis);

%% Plot HFA Correlations
stat_id = 'corr_HFA_ts_filt0to6_S0to1_R5to1_nb1k';
an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
plt_id  = 'corr_mat_ts_S0to15_R5to10_gR';
fig_filetype = 'png';
for sbj_ix = 1:numel(SBJs)
    fprintf('Running for %s\n',SBJs{sbj_ix});
%     SBJ11a_corr_HFA_acE_surr(SBJs{sbj_ix},pipeline_id,stat_id,an_id_s,plt_id)
%     SBJ11a_corr_HFA_acE_surr(SBJs{sbj_ix},pipeline_id,stat_id,an_id_r,plt_id)
    SBJ11b_corr_HFA_acE_plot_SR_mat_ts(SBJs{sbj_ix},pipeline_id,stat_id,an_id_s,an_id_r,plt_id,save_fig,fig_vis,fig_filetype)
end

%% Run ANOVA
an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% submit: SBJ10a_corrRT_regressRT_ANOVA(SBJ,an_id,stat_id)

%% Plot ANOVA with RT correlation
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
for sbj_ix = 1:numel(SBJs)
    fprintf('Plotting for %s\n',SBJs{sbj_ix});
    SBJ10b_ANOVA_plot_SR_RTcorr(SBJs{sbj_ix},pipeline_id,stat_id,an_id_s,an_id_r,plt_id,1,'off','svg');
    close all;
end

%% Proportions of significant effects across ROI
an_list     = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
plt_id      = 'onsets_trl0to15_evnt_roi';
roi_list    = {'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
fig_vis     = 'on';
fig_filetype= 'png';
for roi_ix = 1:numel(roi_list)
    for an_ix = 1:numel(an_list)
        SBJ10c_HFA_GRP_summary_bar_perc_actv_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
            an_list{an_ix},actv_win,roi_list{roi_ix},plt_id,save_fig,fig_vis,fig_filetype);
    end
end

%% # elecs and # sig per ROI per SBJ
% WARNING: This reveals a non-trivial amount of change after last ANOVA
% window added! No idea why the RT analysis results are changing t hough!
% not really using these in the poster
roi_id = 'gROI';
an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)

an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)

%% Plot onsets of ANOVA+RT
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
roi_id  = 'gROI';
an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
plt_id      = 'onsets_trl0to15_evnt_roi';
fig_filetype = 'svg';
SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on',fig_filetype)

SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,...
                                                    plt_id,save_fig,fig_vis,fig_filetype)
                                                
an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
plt_id      = 'onsets_trl5to1_evnt_roi';
roi_id = 'thryROI';
SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on',fig_filetype)

SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,...
                                                    plt_id,save_fig,fig_vis,fig_filetype)

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
    SBJ08b_HFA_plot_SR_stack_cond(SBJs{sbj_ix},conditions,an_id_s,an_id_r,...
        actv_win,plt_id,save_fig,fig_vis,fig_filetype)
    close all;
end

%% Link all sig results
SBJ10b_HFA_link_sig_ANOVA(SBJs,stat_id,an_id_s,an_id_r)

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
