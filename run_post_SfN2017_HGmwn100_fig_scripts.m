%% Set up variables to enter function
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};

stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
% conditions  = 'CI';
pipeline_id = 'main_ft';
actv_win    = '100';
save_fig    = 1;
fig_vis     = 'off';

%% RT behavior group level
% didn't change...
% B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis);

%% Run ANOVA
an_id       = 'HGm_S_zbtS_trl2to151_sm10_wn100_stat15';
an_id       = 'HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1';
% submit: SBJ10a_corrRT_regressRT_ANOVA(SBJ,an_id,stat_id)

%% Plot ANOVA with RT correlation
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
an_id_s = 'HGm_S_zbtS_trl2to151_sm10_wn100_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1';
plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
for sbj_ix = 1:numel(SBJs)
    fprintf('Plotting for %s\n',SBJs{sbj_ix});
    SBJ10b_ANOVA_plot_SR_RTcorr(SBJs{sbj_ix},stat_id,an_id_s,an_id_r,plt_id,1,'off','svg');
    close all;
end

%% Proportions of significant effects across ROI
an_list     = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
plt_id      = 'onsets_trl0to15_evnt_roi';
roi_list    = {'gROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
fig_vis     = 'on';
fig_filetype= 'png';
for roi_ix = 1:numel(roi_list)
    for an_ix = 1:numel(an_list)
        SBJ10c_HFA_GRP_summary_bar_perc_actv_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
            an_list{an_ix},actv_win,roi_list{roi_ix},plt_id,save_fig,fig_vis,fig_filetype);
    end
end

%% Overlaps between analysis IDs
% Compare HGm smoohting 500ms vs. 100ms
an_id1       = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
an_id2       = 'HGm_S_zbtS_trl2to151_sm10_wn100_stat15';
SBJ10c_HFA_GRP_boxcmp_RT_ANOVA_an_ids(SBJs,stat_id,an_id1,an_id2,1)

an_id1       = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
an_id2       = 'HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1';
SBJ10c_HFA_GRP_boxcmp_RT_ANOVA_an_ids(SBJs,stat_id,an_id1,an_id2,1)
        
% Compare smoothing vs. no smoothing
an_id1       = 'HGm_S_zbtS_trl2to151_sm10_wn100_stat15';
an_id2       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
SBJ10c_HFA_GRP_boxcmp_RT_ANOVA_an_ids(SBJs,stat_id,an_id1,an_id2,1)

an_id1       = 'HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1';
an_id2       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
SBJ10c_HFA_GRP_boxcmp_RT_ANOVA_an_ids(SBJs,stat_id,an_id1,an_id2,1)

%% # elecs and # sig per ROI per SBJ
% WARNING: This reveals a non-trivial amount of change after last ANOVA
% window added! No idea why the RT analysis results are changing t hough!
% not really using these in the poster
roi_id = 'gROI';
an_id       = 'HGm_S_zbtS_trl2to151_sm10_wn100_stat15';
SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)

an_id       = 'HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1';
SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)

%% Plot onsets of ANOVA+RT
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
roi_id  = 'gROI';
an_id       = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
plt_id      = 'onsets_trl0to15_evnt_roi';
fig_filetype = 'png';
SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on',fig_filetype)

an_id       = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on',fig_filetype)

%% Plot time series of significant effects by ROI
an_id_s = 'HGm_S_zbtS_trl2to151_sm10_wn100_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1';
roi_id  = 'gROI';
fig_vis = 'off';
plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
SBJ10b_HFA_plot_SR_ROI_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
                                                roi_id,plt_id,save_fig,fig_vis,'png')

plt_id  = 'ts_S0to15_R5to10_evnt_siglineSM';
SBJ10c_HFA_GRP_plot_SR_ROIavg_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
                                                roi_id,plt_id,save_fig,fig_vis,'png')

%% HG active examples
conditions  = 'CNI';
pipeline_id = 'main_ft';
actv_win    = '100';
an_id_s     = 'HGm_S_zbtS_trl2to151_sm10_wn100_stat15';
an_id_r     = 'HGm_R_zbtS_trl5to101_sm10_wn100_stat5to1';
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
