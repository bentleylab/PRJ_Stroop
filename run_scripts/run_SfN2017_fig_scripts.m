%% Set up variables to enter function
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};

stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
% conditions  = 'CI';
pipeline_id = 'main_ft';
actv_win    = '100';
save_fig    = 1;
fig_vis     = 'off';

%% RT behavior group level
B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis);

%% Run ANOVA
an_id       = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
an_id       = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
% submit: SBJ10a_corrRT_regressRT_ANOVA(SBJ,an_id,stat_id)

%% Plot ANOVA with RT correlation
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
an_id_s = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
for sbj_ix = 1:numel(SBJs)
    fprintf('Plotting for %s\n',SBJs{sbj_ix});
    SBJ10b_ANOVA_plot_SR_RTcorr(SBJs{sbj_ix},stat_id,an_id_s,an_id_r,plt_id,1,'off','svg');
    close all;
end

%% Proportions of significant effects across ROI
an_list     = {'HGm_S_zbtS_trl2to151_sm10_wn30_stat15','HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1'};
plt_id      = 'onsets_trl0to15_evnt_roi';
roi_list    = {'gROI'};%,'LPFC','MPFC','INS','OFC'};
fig_vis     = 'on';
for roi_ix = 1:numel(roi_list)
    for an_ix = 1:numel(an_list)
        SBJ10c_HFA_GRP_summary_bar_perc_actv_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
            an_list{an_ix},actv_win,roi_list{roi_ix},plt_id,save_fig,fig_vis,'svg');
    end
end

%% Overlaps between old and new analyses
% noto really using these in the poster
ci_an_id       = 'HGm_S_zbtS_trl2to15_sm10_wn30_stat15';
an_id       = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
SBJ10c_HFA_GRP_boxcmp_RT_vs_CNI_ANOVA(SBJs,stat_id,ci_an_id,an_id,1)
% SBJ10c_HFA_GRP_boxcmp_CI_vs_RT_ANOVA(SBJs,stat_id,an_id,1)

an_id       = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
SBJ10c_HFA_GRP_boxcmp_CI_vs_RT_ANOVA(SBJs,stat_id,an_id,1)

%% # elecs and # sig per ROI per SBJ
% WARNING: This reveals a non-trivial amount of change after last ANOVA
% window added! No idea why the RT analysis results are changing t hough!
% not really using these in the poster
roi_id = 'gROI';
an_id       = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)

an_id       = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
SBJ10c_HFA_GRP_elec_cnt_ROI_RT_ANOVA(SBJs,stat_id,an_id,roi_id,pipeline_id,1)

%% Plot onsets of ANOVA+RT
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
roi_id  = 'gROI';
an_id       = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
plt_id      = 'onsets_trl0to15_evnt_roi';
SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on','svg')

% !!! Something wrong with RT markers in R_locked
%   probably not going to use these anyways...
% an_id       = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
% plt_id      = 'onsets_trl5to1_evnt_roi';
% SBJ10c_HFA_GRP_onsets_ROI_RTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_id,roi_id,plt_id,1,'on','svg')

%% Plot time series of significant effects by ROI
an_id_s = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
roi_id  = 'gROI';
fig_vis = 'on';
plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
SBJ10b_HFA_plot_SR_ROI_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
                                                roi_id,plt_id,save_fig,fig_vis,'svg')

plt_id  = 'ts_S0to15_R5to10_evnt_siglineSM';
SBJ10c_HFA_GRP_plot_SR_ROIavg_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
                                                roi_id,plt_id,save_fig,fig_vis,'svg')
% Showing with error bars or butterfly style looks really shitty
                                            % plt_avg = 1;
% SBJ10c_HFA_GRP_plot_SR_ROIavg_but_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
%                                                 roi_id,plt_id,plot_avg,save_fig,fig_vis,'svg')

% Too much mess to see anything:
% plt_id  = 'ts_S0to15_R5to10_evnt_siglineSM';
% SBJ10c_HFA_GRP_plot_SR_ROI_RTcorr_ANOVA(SBJs,stat_id,pipeline_id,an_id_s,an_id_r,...
%                                                 roi_id,plt_id,save_fig,fig_vis)

%% HG active examples
SBJ = 'IR35';
conditions  = 'CNI';
pipeline_id = 'main_ft';
actv_win    = '100';
an_id_s     = 'HGm_S_zbtS_trl2to151_sm10_wn30_stat15';
an_id_r     = 'HGm_R_zbtS_trl5to101_sm10_wn30_stat5to1';
plt_id      = 'stack_S2to15_R5to10_evnt_c5';
save_fig    = 1;
fig_vis     = 'on';
fig_filetype= 'svg';
SBJ08b_HFA_plot_SR_stack_cond(SBJ,conditions,an_id_s,an_id_r,actv_win,plt_id,save_fig,fig_vis,fig_filetype)


%% HG condition examples
conditions = 'CI';
an_id_s = 'HGm_S_zbtS_trl2to15_sm10_wn30_stat15';
an_id_r = 'HGm_R_zbtS_trl5to1_sm10_wn30_stat5to1';
plt_id  = 'ts_S15R1_errbr_evnt';
% SBJ_elecs = {{'IR35','LOF6-7'},{'IR41','LIN6-7'},{'IR35','LPC5-6'},{'IR35','LAC1-2'},{'IR39','ROF7-8'},...
%              {'IR41','RSM3-4'},{'IR41','RIN3-4'},{'IR35','RIN3-4'},{'IR39','LOF3-4'}};
SBJ08b_HFA_plot_SR_stats_svg('IR35','RIN3-4',conditions,pipeline_id,...
    an_id_s,an_id_r,plt_id,save_fig,fig_vis);

%% ERP examples
% probably not happening here...

% %% Just to get the significant elecs (can't read labels on saved versions)
% an_id       = 'HGm_S_zbootS_trl2to15_sm10_stat15';
% plt_id      = 'summary_trl0to15_evnt_roi';
% for sbj_ix = 1:numel(SBJs)
%     SBJ08c_HFA_summary_plot_actv_cond_ROI(SBJs{sbj_ix},conditions,pipeline_id,an_id,actv_win,plt_id,0,fig_vis)
% end
% 
% an_id       = 'HGm_R_zbootS_trl5to1_sm10_stat5to1';
% plt_id      = 'summary_trl5to1_evnt_roi';
% for sbj_ix = 1:numel(SBJs)
%     SBJ08c_HFA_summary_plot_actv_cond_ROI(SBJs{sbj_ix},conditions,pipeline_id,an_id,actv_win,plt_id,0,fig_vis)
% end
