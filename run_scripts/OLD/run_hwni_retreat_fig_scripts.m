%% Set up variables to enter function
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};

conditions  = 'CI';
pipeline_id = 'main_ft';
actv_win    = '100';
save_fig    = 1;
fig_vis     = 'on';

%% RT behavior group level
B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis);

%% HG active examples
% None! no room on this poster...

%% HG condition examples
an_id_s = 'HGm_S_zbootS_trl2to15_sm10_stat15';
an_id_r = 'HGm_R_zbootS_trl5to1_sm10_stat5to1';
plt_id  = 'ts_S15R1_errbr_evnt';
SBJ_elecs = {{'IR35','LOF6-7'},{'IR41','LIN6-7'},{'IR35','LPC5-6'},{'IR35','LAC1-2'},{'IR39','ROF7-8'},...
             {'IR41','RSM3-4'},{'IR41','RIN3-4'},{'IR35','RIN3-4'},{'IR39','LOF3-4'}};
for ix = 1:numel(SBJ_elecs)
    SBJ08b_HFA_plot_SR_stats_svg(SBJ_elecs{ix}{1},SBJ_elecs{ix}{2},conditions,pipeline_id,...
        an_id_s,an_id_r,plt_id,save_fig,fig_vis);
end

%% ERP examples
% probably not happening here...

%% Condition Onset by gROI
an_id       = 'HGm_S_zbtS_trl2to15_sm10_wn30_stat15';
plt_id      = 'onsets_trl0to15_evnt_roi';
SBJ08c_HFA_GRP_cond_ROI_onset_order_RTout(SBJs,conditions,pipeline_id,an_id,plt_id,1,fig_vis);

%% Proportions of significant effects across ROI
an_id       = 'HGm_S_zbootS_trl2to15_sm10_stat15';
plt_id      = 'onsets_trl0to15_evnt_roi';
SBJ08c_HFA_GRP_summary_bar_perc_actv_cond_ROI(SBJs,conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis);

an_id       = 'HGm_R_zbootS_trl5to1_sm10_stat5to1';
plt_id      = 'onsets_trl5to1_evnt_roi';
SBJ08c_HFA_GRP_summary_bar_perc_actv_cond_ROI(SBJs,conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis);

%% Just to get the significant elecs (can't read labels on saved versions)
an_id       = 'HGm_S_zbootS_trl2to15_sm10_stat15';
plt_id      = 'summary_trl0to15_evnt_roi';
for sbj_ix = 1:numel(SBJs)
    SBJ08c_HFA_summary_plot_actv_cond_ROI(SBJs{sbj_ix},conditions,pipeline_id,an_id,actv_win,plt_id,0,fig_vis)
end

an_id       = 'HGm_R_zbootS_trl5to1_sm10_stat5to1';
plt_id      = 'summary_trl5to1_evnt_roi';
for sbj_ix = 1:numel(SBJs)
    SBJ08c_HFA_summary_plot_actv_cond_ROI(SBJs{sbj_ix},conditions,pipeline_id,an_id,actv_win,plt_id,0,fig_vis)
end
