if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
%%
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Set up variables to enter function
SBJs = {'CP24','CP26','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};
pipeline_id = 'main_ft';
actv_win    = '100';
save_fig    = 1;
fig_vis     = 'on';
fig_type    = 'svg';

%% New tests
conditions  = 'CNI';%'CSE';%
pipeline_id = 'main_ft';
an_id_s     = 'HGh_S_zbtS_trl2to151_fLog_sm0_stat15';
an_id_r     = 'HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1';
plt_id      = 'stack_S2to15_R5to10_evnt_c5';%'ts_S15R1_errbr_evnt';
save_fig    = 1;
fig_vis     = 'off';
fig_ftype   = 'png';
for s = 15%:numel(SBJs)
    SBJ08b_HFA_plot_SR_stack_cond_saved(SBJs{s},conditions,an_id_s,an_id_r,...
                                        plt_id,save_fig,fig_vis,fig_ftype);
%     SBJ08b_HFA_plot_SR_stats(SBJs{s},conditions,an_id_s,an_id_r,plt_id,save_fig,fig_vis,fig_ftype)
    close all;
end

%% ================================================================================
%   BEHAVIOR
%  =================================================================================
%% RT behavior group level
B00_RT_GRP_hist_norm(SBJs,'CNI',save_fig,fig_vis,fig_type);

%% ================================================================================
%   RECONS with EFFECTS
%  =================================================================================
%% Plot group recon with mgROI
% fn_view_recon_atlas_grp(SBJs,pipeline_id,'v',0,'l','Dx','mgROI',0);
% fn_view_recon_atlas_grp(SBJs,pipeline_id,'v',0,'r','Dx','mgROI',0);

atlas_id = 'Dx';
roi_id = 'gROI';
reg_type = 'v';
show_lab = 0;
roi_opts = {{'l','lat'},{'r','lat'},{'l','MPFC'},{'r','MPFC'},{'l','deep'},{'r','deep'},{'b','OFC'}};

for roi_ix = [5 6]%1:numel(roi_opts)
    fn_view_recon_atlas_grp_ROI(SBJs, pipeline_id, reg_type, show_lab,...
                                roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
                                'save_fig', 1, 'fig_ftype', 'fig');
end

%% Plot CNI sig elecs
stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
an_opts     = {'HGh_S_zbtS_trl2to151_fLog_sm0_stat15','HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1'};
% an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
tbin_id    = 'cnts';
hemi_opts = {'r','l'};
roi_id   = 'gROI';
atlas_id = 'Dx';
reg_type = 'v';
plot_out = 0;
show_labels = 0;

for an_ix = 1:numel(an_opts)
    for hemi_ix = 1:numel(hemi_opts)
        fn_view_recon_atlas_grp_stat(SBJs, pipeline_id, stat_id, an_opts{an_ix}, reg_type, show_labels,...
                                 hemi_opts{hemi_ix}, atlas_id, roi_id, plot_out)
%         fn_view_recon_atlas_grp_stat_onset(SBJs, pipeline_id, stat_id, an_opts{an_ix}, reg_type, show_labels,...
%                                  hemi_opts{hemi_ix}, atlas_id, roi_id, tbin_id)
%         fn_view_recon_stat(SBJs{s}, pipeline_id, stat_id, an_opts{an_ix}, 'pat', '', show_labels, hemi_opts{hemi_ix}, plot_out);
    end
end

%% Plot CNI sig onsets on recon
stat_id  = 'corrRT_CNI_pcon_WL200_WS50';
an_id    = 'HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1';
tbin_id  = 'cnts';
reg_type = 'v';
show_labels = 0;
atlas_id = 'Dx';
roi_id        = 'gROI';
plot_roi_opts = {{'l','lat'},{'r','lat'},{'l','MPFC'},{'r','MPFC'},{'l','deep'},{'r','deep'},{'b','OFC'}};

for roi_ix = 1:numel(plot_roi_opts)
    fn_view_recon_atlas_grp_ROI_stat_onset(SBJs, pipeline_id, stat_id, an_id, reg_type, show_labels,...
                                        plot_roi_opts{roi_ix}{1}, atlas_id, roi_id, plot_roi_opts{roi_ix}{2}, tbin_id,...
                                        'save_fig', 1, 'fig_ftype', 'svg');
end

%% HG active examples
SBJ         = 'IR74';
conditions  = 'CNI';
pipeline_id = 'main_ft';
actv_win    = '100';
an_id_s     = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
an_id_r     = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
plt_id      = 'stack_S2to15_R5to10_evnt_c5';
save_fig    = 1;
fig_vis     = 'on';
fig_ftype= 'svg';

SBJ10b_ANOVA_plot_SR_RTcorr(SBJ,stat_id,an_id_s,an_id_r,plt_id,save_fig,fig_vis,fig_ftype)
%% Plot ANOVA with RT correlation
SBJ = 'IR35';
stat_id = 'corrRT_CNI_pcon_WL200_WS50';
an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
SBJ10b_ANOVA_plot_SR_RTcorr(SBJ,stat_id,an_id_s,an_id_r,plt_id,1,'on','svg');

% OLD SHIT:
% SBJ08b_HFA_plot_SR_stack_cond_onset_noLab(SBJ,conditions,an_id_s,an_id_r,...
%     pipeline_id,actv_win,plt_id,save_fig,fig_vis,fig_ftype)
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

%% Proportions of significant effects across ROI
stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
atlas_id    = 'Dx';
an_opts     = {'HGh_S_zbtS_trl2to151_fLog_sm0_stat15','HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1'};
% an_opts     = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
actv_win    = 100;
plt_id      = 'onsets_trl0to15_evnt_roi';
roi_id      = 'gROI';%{'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
gm_thresh   = 0;
plot_out    = 0;
plot_scat   = 0;
fig_vis     = 'on';
fig_ftype   = 'svg';

for an_ix = 1:2
    % with SBJ scatter (for myself to understand)
    SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
        an_opts{an_ix},actv_win,roi_id,atlas_id,gm_thresh,plt_id,plot_out,1,save_fig,fig_vis,fig_ftype);
    % with just errbr (for nice plot)
    SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
        an_opts{an_ix},actv_win,roi_id,atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype);
    % with CSE (no results with HGh! try again later with HGm...)
%     SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_CSE_RT_ANOVA_ROI(SBJs,stat_id,pipeline_id,...
%         an_opts{an_ix},roi_id,atlas_id,gm_thresh,plt_id,plot_out,plot_scat,save_fig,fig_vis,fig_ftype);
end

%% Plot onsets of ANOVA+RT
pipeline_id = 'main_ft';
stat_id     = 'corrRT_CNI_pcon_WL200_WS50';
% clust_id    = 'kmeans_corr_nROI_itr1k';%'kmeans_corr_nCH_itr1k_srCmb',
tbin_id     = 'cnts';
an_opts     = {'HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1'};% stim?
gm_thresh   = 0;
median_yn   = 0;
roi_opts    = {'gROI'};%,'Yeo7'};%{'gROI','thryROI'};%,'LPFC','MPFC','INS','OFC','thryROI'};
atlas_opts  = {'Dx'};
% plt_opts    = {{'onsets_trl5to1_evnt_roi'}};%{'onsets_trl0to15_evnt_roi'},
plt_opts    = {{'onsets_trl0to15_violin_allROI','onsets_trl5to1_violin_all'}};%{{'onsets_trl0to15_violin_allSBJ','onsets_trl0to15_violin_allROI','onsets_trl0to15_violin_avgROI'},...
%                {'onsets_trl5to1_violin_allSBJ','onsets_trl5to1_violin_allROI','onsets_trl5to1_violin_avgROI'}};
fig_ftype = 'svg';%'png';
an_ix = 1;
for roi_ix = 1%:numel(roi_opts)
    for plt_ix = 1:numel(plt_opts)
            fprintf('roi: %s; an: %s; plt: %s\n',roi_opts{roi_ix},an_opts{an_ix},plt_opts{plt_ix}{plt_ix});
%             % Violin Plots
            SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
                                                    atlas_opts{roi_ix},gm_thresh,plt_opts{an_ix}{plt_ix},1,'on',fig_ftype)
            % Pair differences
%             SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,pipeline_id,an_opts{an_ix},roi_opts{roi_ix},...
%                                                     atlas_id,0,plt_opts{1}{2},save_fig,fig_vis)
    end
end


