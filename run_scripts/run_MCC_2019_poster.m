if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
%%
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Set up variables to enter function
SBJs = {'CP24','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};
proc_id = 'main_ft';

%% ================================================================================
%  Design
%  ================================================================================
for s = 1:numel(SBJs)
    plot_design_statistics(SBJs{s});
end

%% ================================================================================
% %   BEHAVIOR
% %  =================================================================================
plt_id   = 'rt_hist';
save_fig = 1;
fig_vis  = 'on';
fig_ftype = 'png';
for s = 1:numel(SBJs)
    B00a_SBJ_RT_violins(SBJs{s},plt_id,fig_vis,save_fig,fig_ftype);
end

% RT behavior group level
B00b_GRP_RT_violins_norm(SBJs,plt_id,save_fig,fig_vis,fig_ftype);

%% Plot stacks with mean traces
% conditions  = 'actv';%'CNI';%
stat_id_b = 'PC_B2t0';%pCNI_
% stat_id_s = 'PCi_S0tmRT_WL1_WS50';%{'CNI_PC_S0tmRT_WL1_WS25',};
% stat_id_sd= 'PCi_S0tmRT_WL1_WS50_D1tRT';
% stat_id_r = 'PCi_R1t5_WL1_WS50';%{'CNI_PC_R1t5_WL1_WS25',};
stat_id_s = 'CNI_PC_S0tmRT_WL1_WS25';%{'CNI_PC_S0tmRT_WL1_WS50',};
stat_id_d = 'CNI_PC_D1tRT';%WS50
stat_id_sd= 'CNI_PC_S0tmRT_WL1_WS25_D1tRT';%WS50
stat_id_r = 'CNI_PC_R1t5_WL1_WS25';%{'CNI_PC_R1t5_WL1_WS50',};
proc_id   = 'main_ft';
an_id_s   = 'HGm_S2t151_zbtA_sm0_wn100';%wn200 --> nope
an_id_d   = 'HGm_S2t251_zbtA_sm0_wn100';%wn200 --> nope
an_id_r   = 'HGm_R5t101_zbtA_sm0_wn100';%wn200 --> nope
actv_id_b = 'actv_B2t0_mn50';
actv_id_s = 'actv_S0tmRT_mn1';
actv_id_d = 'actv_D1tRT';
actv_id_r = 'actv_R1t5_mn1';
% an_id_s   = 'HGm_S_zbtA_trl2to151_sm0_wn100_stat15';%'HGh_S_zbtS_trl2to151_fLog_sm0_stat15';
% an_id_r   = 'HGm_R_zbtA_trl5to101_sm0_wn100_stat5to1';%'HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1';
plt_id_s  = 'ERPstack_S2to15_evnt_c5';%'ts_S15R1_errbr_evnt';
plt_id_sr = 'ERPstack_S2to15_R5to10_evnt_c5';%'ts_S15R1_errbr_evnt';
save_fig  = 1;
fig_vis   = 'off';
fig_ftype = 'png';
for s = 1:numel(SBJs)
%         SBJ08ab_HFA_actv(SBJs{s},an_id_s,actv_id_b);
%         SBJ08ab_HFA_actv(SBJs{s},an_id_s,actv_id_s);
%         SBJ08ab_HFA_actv(SBJs{s},an_id_d,actv_id_d);
%         SBJ08ab_HFA_actv(SBJs{s},an_id_r,actv_id_r);
%     SBJ08b_HFA_plot_SR_ERPstack_cond(SBJs{s},'actv',an_id_s,an_id_r,actv_id_s,actv_id_r,...
%         plt_id_sr,save_fig,fig_vis,fig_ftype)
%     close all;
    SBJ08b_HFA_plot_SR_ERPstack_cond(SBJs{s},'CNI',an_id_s,an_id_r,stat_id_sd,stat_id_r,...
        plt_id_sr,save_fig,fig_vis,fig_ftype);
    close all;
%     SBJ08b_HFA_plot_SR_ERPstack_cond(SBJs{s},'PC',an_id_s,an_id_r,stat_id_b,stat_id_r{2},...
%         plt_id,save_fig,fig_vis,fig_ftype)
%     close all;
%     SBJ08b_HFA_plot_ERPstack_cond(SBJs{s},'pCNI',an_id_s,stat_id_b,...
%         plt_id_s,save_fig,fig_vis,fig_ftype)
%     close all;
    % PC - B
%     SBJ08b_HFA_plot_ERPstack_cond(SBJs{s},'PC',an_id_s,stat_id_b,...
%         plt_id_s,save_fig,fig_vis,fig_ftype)
    close all;
end

%% ================================================================================
%   ANATOMY
%  =================================================================================
%% Plot group recon with mgROI
% fn_view_recon_atlas_grp(SBJs,proc_id,'v',0,'l','Dx','mgROI',0);
% fn_view_recon_atlas_grp(SBJs,proc_id,'v',0,'r','Dx','mgROI',0);

atlas_id = 'Dx';
roi_id = 'gROI';
reg_type = 'v';
show_lab = 0;
roi_opts = {{'l','PAR'},{'r','PAR'}};%{{'l','lat'},{'r','lat'},{'l','MPFC'},{'r','MPFC'},{'l','deep'},{'r','deep'},{'b','OFC'}};

for roi_ix = 1:numel(roi_opts)
    fn_view_recon_atlas_grp_ROI(SBJs, proc_id, reg_type, show_lab,...
                                roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
                                'save_fig', 0, 'fig_ftype', 'fig');
end

%% ================================================================================
%   RECONS with EFFECTS
%  =================================================================================
%% Combine S + D analyses
stat_id_s  = 'CNI_PC_S0tmRT_WL1_WS25';%,'CNI_PC_R1t5_WL1_WS50'};%
stat_id_d  = 'CNI_PC_D1tRT';
stat_id_r  = 'CNI_PC_R1t5_WL1_WS25';
stat_id_sd = 'CNI_PC_S0tmRT_WL1_WS25_D1tRT';
% stat_id_s  = 'PCi_S0tmRT_WL1_WS50';%,'CNI_PC_R1t5_WL1_WS50'};%
% stat_id_d  = 'PCi_D1tRT';
% stat_id_r  = 'PCi_R1t5_WL1_WS50';
% stat_id_sd = 'PCi_S0tmRT_WL1_WS50_D1tRT';
an_id_s    = 'HGm_S2t151_zbtA_sm0_wn100';
an_id_d    = 'HGm_S2t251_zbtA_sm0_wn100';
an_id_r    = 'HGm_R5t101_zbtA_sm0_wn100';
for s = 1:numel(SBJs)
    SBJ10ab_combine_stats(SBJs{s},an_id_s,an_id_d,stat_id_s,stat_id_d);
    SBJ10ab_combine_stats(SBJs{s},an_id_s,an_id_r,stat_id_sd,stat_id_r);
end

%% Plot CNI sig elecs
% stat_id  = 'corrRT_CNI_PC_WL200_WS50';
% % an_opts     = {'HGh_S_zbtS_trl2to151_fLog_sm0_stat15','HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1'};
% an_opts  = {'HGm_S_zbtS_trl2to151_sm0_wn100_stat15','HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
% % an_opts  = {'HGm_S_zbtA_trl2to151_sm0_wn100_stat15','HGm_R_zbtA_trl5to101_sm0_wn100_stat5to1'};
% tbin_id    = 'cnts';
% hemi_opts = {'r','l'};
% roi_id   = 'gROI';
% atlas_id = 'Dx';
% reg_type = 'v';
% plot_out = 0;
% show_labels = 0;
% 
% for an_ix = 1%:numel(an_opts)
%     for hemi_ix = 1:numel(hemi_opts)
%         fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_id, an_opts{an_ix}, reg_type, show_labels,...
%                                  hemi_opts{hemi_ix}, atlas_id, roi_id, plot_out)
%         fn_view_recon_atlas_grp_stat_venn(SBJs, proc_id, stat_id, an_opts{an_ix}, reg_type, show_labels,...
%                                  hemi_opts{hemi_ix}, atlas_id, roi_id, plot_out)
% %         fn_view_recon_atlas_grp_stat_onset(SBJs, proc_id, stat_id, an_opts{an_ix}, reg_type, show_labels,...
% %                                  hemi_opts{hemi_ix}, atlas_id, roi_id, tbin_id)
% %         fn_view_recon_stat(SBJs{s}, proc_id, stat_id, an_opts{an_ix}, 'pat', '', show_labels, hemi_opts{hemi_ix}, plot_out);
%     end
% end

%% Plot overlap recons
stat_conds = {...
    {'PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'},...
    {'PCi_S0tmRT_WL1_WS50_D1tRT','HGm_S2t151_zbtA_sm0_wn100','PC'},...
    {'CNI_PC_S0tmRT_WL1_WS50_D1tRT_R1t5_WL1_WS50','HGm_S2t151_zbtA_sm0_wn100','CNI'},...
    };
%     {'CNI_PC_S0tmRT_WL1_WS50','HGm_S2t151_zbtA_sm0_wn100','CNI'},...
%     {'CNI_PC_D1tRT','HGm_S2t251_zbtA_sm0_wn100','CNI'},...
%     {'CNI_PC_R1t5_WL1_WS50','HGm_R5t101_zbtA_sm0_wn100','CNI'}...
%     {'PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'},...
%     {'CNI_PC_S0tmRT_WL1_WS50_D1tRT_R1t5_WL1_WS50','HGm_S2t151_zbtA_sm0_wn100','CNI'},...
%     {'CNI_PC_S0tmRT_WL1_WS50_D1tRT','HGm_S2t151_zbtA_sm0_wn100','CNI'},...

proc_id   = 'main_ft';
roi_id    = 'gROI';
atlas_id  = 'Dx';
plt_id    = 'venn';
reg_type  = 'v';
plot_out  = 0;
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'png';

% % Plot Venn Diagrams
% SBJ10c_HFA_GRP_venn_stats_ROI(SBJs, proc_id, stat_conds, 'b', atlas_id, roi_id,...
%                               plot_out, plt_id, save_fig, fig_ftype)

% Plot Recons
hemi_opts = {'r','l'};
for hemi_ix = 1:numel(hemi_opts)
    fn_view_recon_atlas_grp_stat_venn(SBJs, proc_id, stat_conds, reg_type, show_lab,...
        hemi_opts{hemi_ix}, atlas_id, roi_id, plot_out)
end

%% Proportions of significant effects across ROI
stat_conds = {...
    {'PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'},...
    {'PCi_S0tmRT_WL1_WS50_D1tRT','HGm_S2t151_zbtA_sm0_wn100','PC'},...
    {'CNI_PC_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'},...
    {'CNI_PC_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','CNI'},...
    {'CNI_PC_R1t5_WL1_WS25','HGm_R5t101_zbtA_sm0_wn100','CNI'}...
%               {'actv_D1tRT','HGm_S2t251_zbtA_sm0_wn100','actv'},...
%     {'CNI_PC_S0tmRT_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'},...
%     {'CNI_PC_D1tRT','HGm_S2t251_zbtA_sm0_wn100','CNI'},...
    };

proc_id     = 'main_ft';
atlas_id    = 'Dx';
plt_id      = 'errbar';
roi_id      = 'MPFC';
z_thresh    = 0;%1.5;
gm_thresh   = 0;
plot_out    = 0;
plot_scat   = 0;
fig_vis     = 'on';
fig_ftype   = 'png';
save_fig    = 1;

% with SBJ scatter (for myself to understand)
SBJ10c_HFA_GRP_errbar_ROI_stat_comb(SBJs,proc_id,stat_conds,atlas_id,roi_id,...
            gm_thresh,z_thresh,plt_id,plot_out,1,save_fig,fig_vis,fig_ftype);
% Pretty plot (no scatter)
SBJ10c_HFA_GRP_errbar_ROI_stat_comb(SBJs,proc_id,stat_conds,atlas_id,roi_id,...
            gm_thresh,z_thresh,plt_id,plot_out,0,save_fig,fig_vis,fig_ftype);

%% Redux (old version of SBJ10c bar plots
% stat_id = 'CNI_PC_S0tmRT_WL1_WS50';
% an_id   = 'HGm_S2t151_zbtA_sm0_wn100';
% actv_id = 'actv_S0tmRT_mn1';
% SBJ10c_HFA_GRP_summary_errbar_perc_GMlim_actv_RT_ANOVA_ROI(SBJs,proc_id,an_id,stat_id,actv_id,atlas_id,roi_id,...
%                                                             gm_thresh,plt_id,plot_out,1,save_fig,fig_vis,fig_ftype)
%% ROLs
% proc_id = 'main_ft';
% an_id = 'HGm_S_zbtS_trl2to251_sm0_wn100_stat0';
% rol_id = 'rol';
% plot_qa_trl = 1;
% plot_qa_summary = 1;
% plot_stack =1;
% fig_ftype = 'png';
% fig_vis = 'off';
% 
% for s = 1:numel(SBJs)
%     SBJ08d_HFA_ROL(SBJs{s}, an_id, rol_id, plot_qa_trl, plot_qa_summary, plot_stack, fig_vis, fig_ftype);
% end

%% Plot CNI sig onsets on recon
% stat_id  = 'corrRT_CNI_PC_WL200_WS50';
% an_id    = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';%,'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1'};
% tbin_id  = 'cnts';
% reg_type = 'v';
% show_labels = 0;
% atlas_id = 'Dx';
% roi_id        = 'gROI';
% plot_roi_opts = {{'l','lat'},{'r','lat'},{'l','MPFC'},{'r','MPFC'},{'l','deep'},{'r','deep'},{'b','OFC'}};
% 
% % individual ROIs
% for roi_ix = 1:numel(plot_roi_opts)
%     fn_view_recon_atlas_grp_ROI_stat_onset(SBJs, proc_id, stat_id, an_id, reg_type, show_labels,...
%                                         plot_roi_opts{roi_ix}{1}, atlas_id, roi_id, plot_roi_opts{roi_ix}{2}, tbin_id,...
%                                         'save_fig', 1, 'fig_ftype', 'svg');
% end
% 
% % all ROIs
% fn_view_recon_atlas_grp_stat_onset(SBJs, proc_id, stat_id, an_id, reg_type, show_labels,...
%                                  'r', atlas_id, roi_id, tbin_id)
% fn_view_recon_atlas_grp_stat_onset(SBJs, proc_id, stat_id, an_id, reg_type, show_labels,...
%                                  'l', atlas_id, roi_id, tbin_id)
                             
%% smANOVA time series by ROI
stat_id_s = 'CNI_PC_S0tmRT_WL1_WS25';%{'CNI_PC_S0tmRT_WL1_WS50',};
stat_id_d = 'CNI_PC_D1tRT';%WS50
stat_id_sd= 'CNI_PC_S0tmRT_WL1_WS25_D1tRT';%WS50
stat_id_r = 'CNI_PC_R1t5_WL1_WS25';%{'CNI_PC_R1t5_WL1_WS50',};
an_id_s   = 'HGm_S2t151_zbtA_sm0_wn100';%wn200 --> nope
an_id_d   = 'HGm_S2t251_zbtA_sm0_wn100';%wn200 --> nope
an_id_r   = 'HGm_R5t101_zbtA_sm0_wn100';%wn200 --> nope
plt_id_s  = 'ERPstack_S2to15_evnt_c5';%'ts_S15R1_errbr_evnt';
plt_id_sr = 'ERPstack_S2to15_R5to10_evnt_c5';%'ts_S15R1_errbr_evnt';
gm_thresh = 0;
z_thresh  = 0;
atlas_id  = 'Dx';
plot_nsig = 0;
save_fig  = 1;
fig_vis   = 'on';
fig_ftype = 'png';

for s = 1:numel(SBJs)
    SBJ10b_ANOVA_plot_SR_RTcorr_ROIcomb(SBJs{s},proc_id,stat_id_sd,stat_id_r,an_id_s,an_id_r,...
                                        atlas_id,gm_thresh,z_thresh,plot_nsig,plt_id_sr,...
                                        save_fig,fig_vis,fig_ftype)
%     SBJ10b_ANOVA_plot_SR_RTcorr_gROIcomb(SBJs{s},proc_id,stat_id_sd,stat_id_r,an_id_s,an_id_r,...
%                                         atlas_id,gm_thresh,z_thresh,plot_out,plt_id_sr,...
%                                         save_fig,fig_vis,fig_ftype)
end

%% HG active examples
% SBJ         = 'IR74';
% conditions  = 'CNI';
% proc_id = 'main_ft';
% actv_win    = '100';
% an_id_s     = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r     = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% plt_id      = 'stack_S2to15_R5to10_evnt_c5';
% save_fig    = 1;
% fig_vis     = 'on';
% fig_ftype= 'svg';
% 
% SBJ10b_ANOVA_plot_SR_RTcorr(SBJ,stat_id,an_id_s,an_id_r,plt_id,save_fig,fig_vis,fig_ftype)
% %% Plot ANOVA with RT correlation
% SBJ = 'IR35';
% stat_id = 'corrRT_CNI_PC_WL200_WS50';
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
% SBJ10b_ANOVA_plot_SR_RTcorr(SBJ,stat_id,an_id_s,an_id_r,plt_id,1,'on','svg');
% 
%% Plot onsets of ANOVA+RT
proc_id = 'main_ft';
stat_opts   = {'CNI_PC_S0tmRT_WL1_WS25_D1tRT','CNI_PC_R1t5_WL1_WS25'};
tbin_id     = 'cnts';
an_opts     = {'HGm_S2t151_zbtA_sm0_wn100','HGm_R5t101_zbtA_sm0_wn100'};%'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1','HGm_R_zbtA_trl5to101_sm0_wn100_stat5to1'};% stim?
z_thresh    = 1.5;
gm_thresh   = 0;
median_yn   = 0;
roi_opts    = {'MPFC'};
atlas_opts  = {'Dx'};
plt_opts    = {'onsets_trl5to1_violin_allROI','onsets_trl5to1_violin_all'};%{'onsets_trl0to15_evnt_roi'},
fig_ftype = 'png';
an_ix = 1;
for roi_ix = 1%:numel(roi_opts)
    for an_ix = 1:numel(an_opts)
            fprintf('roi: %s; an: %s; plt: %s\n',roi_opts{roi_ix},an_opts{an_ix},plt_opts{an_ix});
%             % Violin Plots
            SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_opts{an_ix},proc_id,an_opts{an_ix},roi_opts{roi_ix},...
                                                    atlas_opts{roi_ix},gm_thresh,z_thresh,plt_opts{an_ix},1,'on',fig_ftype)
            % Music staff + histogram of onsets
%             SBJ10c_HFA_GRP_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,proc_id,an_opts{an_ix},median_yn,roi_opts{roi_ix},...
%                                                     atlas_opts{roi_ix},gm_thresh,plt_opts{an_ix},1,'on',fig_ftype)
            % Pair differences
%             SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_id,proc_id,an_opts{an_ix},roi_opts{roi_ix},...
%                                                     atlas_id,0,plt_opts{1}{2},save_fig,fig_vis)
    end
end


