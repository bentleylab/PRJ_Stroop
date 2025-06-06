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

%% ================================================================================
%  Design
%  ================================================================================
% for s = 1:numel(SBJs)
%     plot_design_statistics(SBJs{s});
% end

%% ================================================================================
% %   BEHAVIOR
% %  =================================================================================
plt_id   = 'rt_hist_horz';
save_fig = 1;
fig_vis  = 'on';
fig_ftype = 'svg';
% for s = 1:numel(SBJs)
%     B00a_SBJ_RT_violins(SBJs{s},plt_id,fig_vis,save_fig,fig_ftype);
% end

% RT behavior group level
B00b_GRP_RT_violins_norm(SBJs,plt_id,save_fig,fig_vis,fig_ftype);

%% Plot HFA stacks with mean traces
% stat_id_b = 'PC_B2t0';%pCNI_
% % stat_id_s = 'PCi_S0tmRT_WL1_WS25';%{'CNI_PC_S0tmRT_WL1_WS25',};
% % stat_id_sd= 'PCi_S0tmRT_WL1_WS25_D1tRT';
% % stat_id_r = 'PCi_R1t5_WL1_WS25';%{'CNI_PC_R1t5_WL1_WS25',};
% stat_id_s = 'CNI_p_PC_S0tmRT_WL1_WS25';%{'CNI_PC_S0tmRT_WL1_WS50',};
stat_id_s15= 'CNI_p_PC_S0t15_WL1_WS25';%{'CNI_PC_S0tmRT_WL1_WS50',};
% stat_id_d = 'CNI_p_PC_D1tRT';%WS50
stat_id_sd= 'CNI_p_PC_S0tmRT_WL1_WS25_D1tRT';%WS50
stat_id_r15 = 'CNI_p_PC_R5t1_WL1_WS25';%{'CNI_PC_R1t5_WL1_WS50',};
stat_id_r = 'CNI_p_PC_R1t5_WL1_WS25';%{'CNI_PC_R1t5_WL1_WS50',};
proc_id   = 'main_ft';
an_id_s   = 'HGm_S2t151_zbtA_sm0_wn100';%wn200 --> nope
% an_id_d   = 'HGm_S2t251_zbtA_sm0_wn100';%wn200 --> nope
an_id_r   = 'HGm_R5t101_zbtA_sm0_wn100';%wn200 --> nope
% actv_id_b = 'actv_B2t0_mn50';
% actv_id_s = 'actv_S0tmRT_mn1';
% actv_id_d = 'actv_D1tRT';
% actv_id_r = 'actv_R1t5_mn1';
% % an_id_s   = 'HGm_S_zbtA_trl2to151_sm0_wn100_stat15';%'HGh_S_zbtS_trl2to151_fLog_sm0_stat15';
% % an_id_r   = 'HGm_R_zbtA_trl5to101_sm0_wn100_stat5to1';%'HGh_R_zbtS_trl5to101_fLog_sm0_stat5to1';
% plt_id_s  = 'ERPstack_S2to15_evnt_c5';%'ts_S15R1_errbr_evnt';
plt_id_sr = 'ERPstack_S2to15_R5to10_evnt_c5';%'ts_S15R1_errbr_evnt';%
atlas_id  = 'Dx';
roi_id    = 'gROI';
gm_thresh = 0;
z_thresh  = 0;
save_fig  = 0;
fig_vis   = 'on';
fig_ftype = 'svg';
% SBJ10b_ANOVA_plot_SR('IR35',proc_id,stat_id_s,stat_id_r,an_id_s,an_id_r,...
%                      atlas_id,gm_thresh,z_thresh,roi_id,plt_id_sr,...
%                      save_fig,fig_vis,fig_ftype)
% close all;
SBJ08b_HFA_plot_SR_ERPstack_cond('IR35','CI',an_id_s,an_id_r,stat_id_s15,stat_id_r15,...
                     plt_id_sr,save_fig,fig_vis,fig_ftype);
close all;
% for s = 1:numel(SBJs)
% %         SBJ08ab_HFA_actv(SBJs{s},an_id_s,actv_id_b);
% %         SBJ08ab_HFA_actv(SBJs{s},an_id_s,actv_id_s);
% %         SBJ08ab_HFA_actv(SBJs{s},an_id_d,actv_id_d);
% %         SBJ08ab_HFA_actv(SBJs{s},an_id_r,actv_id_r);
% %     SBJ08b_HFA_plot_SR_ERPstack_cond(SBJs{s},'actv',an_id_s,an_id_r,actv_id_s,actv_id_r,...
% %         plt_id_sr,save_fig,fig_vis,fig_ftype)
% %     close all;
%     SBJ08b_HFA_plot_SR_ERPstack_cond(SBJs{s},'CNI',an_id_s,an_id_r,stat_id_sd,stat_id_r,...
%         plt_id_sr,save_fig,fig_vis,fig_ftype);
%     close all;
% %     SBJ08b_HFA_plot_SR_ERPstack_cond(SBJs{s},'PC',an_id_s,an_id_r,stat_id_b,stat_id_r{2},...
% %         plt_id,save_fig,fig_vis,fig_ftype)
% %     close all;
% %     SBJ08b_HFA_plot_ERPstack_cond(SBJs{s},'pCNI',an_id_s,stat_id_b,...
% %         plt_id_s,save_fig,fig_vis,fig_ftype)
% %     close all;
%     % PC - B
% %     SBJ08b_HFA_plot_ERPstack_cond(SBJs{s},'PC',an_id_s,stat_id_b,...
% %         plt_id_s,save_fig,fig_vis,fig_ftype)
%     close all;
% end

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
% stat_id = 'corrRT_CNI_PC_WL200_WS25';
% an_id_s = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
% an_id_r = 'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1';
% plt_id  = 'ts_S0to15_R5to10_evnt_sigline';
% SBJ10b_ANOVA_plot_SR_RTcorr(SBJ,stat_id,an_id_s,an_id_r,plt_id,1,'on','svg');
% 
%% ================================================================================
%   ANATOMY
%  =================================================================================
%% Plot group recon with mgROI
% fn_view_recon_atlas_grp(SBJs,proc_id,'v',0,'l','Dx','mgROI',0);
% fn_view_recon_atlas_grp(SBJs,proc_id,'v',0,'r','Dx','mgROI',0);

roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
%,{'r','deep'},{'r','lat'},{'r','MPFC'}
proc_id   = 'main_ft';
atlas_id  = 'Dx';
reg_type  = 'v';
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'png';

for roi_ix = 1:numel(roi_opts)
    fn_view_recon_atlas_grp_ROI(SBJs, proc_id, reg_type, show_lab,...
                                roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
                                roi_opts{roi_ix}{3},'save_fig', save_fig, 'fig_ftype', fig_ftype);
end

%% Example HFA electrodes
SBJ = 'IR35';
plot_elecs = {'LAC2-3','RIN3-4'};
roi_opts  = {{'l','MPFC'},{'r','deep'}};
proc_id   = 'main_ft';
atlas_id  = 'Dx';

SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);
for e = 1:numel(plot_elecs)
    cfgs = []; cfgs.channel = plot_elecs(e);
    ch_elec = fn_select_elec(cfgs,elec);
    ch_elec.color = fn_roi2color(ch_elec.gROI);
    view_angle = fn_get_view_angle(roi_opts{e}{1},roi_opts{e}{2});
    
    atlas = fn_load_recon_atlas([],atlas_id);
    atlas_labels = fn_atlas_roi_select_mesh(atlas_id, roi_opts{e}{2}, roi_opts{e}{1});
    cfg = [];
    cfg.inputcoord = atlas.coordsys;
    cfg.atlas = atlas;
    cfg.roi = atlas_labels;
    roi_mask = ft_volumelookup(cfg,atlas);
    seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
    seg.brain = roi_mask;
    cfg = [];
    cfg.method      = 'iso2mesh';   % surface toolbox Arjen found
    cfg.radbound    = 2;            % scalar indicating the radius of the target surface mesh element bounding sphere
    cfg.maxsurf     = 0;
    cfg.tissue      = 'brain';
    cfg.numvertices = 100000;
    cfg.smooth      = 3;
    cfg.spmversion  = 'spm12';
    roi_mesh = ft_prepare_mesh(cfg, seg);
    figure;
    ft_plot_mesh(roi_mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', 0.3);
    ft_plot_sens(ch_elec, 'elecshape', 'sphere', 'facecolor', ch_elec.color, 'label', 'off');
    view(view_angle); material dull; lighting gouraud;
    l = camlight;
end

%% ================================================================================
%   RECONS with EFFECTS
%  =================================================================================
%% Combine S + D analyses
stat_id_s  = 'crCNI_p_PC_S0tmRT_WL1_WS25';%,'CNI_PC_R1t5_WL1_WS50'};%
stat_id_d  = 'crCNI_p_PC_D1tRT';
stat_id_r  = 'crCNI_p_PC_R1t5_WL1_WS25';
stat_id_sd = 'crCNI_p_PC_S0tmRT_WL1_WS25_D1tRT';
% stat_id_s  = 'pCI_PCn_S0tmRT_WL1_WS25';
% stat_id_d  = 'pCI_PCn_D1tRT';
% stat_id_r  = 'pCI_PCn_R1t5_WL1_WS25';
% stat_id_sd = 'pCI_PCn_S0tmRT_WL1_WS25_D1tRT';
an_id_s    = 'HGm_S2t151_zbtA_sm0_wn100';
an_id_d    = 'HGm_S2t251_zbtA_sm0_wn100';
an_id_r    = 'HGm_R5t101_zbtA_sm0_wn100';
for s = 1:numel(SBJs)
    SBJ10ab_combine_stats(SBJs{s},an_id_s,an_id_d,stat_id_s,stat_id_d);
    SBJ10ab_combine_stats(SBJs{s},an_id_s,an_id_r,stat_id_sd,stat_id_r);
end

%% Plot CNI sig elecs
proc_id   = 'main_ft';
stat_opts   = {
%                'pCNI_PC_B2t0'
%                'CNI_p_PC_S0tmRT_WL1_WS25'
%                'CNI_p_PC_D1tRT'
%                'CNI_p_PC_R1t5_WL1_WS25'
%                'CNI_p_PC_S0tmRT_WL1_WS25_D1tRT'
               'CNI_p_PC_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25'
               };
an_opts     = {'HGm_S2t151_zbtA_sm0_wn100'
%                'HGm_S2t151_zbtA_sm0_wn100'
%                'HGm_S2t251_zbtA_sm0_wn100'
%                'HGm_R5t101_zbtA_sm0_wn100'
% %                'HGm_S2t151_zbtA_sm0_wn100'
%                'HGm_S2t151_zbtA_sm0_wn100'
               };
hemi        = 'l';
roi_id      = 'gROI';
atlas_id    = 'Dx';
reg_type    = 'v';
mirror      = 1;
plot_out    = 0;
show_labels = 0;
save_fig    = 1;
fig_ftype   = 'png';

roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
for st_ix = 1:numel(stat_opts)
    %     % Lateral and medial stat recons
    %     fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_opts{st_ix}, an_opts{st_ix}, reg_type, show_labels,...
    %         hemi, atlas_id, roi_id, mirror, plot_out, save_fig, fig_ftype, 'view_angle', 'lat');
    %     fn_view_recon_atlas_grp_stat(SBJs, proc_id, stat_opts{st_ix}, an_opts{st_ix}, reg_type, show_labels,...
    %         hemi, atlas_id, roi_id, mirror, plot_out, save_fig, fig_ftype, 'view_angle', 'med');
    %     close all;
    for roi_ix = 1:numel(roi_opts)
        fn_view_recon_atlas_grp_stat_ROI(SBJs, proc_id, stat_opts{st_ix}, an_opts{st_ix}, ...
            reg_type, show_labels, roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2},...
            roi_opts{roi_ix}{3}, save_fig, fig_ftype);
    end
end

%         fn_view_recon_atlas_grp_stat_onset(SBJs, proc_id, stat_id, an_opts{an_ix}, reg_type, show_labels,...
%                                  hemi_opts{hemi_ix}, atlas_id, roi_id, tbin_id)
%         fn_view_recon_stat(SBJs{s}, proc_id, stat_id, an_opts{an_ix}, 'pat', '', show_labels, hemi_opts{hemi_ix}, plot_out);

%% Plot overlap recons
stat_conds = {...
    %     {'CNI_p_PC_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'}
    %     {'pCNI_PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','pCNI'}
%     {'pCNI_PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'}
    {'crCNI_p_PC_S0tmRT_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'}
    {'crCNI_p_PC_D1tRT','HGm_S2t251_zbtA_sm0_wn100','CNI'}
    {'crCNI_p_PC_R1t5_WL1_WS25','HGm_R5t101_zbtA_sm0_wn100','CNI'}
    };
%     {'pCI_PCn_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','pCI'}
proc_id   = 'main_ft';
roi_id    = 'gROI';
atlas_id  = 'Dx';
hemi      = 'b';
plt_id    = 'venn';
reg_type  = 'v';
plot_out  = 0;
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'svg';

% Plot Venn Diagrams
SBJ10c_HFA_GRP_venn_stats_ROI(SBJs, proc_id, stat_conds, hemi, atlas_id, roi_id,...
                              plot_out, plt_id, save_fig, fig_ftype);

% Plot lateral and medial stat venn recons
hemi      = 'l';
mirror    = 1;
plot_out  = 0;
fig_ftype = 'png';

fn_view_recon_atlas_grp_stat_venn(SBJs, proc_id, stat_conds, reg_type, show_lab,...
    hemi, atlas_id, roi_id, mirror, plot_out, save_fig, fig_ftype, 'view_angle', 'lat');
fn_view_recon_atlas_grp_stat_venn(SBJs, proc_id, stat_conds, reg_type, show_lab,...
    hemi, atlas_id, roi_id, mirror, plot_out, save_fig, fig_ftype, 'view_angle', 'med');

% Plot stat venn by specific ROI mesh
roi_opts  = {{'l','deep',1},{'l','lat',1},{'l','MPFC',1},{'b','OFC',0}};
proc_id   = 'main_ft';
atlas_id  = 'Dx';
reg_type  = 'v';
show_lab  = 0;
save_fig  = 1;
fig_ftype = 'png';

for roi_ix = 1:numel(roi_opts)
    fn_view_recon_atlas_grp_stat_venn_ROI(SBJs, proc_id, stat_conds, reg_type, show_lab,...
                                 roi_opts{roi_ix}{1}, atlas_id, roi_id, roi_opts{roi_ix}{2}, roi_opts{roi_ix}{3},...
                                 save_fig, fig_ftype)
end

%% Proportions of significant effects across ROI
stat_conds = {...
%     {'CNI_p_PC_S0tmRT_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'}...  % rCNI_p_PC for S
%     {'CNI_p_PC_D1tRT','HGm_S2t251_zbtA_sm0_wn100','CNI'}...            % rCNI_p_PC for D
%     {'CNI_p_PC_R1t5_WL1_WS25','HGm_R5t101_zbtA_sm0_wn100','CNI'}...    % rCNI_p_PC for R
    {'crCNI_p_PC_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'}
    {'crCNI_p_PC_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','pCNI'}
    {'crCNI_p_PC_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','PC'}
%     {'pCI_PCn_S0tmRT_WL1_WS25_D1tRT_R1t5_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','pCI'}
    {'pCNI_PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','pCNI'}
    {'pCNI_PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'}
    };
%     {'PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'}...
%     {'PCi_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','PC'}...
%     {'CNI_PC_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','CNI'}...
%     {'CNI_PC_R1t5_WL1_WS25','HGm_R5t101_zbtA_sm0_wn100','CNI'}...
%               {'actv_D1tRT','HGm_S2t251_zbtA_sm0_wn100','actv'}...
%     {'CNI_PC_S0tmRT_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'}...
%     {'CNI_PC_D1tRT','HGm_S2t251_zbtA_sm0_wn100','CNI'}...

proc_id     = 'main_ft';
atlas_id    = 'Dx';
plt_id      = 'errbar_y5';
roi_id      = 'gROI';
z_thresh    = 0;%1.5;
gm_thresh   = 0;
plot_out    = 0;
plot_scat   = 0;
fig_vis     = 'on';
fig_ftype   = 'svg';
save_fig    = 1;

% with SBJ scatter (for myself to understand)
% SBJ10c_HFA_GRP_errbar_ROI_stat_comb(SBJs,proc_id,stat_conds,atlas_id,roi_id,...
%             gm_thresh,z_thresh,plt_id,plot_out,1,save_fig,fig_vis,fig_ftype);
% Pretty plot (no scatter)
SBJ10c_HFA_GRP_errbar_ROI_stat_comb(SBJs,proc_id,stat_conds,atlas_id,roi_id,...
            gm_thresh,z_thresh,plt_id,plot_out,0,save_fig,fig_vis,fig_ftype);

%% Movies of significant effects
stat_conds = {...
    {'CNI_p_PC_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','CNI'}...
    {'CNI_p_PC_R1t5_WL1_WS25','HGm_R5t101_zbtA_sm0_wn100','CNI'}...
    };
%     actv:   {'actv_S0t150_mn1','HGm_S2t151_zbtA_sm0_wn100','actv'};
%     CNI.S:  {'CNI_PC_S0tmRT_WL1_WS25','HGm_S2t151_zbtA_sm0_wn100','CNI'};
%     CNI.SD: {'CNI_PC_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','CNI'};
%     CNI.R:  {'CNI_PC_R1t5_WL1_WS25','HGm_R5t101_zbtA_sm0_wn100','CNI'};
%     PB.B:   {'pCNI_PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'};

%     {'pCNI_PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','pCNI'},...
%     {'PC_B2t0','HGm_S2t151_zbtA_sm0_wn100','PC'},...
%     {'PCi_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','PC'},...
%     {'CNI_PC_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','CNI'},...
%               {'actv_D1tRT','HGm_S2t251_zbtA_sm0_wn100','actv'},...
%     ,...
%     {'CNI_PC_D1tRT','HGm_S2t251_zbtA_sm0_wn100','CNI'},...

proc_id     = 'main_ft';
atlas_id    = 'Dx';
plt_opts    = {'movie_S0tD' 'movie_R1t5'};%'movie_S0tD_sp20';%'movie_S0t15_sp20';
view_angles = {'med', 'lat'};
roi_id      = 'gROI';
hemi        = 'l';
mirror      = 1;
plot_out    = 0;
% z_thresh    = 0;%1.5;
% gm_thresh   = 0;
fig_vis     = 'on';

for s = [6 11 14 15 16 17]% 1 2 3 4 5 7 8 9 10 12 13]%1:numel(SBJs)
    for stat_ix = 1:numel(stat_conds)
        for v_ix = 1:numel(view_angles)
            fn_view_recon_stat_movie(SBJs{s}, proc_id, stat_conds{stat_ix},...
                atlas_id, roi_id, hemi, mirror, plot_out, plt_opts{stat_ix}, fig_vis,...
                'view_angle', view_angles{v_ix});
            close all;
        end
    end
end

stat_conds = {...
    {'rCNI_p_PC_S0tmRT_WL1_WS25_D1tRT','HGm_S2t151_zbtA_sm0_wn100','CNI'}...
    {'rCNI_p_PC_R1t5_WL1_WS25','HGm_R5t101_zbtA_sm0_wn100','CNI'}...
    };
for s = [6 11 14 15 16 17 1 2 3 4 5 7 8 9 10 12 13]%1:numel(SBJs)
    for stat_ix = 1:numel(stat_conds)
        for v_ix = 1:numel(view_angles)
            fn_view_recon_stat_movie(SBJs{s}, proc_id, stat_conds{stat_ix},...
                atlas_id, roi_id, hemi, mirror, plot_out, plt_opts{stat_ix}, fig_vis,...
                'view_angle', view_angles{v_ix});
            close all;
        end
    end
end
stat_conds = {...
    {'actv_S0t150_mn1','HGm_S2t151_zbtA_sm0_wn100','actv'}
    };
plt_opts = {'movie_S0t15'};
for s = [6 11 14 15 16 17]% 1 2 3 4 5 7 8 9 10 12 13]%1:numel(SBJs)
    for stat_ix = 1:numel(stat_conds)
        for v_ix = 1:numel(view_angles)
            fn_view_recon_stat_movie(SBJs{s}, proc_id, stat_conds{stat_ix},...
                atlas_id, roi_id, hemi, mirror, plot_out, plt_opts{stat_ix}, fig_vis,...
                'view_angle', view_angles{v_ix});
            close all;
        end
    end
end

%% Plot onsets of ANOVA+RT
proc_id = 'main_ft';
stat_opts   = {'CNI_p_PC_S0tmRT_WL1_WS25_D1tRT','CNI_p_PC_R1t5_WL1_WS25'};%{'pCI_PCn_S0tmRT_WL1_WS25_D1tRT','pCI_PCn_R1t5_WL1_WS25'};%
tbin_id     = 'cnts';
an_opts     = {'HGm_S2t151_zbtA_sm0_wn100','HGm_R5t101_zbtA_sm0_wn100'};%'HGm_R_zbtS_trl5to101_sm0_wn100_stat5to1','HGm_R_zbtA_trl5to101_sm0_wn100_stat5to1'};% stim?
z_thresh    = 0;
gm_thresh   = 0;
median_yn   = 0;
roi_opts    = {'gROI','PFC'};
atlas_id    = 'Dx';
plt_id      = 'onsets_trl5to1_violin_allROI_mirror';%,'onsets_trl5to1_violin_allROI'};%{'onsets_trl0to15_evnt_roi'},
fig_ftype   = 'svg';

for roi_ix = 1:numel(roi_opts)
    for st_ix = 1:numel(an_opts)
            fprintf('roi: %s; an: %s; plt: %s\n',roi_opts{roi_ix},an_opts{st_ix},plt_id);
            % Violin Plots
%             SBJ10c_HFA_GRPavg_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_opts{st_ix},proc_id,an_opts{st_ix},roi_opts{roi_ix},...
%                                                     atlas_id,gm_thresh,z_thresh,plt_id,1,'on',fig_ftype)
            % Pair differences
            SBJ10c_HFA_GRP_onsets_ROI_pairdiffs_ANOVA(SBJs,stat_opts{st_ix},proc_id,an_opts{st_ix},roi_opts{roi_ix},...
                                                    atlas_id)
    end
end

            % Music staff + histogram of onsets
%             SBJ10c_HFA_GRP_onsets_ROI_normRTout_RT_ANOVA(SBJs,stat_id,proc_id,an_opts{an_ix},median_yn,roi_opts{roi_ix},...
%                                                     atlas_opts{roi_ix},gm_thresh,plt_opts{an_ix},1,'on',fig_ftype)

%% ========================================================================
%   NOT IN MCC POSTER OR OXFORD TALK
%  ========================================================================

%% smANOVA time series by ROI
% stat_id_s = 'CNI_PC_S0tmRT_WL1_WS25';%{'CNI_PC_S0tmRT_WL1_WS25',};
% stat_id_d = 'CNI_PC_D1tRT';%WS50
% stat_id_sd= 'CNI_PC_S0tmRT_WL1_WS25_D1tRT';%WS50
% stat_id_r = 'CNI_PC_R1t5_WL1_WS25';%{'CNI_PC_R1t5_WL1_WS25',};
% an_id_s   = 'HGm_S2t151_zbtA_sm0_wn100';%wn200 --> nope
% an_id_d   = 'HGm_S2t251_zbtA_sm0_wn100';%wn200 --> nope
% an_id_r   = 'HGm_R5t101_zbtA_sm0_wn100';%wn200 --> nope
% plt_id_sr = 'ERPstack_S0to10_R1to5_evnt_c5';%'ts_S15R1_errbr_evnt';
% gm_thresh = 0;
% z_thresh  = 0;
% atlas_id  = 'Dx';
% plot_nsig = 0;
% save_fig  = 1;
% fig_vis   = 'off';
% fig_ftype = 'png';
% 
% for s = 1:numel(SBJs)
%     SBJ10b_ANOVA_plot_SR_RTcorr_ROIcomb(SBJs{s},proc_id,stat_id_sd,stat_id_r,an_id_s,an_id_r,...
%                                         atlas_id,gm_thresh,z_thresh,plot_nsig,plt_id_sr,...
%                                         save_fig,fig_vis,fig_ftype)
%     close all;
%     SBJ10b_ANOVA_plot_SR_RTcorr_gROIcomb(SBJs{s},proc_id,stat_id_sd,stat_id_r,an_id_s,an_id_r,...
%                                         atlas_id,gm_thresh,z_thresh,plot_nsig,plt_id_sr,...
%                                         save_fig,fig_vis,fig_ftype)
%     close all;
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
                             
%% ===============================================================================================================
% OLDER:
%=================================================================================================================
%% Redux (old version of SBJ10c bar plots
% stat_id = 'CNI_PC_S0tmRT_WL1_WS25';
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


