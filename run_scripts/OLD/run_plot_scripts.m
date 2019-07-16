%% Set up variables to enter function
SBJs = {'IR21','IR31','IR32','IR35','IR39','IR41','IR54','IR57'};

SBJ         = 'IR41';
conditions  = 'CI';
pipeline_id = 'main_ft';
an_id       = 'HGm_S_zbootS_trl2to15_sm10_stat15';%'HGm_R_zbootS_trl5to1_sm10_stat5to1';
actv_win    = '100';
plt_id      = 'onsets_trl0to15_evnt_roi';%'summary_trl5to1_evnt_roi';
save_fig    = 1;
fig_vis     = 'on';

%% SBJ08c
for sbj_ix = 1:numel(SBJs)
    SBJ08c_HFA_summary_hist_actv_cond_ROI(SBJs{sbj_ix},conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis);
end
SBJ08c_HFA_GRP_summary_bar_perc_actv_cond_ROI(SBJs,conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis);

%%
% roi_list = {};
% for sbj_ix = 1:numel(SBJs)
%     einfo_filename = ['~/PRJ_Stroop/data/' SBJs{sbj_ix} '/02_preproc/' SBJs{sbj_ix} '_einfo_' pipeline_id '.mat'];
%     load(einfo_filename);
%     roi_list{sbj_ix} = unique(einfo(:,2));
% %     clear einfo
% %     SBJ08c_HFA_summary_plot_actv_cond_ROI(SBJs{sbj_ix},conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis);
% end
% all_rois = unique(vertcat(roi_list{:}));
all_rois = {'FPC','DLPFC','VLPFC','PM','M1','S1','ACC','preSMA','aMCC',...
            'SMA','pMCC','Precuneus','vaINS','daINS','FO','mINS','pINS',...
            'mOFC','lOFC','FWM','','OUT'};
% groi_list = {};
% for sbj_ix = 1:numel(SBJs)
%     einfo_filename = ['~/PRJ_Stroop/data/' SBJs{sbj_ix} '/02_preproc/' SBJs{sbj_ix} '_einfo_' pipeline_id '.mat'];
%     load(einfo_filename);
%     groi_list{sbj_ix} = unique(einfo(:,3));
% %     clear einfo
% %     SBJ08c_HFA_summary_plot_actv_cond_ROI(SBJs{sbj_ix},conditions,pipeline_id,an_id,actv_win,plt_id,save_fig,fig_vis);
% end
% all_grois = unique(vertcat(groi_list{:}));
all_grois = {'LPFC','MPFC','INS','OFC','FWM','OUT'};