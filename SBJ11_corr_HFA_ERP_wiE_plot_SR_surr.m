function SBJ11_corr_HFA_acE_plot_mat_SR_surr(SBJ,pipeline_id,stat_id,an_id_s,an_id_r,plt_id,save_fig,fig_vis,fig_filetype)
% Plot times series of two effects in same electrode, with significant
% epochs of correlation highlighted
%   stats done via generation of a null distribution
error('never wrote this script, would require multiple measures within elec');
if ischar(save_fig); save_fig = str2num(save_fig); end

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);

[grp_lab, grp_colors, grp_style] = fn_group_label_styles(model_lab);
[rt_lab, rt_color, rt_style]     = fn_group_label_styles('RT');
event_lab = {'stim', 'resp'};

% Load RTs
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'),'trial_info');

% Load HFA data
actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id_s,'_mn',actv_win,'.mat');
tmp = load(actv_filename); hfa{1} = tmp.hfa; actv_ch{1} = tmp.actv_ch;
actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id_r,'_mn',actv_win,'.mat');
tmp = load(actv_filename); hfa{2} = tmp.hfa; actv_ch{2} = tmp.actv_ch;
clear tmp;

% Load ERP data
stats_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_s,'.mat');
stats_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_',conditions,'_ROI_',an_id_r,'.mat');
tmp = load(stats_filename1,'stat'); stat{1} = tmp.stat;
tmp = load(stats_filename2,'stat'); stat{2} = tmp.stat;
tmp = load(stats_filename1,'roi_erp');
for cond_ix = 1:numel(cond_lab)
    erp{1,cond_ix} = tmp.roi_erp{cond_ix};
end
tmp = load(stats_filename2,'roi_erp');
for cond_ix = 1:numel(cond_lab)
    erp{2,cond_ix} = tmp.roi_erp{cond_ix};
end
clear tmp

%!!! is this the best way to do this??? Maybe not...
sample_rate = (numel(stat{1}.time)-1)/(stat{1}.time(end)-stat{1}.time(1));

% Load ROI and GM/WM info
einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
load(einfo_filename);
% Electrode Info Table:
%   label- name of electrode
%   ROI- specific region
%   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%   ROI2- specific region of second electrode
%   tissue- primary tissue type
%   GM weight- percentage of electrode pair in GM
%   Out- 0/1 flag for whether this is partially out of the brain

%% Correlation

end
