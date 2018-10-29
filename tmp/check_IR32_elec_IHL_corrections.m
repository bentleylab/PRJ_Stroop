[root_dir, app_dir] = fn_get_root_dir(); ft_dir = [app_dir 'fieldtrip/'];
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
view_angle = [-90 0];
pipeline_id = 'main_ft';
show_labels = 1;
hemi = 'l';

mesh = ft_read_headshape([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_cortex_lh.mat']);
%mni: load([ft_dir 'template/anatomy/surface_pial_left.mat']);
%fsavg: mesh = ft_read_headshape([root_dir 'PRJ_Stroop/data/atlases/freesurfer/fsaverage/' hemi 'h.pial']);
%%
elec_fname = 'IR32_elec_main_ft_mni_v.mat';
% view_space = 'pat';
% reg_type = '';
% atlas_name = 'Dx';
% roi_id = 'ROI';

%%
load(elec_fname);
cfgs = [];
cfgs.channel = {'IHL*'};
elec = fn_select_elec(cfgs,elec);
% if strcmp(reg_type,'v') || strcmp(reg_type,'s')
%     reg_suffix = ['_' reg_type];
% else
%     reg_suffix = '';
% end
% if strcmp(roi_id,'tissue') || strcmp(roi_id,'tissueC')
%     tis_suffix = '_tis';
% else
%     tis_suffix = '';
% end
% [roi_list, roi_colors] = fn_roi_label_styles(roi_id);

%%
h = figure;
mesh_alpha = 0.2;
ft_plot_mesh(mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
ft_plot_sens(elec, 'elecshape', 'point', 'marker', '.','elecsize', 10, 'label', 'label');
view(view_angle); material dull; lighting gouraud;
l = camlight;
set(h, 'windowkeypressfcn',   @cb_keyboard);

%% Fix the bad ones:
fn_compile_elec_struct(SBJ,'main_ft','pat','');
fn_compile_elec_struct(SBJ,'main_ft','mni','v');
fn_compile_elec_struct(SBJ,'main_ft','mni','s');
% fn_save_elec_atlas(SBJ,'main_ft','pat','','DK');
% fn_save_elec_atlas(SBJ,'main_ft','pat','','Dx');
fn_save_elec_tissue(SBJ,'main_ft','pat','','Dx');
fn_save_elec_tissue(SBJ,'main_ft','pat','','DK');
fn_save_elec_atlas(SBJ,'main_ft','mni','v','Yeo7');
fn_save_elec_atlas(SBJ,'main_ft','mni','v','Yeo17');
