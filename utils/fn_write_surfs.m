function fn_write_surfs(SBJ)
%% Write freesurfer inflated files to .mat

% Set up directories
[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if strcmp(SBJ(1:2),'IR')
    hosp_code = 'Irvine';
elseif strcmp(SBJ(1:2),'CP')
    hosp_code = 'CPMC';
else
    error('unknown hospital code');
end

img_dir = ['/home/knight/ecog/DATA_FOLDER/' hosp_code '/' SBJ '/3D_Images/']; 
recon_str = dir([img_dir 'Recon_*']);
ft_pipe_dir = [img_dir '/' recon_str.name '/FT_Pipeline/'];

fshome = '/usr/local/freesurfer_x86_64-5.3.0'; % do NOT edit

addpath([fshome '/matlab']);
setenv('FREESURFER_HOME', fshome);
PATH = getenv('PATH');
setenv('PATH', [PATH ':' fshome '/bin']);

if exist([ft_pipe_dir 'Scans/' SBJ '_fsMR_pre_acpc.nii']) == 2
  scan_type = '_pre-op';
elseif exist([ft_pipe_dir 'Scans/' SBJ '_fsMR_post_acpc.nii']) == 2
  scan_type = '_post-op';
end

% Create inflated files
lh_infl_path = [img_dir 'freesurfer' scan_type '/freesurfer/surf/lh.inflated'];
rh_infl_path = [img_dir 'freesurfer' scan_type '/freesurfer/surf/rh.inflated'];

% Prepare Left and Right Cortical Surfaces
infl_lh = ft_read_headshape(lh_infl_path);
infl_rh = ft_read_headshape(rh_infl_path);

% Save Left and Right Cortical Surfaces
save([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_infl_lh.mat'], '-v7.3', 'infl_lh');
save([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_infl_rh.mat'], '-v7.3', 'infl_rh');

% Create smoothwm files
lh_wm_path = [img_dir 'freesurfer' scan_type '/freesurfer/surf/lh.smoothwm'];
rh_wm_path = [img_dir 'freesurfer' scan_type '/freesurfer/surf/rh.smoothwm'];

% Prepare Left and Right Cortical Surfaces
wm_lh = ft_read_headshape(lh_wm_path);
wm_rh = ft_read_headshape(rh_wm_path);

% Save Left and Right Cortical Surfaces
save([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_wm_lh.mat'], '-v7.3', 'wm_lh');
save([SBJ_vars.dirs.recon 'Surfaces/' SBJ '_wm_rh.mat'], '-v7.3', 'wm_rh');

%% Plot to check
% figure;
% ft_plot_mesh(inflated_lh);
% material dull; lighting gouraud;
% figure;
% ft_plot_mesh(inflated_rh,'vertexcolor','curv');
% material dull; lighting gouraud;
% l = camlight;

end
