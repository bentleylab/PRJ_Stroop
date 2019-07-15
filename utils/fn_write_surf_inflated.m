function fn_write_surf_inflated(SBJ)
%% Write freesurfer inflated files to .mat
subdir = ['/home/knight/ecog/DATA_FOLDER/Irvine/' SBJ '/']; 
recondir = [subdir '3D_Images/Recon_Apr_2018/FT_Pipeline/'];

fshome = '/usr/local/freesurfer_x86_64-5.3.0'; % do NOT edit

addpath([fshome '/matlab']);
setenv('FREESURFER_HOME', fshome);
PATH = getenv('PATH');
setenv('PATH', [PATH ':' fshome '/bin']);

if exist([recondir 'Scans/' SBJ '_fsMR_pre_acpc.nii']) == 2
  scantype = '_pre-op';
elseif exist([recondir 'Scans/' SBJ '_fsMR_post_acpc.nii']) == 2
  scantype = '_post-op';
end

lh_inflated_path = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/lh.inflated'];
rh_inflated_path = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/rh.inflated'];

% Prepare Left and Right Cortical Surfaces
inflated_lh = ft_read_headshape(lh_inflated_path);
inflated_rh = ft_read_headshape(rh_inflated_path);

% Plot to check
figure;
ft_plot_mesh(inflated_lh);
material dull; lighting gouraud;
figure;
ft_plot_mesh(inflated_rh,'vertexcolor',repmat(inflated_rh.curv,[1 3]));
material dull; lighting gouraud;
l = camlight;

% Save Left and Right Cortical Surfaces
save([recondir 'Surfaces/' SBJ '_inflated_lh.mat'], 'inflated_lh');
save([recondir 'Surfaces/' SBJ '_inflated_rh.mat'], 'inflated_rh');

end
