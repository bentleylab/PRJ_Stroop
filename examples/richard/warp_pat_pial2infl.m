% Warp patient elec file to patient inflated space
cfg = [];
cfg.fshome = fshome;
cfg.elec = elec_acpc_r;
cfg.method = 'headshape';
cfg.headshape = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/lh.pial']; %FS path to hull
cfg.warp = 'fsinflated'; 
elec_acpc_ri = ft_electroderealign(cfg, elec_acpc_r); 

% Plot patient inflated space with elecs
cfg = [];
cfg.fshome = fshome;
cfg.elec = elec_acpc_r;
cfg.method = 'headshape';
cfg.headshape = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/lh.pial']; %FS path to hull
cfg.warp = 'fsinflated';
elec_acpc_ri = ft_electroderealign(cfg, elec_acpc_r);
