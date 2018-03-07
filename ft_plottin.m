%% Viewing Recons
% Read in MRI:
%   mri = ft_read_mri(<SBJ>_fsMR_<pre or post>tal.nii);
% Grids only: read in the surface mesh ending in .pial
%   cortex = ft_read_headshape('<SBJ>_xh.pial'); % x = l or r
% Load elec file: '_elec_tal_r.mat'
%     load(<SBJ>_elec_tal_r.mat);
%     NOTE: if no grids, '<SBJ>_elec_tal_r,mat' will not exist; use '...tal_f.mat', or'..._tal.mat'
% Depth view (3 orthogonal slices):
%     cfg = [];
%     cfg.elec = elec_struct;
%     ft_electrodeplacement(cfg, mri)
% Surface view:
%     ft_plot_mesh(cortex, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
%     view([-90 10]); lighting gouraud; camlight;
%     hs = ft_plot_sens(elec_struct, 'shape', 'sphere');

%% Interactive 2D plot

% Compute tfr
%TFR commands

% pair down the tfr
cfg = [];
cfg.frequency = [5.9980 10.9963];
cfg.avgoverfreq = 'yes';
cfg.latency = [0.2 0.5];
cfg.avgovertime = 'yes';
cfg.avgoverrpt = 'yes';
tfr_sel = ft_selectdata(cfg,tfr);

% Load mesh
mesh = ft_read_headshape(...
    '/home/knight/ecog/DATA_FOLDER/Irvine/IR57/3D_Images/freesurfer_pre-op/freesurfer/surf/lh.pial');
mesh.coorsys = 'acpc';

% create layout
cfg =[];
cfg.headshape = mesh;
cfg.projection = 'orthographic';
cfg.channel = {'LAC*'};
cfg.viewpoint = 'left';
cfg.mask = 'convex';
cfg.boxchannel = {'LAC1','LAC2'};
lay = ft_prepare_layout(cfg,tfr);

cfg = [];
cfg.layout = lay;
cfg.showoutline = 'yes';
cfg.baseline = [-0.2 -0.05];
cfg. baselinetype = 'relchange';
ft_multiplotTFR(cfg,tfr);

%% 3D Mesh Plotting (for depths)
%Create meask
atlas = ft_read_atlas(...
    '/home/knight/ecog/DATA_FOLDER/Irvine/IR57/3D_Images/freesurfer_pre-op/freesurfer/mri/aparc+aseg.mgz');
atlas.coordsys = 'acpc';
cfg = [];
cfg.inputcoord = 'acpc';
cfg.atlas      = atlas;
cfg.roi        = {'Right-Hippocampus','Right-Amygdala'};
mask = ft_volumelookup(cfg,atlas);

% Convert to mesh
seg = keepfields(atlas,{'dim','unit','coordsys','transform'});
seg.brain = mask;

cfg = [];
cfg.method  = 'iso2mesh';
cfg.radbound = 2;
cfg.maxsurf  = 0;
cfg.tissue   = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3; 
mesh = ft_prepare_mesh(cfg,seg);

% plot mesh
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim  = 'maxabs';
cfg.method       = 'cloud';
% cfg.slice     
% cfg.nslices
cfg.facealpha    = 0.25;
ft_sourceplot(cfg, tfr_sel, mesh);

% select_elec function in ecogtmp/reconfunciton will do ft_selectdata for
% elec structures
% gets the shift for bipolar electrodes iff the .elec field is there in the
% data when you do the rereferencing with ft_preprocessing
