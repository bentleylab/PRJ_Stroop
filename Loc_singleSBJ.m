%% Viewing Recons
SBJ = 'IR35';
view_space = 'patient';     % {'patient', 'MNI_vol', 'MNI_srf'}
show_labels = true;         % {true, false}
plt_id = 'loc_SBJ_ROI';
an_id = 'HGm_S_zbtS_trl2to151_sm0_wn100_stat15';
actv_win = '100';
epoch_lim = [0.2 0.4];

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Set Up Variable
SBJ_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
% plt_vars_cmd = ['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m'];
% eval(plt_vars_cmd);

if strcmp(SBJ(1:2),'IR')
    im_dir = ['/home/knight/ecog/DATA_FOLDER/Irvine/' SBJ '/3D_Images/'];
    recon_name = dir([im_dir 'Recon_*']);
    recon_dir = [im_dir recon_name.name '/FT_Pipeline/'];
elseif strcmp(SBJ(1:2),'CP')
    im_dir = ['/home/knight/ecog/DATA_FOLDER/CPMC/' SBJ '/3D_Images/'];
    recon_name = dir([im_dir 'Recon_*']);
    recon_dir = [im_dir recon_name.name '/FT_Pipeline/'];
else
    error('what SBJ is this? no recon dir');
end
if ~exist([recon_dir '/Electrodes/' SBJ '_elec_acpc.mat']) && ~exist([recon_dir '/Electrodes/' SBJ '_elec_acpc_f.mat'])
  error('The recon has been started, but is not yet finished')
end
warning(['To adjust specific aspects of the 3d figure, such as transparency or adding electrode labels,'...
    'please explicitly use ft_plot_mesh (to plot the surface) and ft_plot_sens (to plot the electrodes)'...
    'and refer to the documentation for those functions.']);

%% Load elec
load IR35_elec_acpc_f

%% Load and prepare functional data
actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',an_id,'_mn',actv_win,'.mat');
load(actv_filename);
hfa.elec = elec_acpc_f;

% pair down the tfr
cfg = [];
cfg.frequency = hfa.freq;
cfg.avgoverfreq = 'yes';
cfg.latency = epoch_lim;
cfg.avgovertime = 'yes';
cfg.avgoverrpt = 'yes';
hfa_sel = ft_selectdata(cfg,hfa);

%% 3D depth plotting
roi_labels = {...
'ctx-caudalanteriorcingulate',...
'ctx-lh-caudalmiddlefrontal',...
'ctx-lh-isthmuscingulate',...
'ctx-lh-posteriorcingulate',...
'ctx-lh-paracentral',...
'ctx-lh-parsopercularis',...
'ctx-lh-parsorbitalis',...
'ctx-lh-parstriangularis',...
'ctx-lh-precentral',...
'ctx-lh-rostralmiddlefrontal',...
'ctx-lh-superiorfrontal',...
'ctx-lh-insula','ctx-lh-rostralanteriorcingulate','ctx-lh-frontalpole','ctx-lh-lateralorbitofrontal',...
'ctx-lh-medialorbitofrontal'};

%Create mask for atlas ROIs
atlas = ft_read_atlas([im_dir 'freesurfer_post-op/freesurfer/mri/aparc+aseg.mgz']);
atlas.coordsys = 'acpc';
cfg = [];
cfg.inputcoord = 'acpc';
cfg.atlas      = atlas;
cfg.roi        = roi_labels;
mask = ft_volumelookup(cfg,atlas);

% Convert volumetric mask to mesh
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
% cfg.slice     = '3d';
% cfg.nslices    = 2;
cfg.facealpha    = 0.25;
ft_sourceplot(cfg, hfa_sel, mesh);

% select_elec function in ecogtmp/reconfunciton will do ft_selectdata for
% elec structures
% gets the shift for bipolar electrodes iff the .elec field is there in the
% data when you do the rereferencing with ft_preprocessing


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

%% Mesh and Elec loading
mesh_alpha = 1;
if strcmp(view_space, 'patient')
    tmp = load([recon_dir 'Surfaces/' SBJ '_cortex_lh.mat']); mesh_l = tmp.cortex_lh;
    tmp = load([recon_dir 'Surfaces/' SBJ '_cortex_rh.mat']); mesh_r = tmp.cortex_rh;
    if exist([recon_dir 'Electrodes/' SBJ '_elec_acpc_fr.mat'])
        tmp = load([recon_dir 'Electrodes/' SBJ '_elec_acpc_fr.mat']); elec = tmp.elec_acpc_fr;
    elseif exist([recon_dir 'Electrodes/' SBJ '_elec_acpc_r.mat'])
        tmp = load([recon_dir 'Electrodes/' SBJ '_elec_acpc_r.mat']); elec = tmp.elec_acpc_r;
    elseif exist([recon_dir 'Electrodes/' SBJ '_elec_acpc.mat'])
        tmp = load([recon_dir 'Electrodes/' SBJ '_elec_acpc.mat']); elec = tmp.elec_acpc;
        mesh_alpha = 0.3;
    else
        error(['The recon has been started but the electrodes have not'...
            'been projected to the surface yet or there are no surface electrodes']);
    end
elseif strcmp(view_space, 'MNI_vol')
    tmp = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_left.mat'); mesh_l = tmp.mesh;
    tmp = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_right.mat'); mesh_r = tmp.mesh;
    if exist([recon_dir 'Electrodes/' SBJ '_elec_mni_v.mat'])
        tmp = load([recon_dir 'Electrodes/' SBJ '_elec_mni_v.mat']); elec = tmp.elec_mni_v;
    elseif exist([recon_dir 'Electrodes/' SBJ '_elec_mni_frv.mat'])
        tmp = load([recon_dir 'Electrodes/' SBJ '_elec_mni_frv.mat']); elec = tmp.elec_mni_frv;
    else
        error('The recon has been started but the electrodes have not been volume-based normalized yet')
    end
elseif strcmp(view_space, 'MNI_srf')
    fs_home = '/usr/local/freesurfer_x86_64-5.3.0/';
    mesh_l = ft_read_headshape([fs_home 'subjects/fsaverage/surf/lh.pial']);
    mesh_l.coordsys = 'tal';
    mesh_r = ft_read_headshape([fs_home 'subjects/fsaverage/surf/rh.pial']);
    mesh_r.coordsys = 'tal';
    if exist([recon_dir 'Electrodes/' SBJ '_elec_mni_s.mat'])
        tmp = load([recon_dir 'Electrodes/' SBJ '_elec_mni_s.mat']); elec = tmp.elec_mni_s;
    elseif exist([recon_dir 'Electrodes/' SBJ '_elec_fsavg_frs.mat'])
        tmp = load([recon_dir 'Electrodes/' SBJ '_elec_fsavg_frs.mat']); elec = tmp.elec_fsavg_frs;
    else
        error(['The recon has been started but the electrodes have either not been surface-based normalized'...
            'or there are only depth, in which case surface-based normalization is not possible']);
    end
else
    error(['unknown view_space: ' view_space]);
end

%% Plot mesh and electrodes
h = figure;

ft_plot_mesh(mesh_r, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
ft_plot_mesh(mesh_l, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'facealpha', mesh_alpha);
view([0 90]); lighting gouraud; camlight;      %[-55 10]
if show_labels
    ft_plot_sens(elec, 'elecshape', 'sphere', 'label', 'label');
else
    ft_plot_sens(elec, 'elecshape', 'sphere');
end

l = camlight; lighting gouraud; material dull;
fprintf(['To reset the position of the camera light after rotating the figure,\n' ...
    'make sure none of the figure adjustment tools (e.g., zoom, rotate) are active\n' ...
    '(i.e., uncheck them within the figure), and then hit ''l'' on the keyboard\n'])
set(h, 'windowkeypressfcn',   @cb_keyboard);

%% Interactive 2D plot
% Load mesh
mesh = ft_read_headshape(...
    [irvine_dir SBJ '/3D_Images/freesurfer_pre-op/freesurfer/surf/lh.pial']);
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
%% OLD STUFF

% Compute tfr
%TFR commands



