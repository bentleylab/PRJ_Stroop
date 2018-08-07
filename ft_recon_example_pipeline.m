%% FT Recon Example (as in Stolk et al., 2018 Nat. Protocols)
SBJ_dir = '/home/knight/hoycw/ft_recon_example/';
SBJ = 'SubjectUCI29';

fs_home_dir = '/usr/local/freesurfer_x86_64-5.3.0/';
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% View SBJ freesurfer surfaces
pial_lh = ft_read_headshape([SBJ_dir 'freesurfer/surf/lh.pial']);
pial_lh.coordsys = 'acpc';
ft_plot_mesh(pial_lh);
lighting gouraud; camlight;

%% Electrode placement (click CT/MR image to get coordinates)
% Load fused CT
ct_acpc_f = ft_read_mri([SBJ_dir SBJ '_CT_acpc_f.nii']);
% load freesurfer processsed mri
fsmri_acpc = ft_read_mri([SBJ_dir 'freesurfer/mri/T1.mgz']);
fsmri_acpc.coordsys = 'acpc';
% load header
hdr = ft_read_header(raw_file);

% click to place electrodes - not doing this now....
% cfg = [];
% cfg.channel = hdr.label;
% elec_acpc_f = ft_electrodeplacement(cfg,ct_acpc_f, fsmri_acpc);

%% View MRI + Electrodes in ortho (3 slice) plot
load([SBJ_dir SBJ '_elec_acpc_f.mat']);
ft_plot_ortho(fsmri_acpc.anatomy, 'transform', fsmri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');

%% Brain shift compensation (separate for each grid)
% Create a hull to project the electrodes onto
% cfg = [];
% cfg.method = 'cortexhull';
% cfg.headshape = [SBJ_dir 'freesurfer/surf/lh.pial']);
% cfg.fshome = fs_home_dir;
% hull_lh = ft_prepare_mesh(cfg);
% then save(mesh);

% Apply warping algorithm
% elec_acpc_fr = elec_acpc_f;
% grids = {'LPG*','LTG*'};
% for g = 1:numel(grids)
%     cfg = [];
%     cfg.channel = grids{g};
%     cfg.keepchannel = 'yes';
%     cfg.elec = elec_acpc_fr;
%     cfg.method = 'headshape';
%     cfg.headshape = hull_lh;
%     cfg.warp = 'dykstra2012';
%     cfg.feedback = 'yes';
%     elec_acpc_fr = ft_electroderealign(cfg);
% end

% View pre-brain shift compensation
ft_plot_mesh(pial_lh);
ft_plot_sens(elec_acpc_f);
view([-55 10]); material dull; lighting gouraud; camlight;

% View the result after brain shift compensation
ft_plot_mesh(pial_lh);
ft_plot_sens(elec_acpc_fr);
view([-55 10]); material dull; lighting gouraud; camlight;
% save elec_acpc_fr

%% Volume based registration
% see manuscript
% fsmri_mni = ft_volumenormalise(cfg, fsmri_acpc);
% elec_minnfrv = ft_warp_apply(...);

% View the results in MNI space
load([ft_dir 'template/anatomy/surface_pial_left.mat']);    %also surface_pial_both.mat, etc.
ft_plot_mesh(mesh);
ft_plot_sens(elec_mni_frv);
view([-90 20]); material dull; lighting gouraud; camlight;

%% Surface based registration
% Map electrodes onto fsaverage brain surface
cfg = [];
cfg.channel = {'LPG*', 'LTG*'};
cfg.elec = elec_acpc_fr;
cfg.method = 'headshape';
cfg.headshape = [SBJ_dir 'freesurfer/surf/lh.pial'];
cfg.warp = 'fsaverage';
cfg.fshome = fs_home_dir;
elec_fsavg_frs = ft_electroderealign(cfg);

% View fsaverage surface with electrodes now aligned
fspial_lh = ft_read_headshape([fs_home_dir 'subjects/fsaverage/surf/lh.pial']);
fspial_lh.coordsys = 'fsaverage';
ft_plot_mesh(fspial_lh);
ft_plot_sens(elec_fsavg_frs);
view([-90 20]); material dull; lighting gouraud; camlight;

%% Find ROI Labels of electrodes (Create elec table)
% edit generate_electable.m % more extensive version

% Import atlas of interest
atlas = ft_read_atlas([ft_dir 'template/atlas/aal/ROI_M- NI_V4.nii']);
% look up label of electrode location
cfg = [];
cfg.roi = elec_mni_frv.chanpos(match_str(elec_mni_frv.label, 'LHH1'),:);
cfg.atlas = atlas;
cfg.inputcoord = 'mni';
cfg.output = 'label';
labels = ft_volumelookup(cfg, atlas);
[~, indx] = max(labels.count);
labels.name(indx)

%% Load and preprocess freq data
% This data has been rereferenced and preprocessed and converted to powerspectra
load([SBJ_dir SBJ '_freq.mat']);

% baseline correction
cfg = [];
cfg.baseline = [-.3 -.1];
cfg.baselinetype = 'relchange';
freq_blc = ft_freqbaseline(cfg, freq);

% Average for plotting
cfg = [];
cfg.frequency = [70 150];
cfg.avgoverfreq = 'yes';
cfg.latency = [0 0.8];
cfg.avgovertime = 'yes';
freq_sel = ft_selectdata(cfg, freq_blc);

%% Prepare and plot 2D layout
% amke layout
cfg = [];
cfg.headshape = pial_lh;
cfg.projection = 'orthographic';
cfg.channel = {'LPG*', 'LTG*'};
cfg.viewpoint = 'left';
cfg.mask = 'convex';
cfg.boxchannel = {'LTG30', 'LTG31'};
lay = ft_prepare_layout(cfg, freq);

% Plot interactive
cfg = [];
cfg.layout = lay;
cfg.showoutline = 'yes';
ft_multiplotTFR(cfg, freq_blc);

%% View grid activity on cortical mesh
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'surface';
cfg.interpmethod = 'sphere_weighteddistance';
cfg.sphereradius = 8;
cfg.camlight = 'no';
ft_sourceplot(cfg, freq_sel, pial_lh);
view([-90 20]); material dull; lighting gouraud; camlight;
% Add electrodes:
ft_plot_sens(elec_acpc_fr);

%% Plot SEEG data in 3D
% Create volumetric mask of ROIs from fs parcellation/segmentation
atlas = ft_read_atlas([SBJ_dir 'freesurfer/mri/aparc+aseg.mgz']);
atlas.coordsys = 'acpc';
cfg = [];
cfg.inputcoord = 'acpc';
cfg.atlas = atlas;
cfg.roi = {'Right-Hippocampus', 'Right-Amygdala'};
mask_rha = ft_volumelookup(cfg, atlas);

% Select electrodes of interest
cfg = [];
cfg.channel = {'RAM*', 'RTH*', 'RHH*'};
freq_sel2 = ft_selectdata(cfg, freq_sel);

% Plot HFA from bipolar channels via clouds around electrode positions
cfg = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim = [-.5 .5];
cfg.method = 'cloud';
cfg.slice = '3d';
cfg.nslices = 2;
cfg.facealpha = .25;
ft_sourceplot(cfg, freq_sel2, mesh_rha);
view([120 40]); lighting gouraud; camlight;

% 2D slice version:
cfg.slice = '2d';
ft_sourceplot(cfg, freq_sel2, mesh_rha);

% repeat this for different time points, and generate moveis!
