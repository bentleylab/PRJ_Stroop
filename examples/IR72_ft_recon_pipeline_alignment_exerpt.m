%function [warped] = fn_warp_brains()
%% Atlas normalization from ft pipeline
% Below code copied from IR72 recon pipeline script

% Normalization
cfg = [];
cfg.nonlinear = 'yes';
cfg.template = '/home/knight/smg211/fieldtrip/template/anatomy/single_subj_T1_1mm.nii';
cfg.spmversion='spm12';
colin27 = ft_read_mri(cfg.template);
if exist([recondir 'Scans/' patient_id '_fsMR_pre_acpc.nii']) == 2
  fsmri_pre_acpc.coordsys = 'acpc';
  fsmri_mni_v = ft_volumenormalise(cfg, fsmri_pre_acpc);
elseif exist([recondir 'Scans/' patient_id '_fsMR_post_acpc.nii']) == 2
  fsmri_post_acpc.coordsys = 'acpc';
  fsmri_mni_v = ft_volumenormalise(cfg, fsmri_post_acpc);
end

% Warp Electrode Positions to Normalized Scan
if exist([recondir 'Electrodes/' patient_id '_elec_acpc_r.mat']) == 2
  elec_mni_frv = elec_acpc_r;
  elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni_v.params, elec_acpc_r.elecpos, 'individual2sn');
elseif exist([recondir 'Electrodes/' patient_id '_elec_acpc_f.mat']) == 2
  elec_mni_frv = elec_acpc_f;
  elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni_v.params, elec_acpc_f.elecpos, 'individual2sn');
else
  elec_mni_frv = elec_acpc_f;
  elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni_v.params, elec_acpc_f.elecpos, 'individual2sn');
end
elec_mni_frv.chanpos = elec_mni_frv.elecpos;
elec_mni_frv.coordsys = 'mni';


% %% 14b. Snap the Electrodes to the Template Brain (grid cases only)
% % Projection via Dykstra method
% ngrids = numel(elec_acpc_f.surface);
% template_brain = load('/home/knight/ecog_tmp/Recon_Functions/template_hull.mat');
% for g = 1:ngrids
%   cfg = [];
%   cfg.channel = elec_acpc_f.surface{g};
%   cfg.elec = elec_mni_frv;
%   cfg.casesensitive = 'yes';
%   cfg.method = 'headshape';
%   if strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.label, cfg.channel))), 'right');
%     cfg.headshape = template_brain.hull_rh;
%   elseif strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.label, cfg.channel))), 'left');
%     cfg.headshape = template_brain.hull_lh;
%   elseif strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.label, cfg.channel))), 'left-interhemispheric');
%     cfg.headshape = template_brain.hull_lh;
%   end
%   cfg.warp = 'hermes2010'; % hermes2010 or dykstra2012
%   cfg.pairmethod = 'label'; % pos or label
%   cfg.deformweight = 0.2;
%   ecog_mni{g} = ft_electroderealign(cfg);
%   ecog_mni{g}.unit = 'mm';
% end

%% 14c. View Result (do not edit)
figure; hold on;
surface_template_l = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_left.mat');
surface_template_r = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_right.mat');
if any(strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, {'ieeg_strip', 'ieeg_grid'}))), 'right')) && ...
    any(strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, {'ieeg_strip', 'ieeg_grid'}))), 'left'));
  % electrodes on both hemispheres
  ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
  ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
elseif strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, {'ieeg_strip', 'ieeg_grid'}))), 'left'); % electrodes only on left cortex
  ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
  view([-90 0]);
elseif strcmp(unique(elec_acpc_f.chanside(match_str(elec_acpc_f.chantype, {'ieeg_strip', 'ieeg_grid'}))), 'right'); % electrodes only on right cortex
  ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
  view([90 0]);
end
lighting gouraud; camlight; material dull;
ft_plot_sens(elec_mni_frv, 'elecshape', 'sphere', 'label', 'label');

%% 14d. Save Normalied MRI and Elec File
% Save Normalized MRI
cfg = [];
cfg.filename = [recondir 'Scans/' patient_id  '_fsMR_mni_v'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, fsmri_mni_v);

% Save Elec File
elec_mni_frv.surface = elec_acpc_f.surface;
elec_mni_frv.depths = elec_acpc_f.depths;
elec_mni_frv.grid_dim = elec_acpc_f.grid_dim;
save([recondir 'Electrodes/' patient_id '_elec_mni_frv.mat'], 'elec_mni_frv');

%% 15. Surface Based Normalization (only for grid cases; do not edit)
cfg = [];
cfg.detchantype = 'auto';
if exist([recondir 'Electrodes/' patient_id '_elec_acpc_r.mat']) == 2
  if ~isfield(elec_acpc_r, 'surface') || ~isfield(elec_acpc_r, 'chanside')
    elec = describe_elec(cfg, elec_acpc_r);
  else
    elec = elec_acpc_r;
  end
else
  if ~isfield(elec_acpc_f, 'surface')
    elec = describe_elec(cfg, elec_acpc_f);
  else
    elec = elec_acpc_f;
  end
end

% Make a structure of only surface electrodes
cfg = [];
cfg.channel = elec.label(match_str(elec.chantype, {'ieeg_strip', 'ieeg_grid'}));
elec_surf = select_elec(cfg, elec);

cfg = [];
cfg.elec = elec_surf;
cfg.method = 'headshape';
cfg.warp = 'fsaverage';
cfg.fshome = fshome;
% If there are electrodes on the left and right, they need to be normalized
% in two separate steps
if any(strcmp(unique(elec_surf.chanside(:)), 'right')) && any(strcmp(unique(elec_surf.chanside(:)), 'left'));
  elec_fsavg_frs_tmp = {};
  cfg.channel = elec_surf.label(match_str(elec_surf.chanside, {'right'}));
  cfg.headshape = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/rh.pial'];
  elec_fsavg_frs_tmp{1} = ft_electroderealign(cfg);

  cfg.channel = elec_surf.label(match_str(elec_surf.chanside, {'left'}));
  cfg.headshape = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/lh.pial'];
  elec_fsavg_frs_tmp{2} = ft_electroderealign(cfg);

  elec_fsavg_frs = ft_appendsens([], elec_fsavg_frs_tmp{:});
elseif strcmp(unique(elec_surf.chanside(:)), 'right') % electrodes only on the right
  cfg.channel = elec_surf.label;
  cfg.headshape = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/rh.pial'];
  elec_fsavg_frs = ft_electroderealign(cfg);
elseif strcmp(unique(elec_surf.chanside(:)), 'left') % electrodes only on the left
  cfg.channel = elec_surf.label;
  cfg.headshape = [subdir '3D_Images/freesurfer' scantype '/freesurfer/surf/lh.pial'];
  elec_fsavg_frs = ft_electroderealign(cfg);
end

% View Result
figure; hold on;
fs_surface_template_l = ft_read_headshape([fshome '/subjects/fsaverage/surf/lh.pial']);
fs_surface_template_l.coordsys = 'tal';
fs_surface_template_r = ft_read_headshape([fshome '/subjects/fsaverage/surf/rh.pial']);
fs_surface_template_r.coordsys = 'tal';
if any(strcmp(unique(elec_surf.chanside(:)), 'left')); % electrodes only on left cortex
  ft_plot_mesh(fs_surface_template_l, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'label','label');
  view([-90 0]);
end
if any(strcmp(unique(elec_surf.chanside(:)), 'right')); % electrodes only on right cortex
  ft_plot_mesh(fs_surface_template_r, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none', 'label', 'label');
  view([90 0]);
end
lighting gouraud; camlight; material dull;
hs = ft_plot_sens(elec_fsavg_frs, 'elecshape', 'sphere', 'label', 'label');

% Save Surface-Based Normalized Electrode Positions
save([recondir 'Electrodes/' patient_id '_elec_fsavg_frs.mat'], 'elec_fsavg_frs');

