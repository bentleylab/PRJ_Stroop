function elec = loadrecon(ptid)

%% LOADRECON loads the reconstruction for a given patient. Electrodes can be
% viewed in patient space or normalized space, and as an orthoplot or in 3D.
% If the recon is unfinished or in progress,LOADRECON should tell you how
% far along the recon is
%
% Use as loadrecon(patient_ID) where patient_ID is a string
addpath(genpath('/home/knight/ecog_tmp/'));

if any(strfind(ptid(1:2), 'OS')); hosp_code = 'Oslo';
elseif any(strfind(ptid(1:2), 'IR')); hosp_code = 'Irvine';
elseif any(strfind(ptid(1:2), 'CH')); hosp_code = 'Childrens';
elseif any(strfind(ptid(1:2), 'AM')); hosp_code = 'Albany';
elseif any(strfind(ptid(1:2), 'WV')); hosp_code = 'Wadsworth';
elseif any(strfind(ptid(1:2), 'GP')); hosp_code = 'UCSF';
elseif any(strfind(ptid(1:2), 'SF')); hosp_code = 'SFC';
elseif any(strfind(ptid(1:2), 'CH')); hosp_code = 'Childrens';
elseif any(strfind(ptid(1:2), 'JH')); hosp_code = 'Hopkins';
elseif strfind(ptid(1:2), 'P' == 1); hosp_code = 'Germany';
elseif any(strfind(ptid(1:2), 'CP')); hosp_code = 'CPMC';
elseif any(strfind(ptid(1:2), 'NY')); hosp_code = 'New_York';
end
subdir = ['/home/knight/ecog/DATA_FOLDER/' hosp_code '/' ptid '/'];
ogrecondir = dir([subdir '3D_Images/Recon_*']);
if ~isempty(ogrecondir)
  recon_date = ogrecondir.name(end-7:end); % 3 letter month code + year
elseif ~isempty(dir([subdir '3D_Images/DICOM_Anon/CT*'])) && ~isempty(dir([subdir '3D_Images/DICOM_Anon/MR*']));
  error('No recon has been started yet for this subject, but we do have CT and MR scans so one should be started soon');
else
  error('No CT and MR scans were found');
end
recondir = [subdir '3D_Images/Recon_' recon_date '/FT_Pipeline/'];
if exist([recondir '/Electrodes/' ptid '_elec_acpc.mat'])
  load([recondir 'Electrodes/' ptid '_elec_acpc.mat'])
else
  error('The recon has been started, but is not yet finished')
end

method = input('What do you want to see? Enter "orthoplot" or "3d"\n', 's');
if strcmp(method, 'orthoplot')
  methodortho = input('Do you want to see the electrodes in patient space or normalized space?\nEnter "patient" or "normalized"\n', 's');
  if strcmp(methodortho, 'patient')
    ct_acpc_f = ft_read_mri([recondir 'Scans/' ptid '_CT_acpc_f.nii']);
    ct_acpc_f.coordsys = 'acpc';
    % fsmri_acpc = ft_read_mri(['Scans/' patient_id '_fsMR_pre_acpc.nii']);
    if exist([recondir 'Scans/' ptid '_fsMR_post_acpc.nii']) && exist([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']);
      scantype = input('There is a post-op and pre-op MR for this subject, which would you like to see?\nEnter "post" or "pre"\n', 's');
      if strcmp(scantype, 'post')
        fsmri_post_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_post_acpc.nii']);
        fsmri_post_acpc.coordsys = 'acpc';

        cfg = [];
        cfg.elec = elec_acpc;
        ft_electrodeplacement(cfg, ct_acpc_f, fsmri_post_acpc);
      elseif strcmp(scantype, 'pre')
        if exist([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
        load([recondir '/Electrodes/' ptid '_elec_acpc_f.mat'])
        else
          error('The post-op MRI has not been snapped to the pre-op MRI yet, so electrodes can only be viewed in the post-op MRI')
        end
        fsmri_pre_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']);
        fsmri_pre_acpc.coordsys = 'acpc';

        cfg = [];
        cfg.elec = elec_acpc;
        ft_electrodeplacement(cfg, ct_acpc_f, fsmri_pre_acpc);
      else
        error('You must enter "post" or "pre"')
      end
    elseif exist([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']) == 2
      fsmri_pre_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']);
      fsmri_pre_acpc.coordsys = 'acpc';
      
      cfg = [];
      cfg.elec = elec_acpc;
      ft_electrodeplacement(cfg, ct_acpc_f, fsmri_pre_acpc);
    elseif exist([recondir 'Scans/' ptid '_fsMR_post_acpc.nii']) == 2
      fsmri_post_acpc = ft_read_mri([recondir 'Scans/' ptid '_fsMR_post_acpc.nii']);
      fsmri_post_acpc.coordsys = 'acpc';
      
      cfg = [];
      cfg.elec = elec_acpc;
      ft_electrodeplacement(cfg, ct_acpc_f, fsmri_post_acpc);
    elseif exist([recondir 'Scans/' ptid '_fsMR_pre_acpc.nii']) == 2
      fsmri_pre_acpc = ft_read_mri(['Scans/' ptid '_fsMR_pre_acpc.nii']);
      fsmri_pre_acpc.coordsys = 'acpc';
      
      cfg = [];
      cfg.elec = elec_acpc;
      ft_electrodeplacement(cfg, ct_acpc_f, fsmri_pre_acpc);
    end
  elseif strcmp(methodortho, 'normalized')
    mri_template = ft_read_mri('/home/knight/smg211/fieldtrip/template/anatomy/single_subj_T1_1mm.nii');
    mri_template.coordsys = 'mni';
    if exist([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
      load([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
    else
      error('The recon has been started but the electrodes have not been volume-based normalized yet')
    end
    
    cfg = [];
    cfg.elec = elec_mni_v;
    ft_electrodeplacement(cfg, mri_template);
  else
    error('you must enter "patient" or "normalized"')
  end
elseif strcmp(method, '3d')
  figure;
  method3d = input('Do you want to see the electrodes in patient space or normalized space?\nEnter "patient" or "normalized"\n', 's');
  if strcmp(method3d, 'patient')
    load([recondir 'Surfaces/' ptid '_cortex_lh.mat']);
    load([recondir 'Surfaces/' ptid '_cortex_rh.mat']);
    ft_plot_mesh(cortex_rh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
    ft_plot_mesh(cortex_lh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
    view([-55 10]); lighting gouraud; camlight;
    
    if exist([recondir 'Electrodes/' ptid '_elec_acpc_r.mat'])
      load([recondir 'Electrodes/' ptid '_elec_acpc_r.mat'])
    else
      error('The recon has been started but the electrodes have not been projected to the surface yet or there are no surface electrodes')
    end
    ft_plot_sens(elec_acpc_r, 'elecshape', 'sphere');
    
  elseif strcmp(method3d, 'normalized')
    method3dnorm = input('Do you want to see the volume normalized-based or surface-based normalized electrodes?\nEnter "volume" or "surface"\n', 's');
    if strcmp(method3dnorm, 'volume')
      surface_template_l = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_left.mat');
      surface_template_r = load('/home/knight/smg211/fieldtrip/template/anatomy/surface_pial_right.mat');
      ft_plot_mesh(surface_template_l.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      ft_plot_mesh(surface_template_r.mesh, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      view([-55 10]); lighting gouraud; camlight;
      
      if exist([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
        load([recondir 'Electrodes/' ptid '_elec_mni_v.mat'])
      else
        error('The recon has been started but the electrodes have not been volume-based normalized yet')
      end
      ft_plot_sens(elec_mni_v, 'elecshape', 'sphere');
      
    elseif strcmp(method3dnorm, 'surface')
      fshome = '/usr/local/freesurfer_x86_64-5.3.0';
      fs_surface_template_l = ft_read_headshape([fshome '/subjects/fsaverage/surf/lh.pial']);
      fs_surface_template_l.coordsys = 'tal';
      fs_surface_template_r = ft_read_headshape([fshome '/subjects/fsaverage/surf/rh.pial']);
      fs_surface_template_r.coordsys = 'tal';
      ft_plot_mesh(fs_surface_template_l, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      ft_plot_mesh(fs_surface_template_r, 'facecolor', [0.781 0.762 0.664], 'EdgeColor', 'none');
      view([-55 10]); lighting gouraud; camlight;
      
      if exist([recondir 'Electrodes/' ptid '_elec_mni_s.mat'])
        load([recondir 'Electrodes/' ptid '_elec_mni_s.mat'])
      else
        error('The recon has been started but the electrodes have either not been surface-based normalized or there are only depth, in which case surface-based normalization is not possible')
      end
      ft_plot_sens(elec_mni_s, 'elecshape', 'sphere');
      
    else
      error('you must enter "volume" or "surface"')
    end
  else
    error('you must enter "patient" or "normalized"')
  end
else
  error('you must enter "orthoplot" or "3d"')
end
