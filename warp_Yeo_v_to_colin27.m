%% Atlas alignment
% Yeo atlases need to be warped to the Colin27 brain (all patients warped there)

%% Set up directories
[root_dir, ft_dir] = fn_get_root_dir();
yeo_v_dir = [root_dir 'PRJ_Stroop/data/atlases/Yeo/Yeo_JNeurophysiol11_MNI152/'];

%% Load Yeo's template brain
yeo_brain = ft_read_mri([yeo_v_dir 'FSL_MNI152_FreeSurferConformed_1mm.nii.gz']);
yeo_brain.coordsys = 'mni';

%% Warp Yeo template to FT template (Colin27)
cfg = [];
cfg.nonlinear  = 'yes';
cfg.template   = [ft_dir 'template/anatomy/single_subj_T1_1mm.nii'];
cfg.spmversion = 'spm12';
yeo_colin = ft_volumenormalise(cfg, yeo_brain);

cfg = [];
cfg.filename = [yeo_v_dir 'FSL_MNI152_FreeSurferConformed_1mm_colin27'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, yeo_colin);


%% Apply brain warps to atlas masks

% elec_mni_frv = elec_acpc_r;
% elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni_v.params, elec_acpc_r.elecpos, 'individual2sn');

%% Load atlases
yeo7  = ft_read_mri([yeo_v_dir 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz']);
yeo7.coordsys = 'mni';
yeo17 = ft_read_mri([yeo_v_dir 'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz']);
yeo17.coordsys = 'mni';

% NOTE: Need to use ft_warp_apply, rather than trying to warp the atlases
% straight to a brain...
% cfg = [];
% cfg.nonlinear  = 'yes';
% cfg.template   = [ft_dir 'template/anatomy/single_subj_T1_1mm.nii'];
% cfg.spmversion = 'spm12';
% cfg.write      = 'yes';
% yeo7_colin = ft_volumenormalise(cfg, yeo7);
% yeo17_colin = ft_volumenormalise(cfg, yeo17);

% Save Warped MRIs
cfg = [];
cfg.filename = [yeo_v_dir 'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, yeo7_colin);
cfg.filename = [yeo_v_dir 'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27'];
ft_volumewrite(cfg, yeo17_colin);


