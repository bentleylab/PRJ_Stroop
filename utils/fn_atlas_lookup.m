function [roi_labels] = fn_atlas_lookup(elec, atlas_name)
% generate a list of roi_labels for a given elec struct and atlas
% INPUTS:
%   elec [struct] - elec struct
%   atlas_name [str] - {'fs',...?}

%% Find ROI Labels of electrodes (Create elec table)
% edit generate_electable.m % more extensive version
if strcmp(atlas_name,'aparc')
    atlas      = ft_read_atlas([fsdir ...
        '/mri/aparc+aseg.mgz']); % Desikan-Killiany (+volumetric)
    atlas.coordsys = 'mni';
    name       = 'Desikan-Killiany';
elseif strcmp(atlas_name,'aparc.a2009s')
    atlas      = ft_read_atlas([fsdir ...
        '/mri/aparc.a2009s+aseg.mgz']); % Destrieux (+volumetric)
    atlas.coordsys = 'mni';
    name       = 'Destrieux';
else
    error(['atlas_name unknown: ' atlas_name]);
end
% elec.elecpos_fs   = elec.elecpos;
fprintf('Using atlas: %s\n',name);

% Find ROI labels for chanpos
cfg = [];
cfg.roi = elec.chanpos;
cfg.atlas = atlas;
cfg.inputcoord = elec.coordsys;
cfg.output = 'label';
roi_labels = ft_volumelookup(cfg, atlas);

% % Import atlas of interest
% atlas = ft_read_atlas([ft_dir 'template/atlas/aal/ROI_MNI_V4.nii']);
% [~, indx] = max(labels.count);
% labels.name(indx)


end


