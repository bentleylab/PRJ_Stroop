function [roi_labels] = fn_atlas_lookup(elec, atlas_name)


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


end