function fn_generate_full_ROI_lists(atlas_id)
%% Create full_roi_lists.mat based on an atlas_id
%   These lists are used to cover all possibel ROIs/gROIs

if ~strcmp(atlas_id,'Dx')
    error('Why use anything but Dx here?');
end

[root_dir, ~] = fn_get_root_dir();

% Load atlas ROI lists
tsv_fname = [root_dir 'PRJ_Stroop/data/atlases/atlas_mappings/atlas_ROI_mappings_' atlas_id '.tsv'];
roi_file = fopen(tsv_fname, 'r');
% roi.csv contents:
%   atlas_label, gROI_label, ROI_label, Tissue, Known Variablility, Alternative ROIs, Notes (not read in)
roi_map = textscan(roi_file, '%s %s %s %s %d %s', 'HeaderLines', 1,...
    'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
fclose(roi_file);

all_grois = unique(roi_map{2});
all_rois  = unique(roi_map{3});
% Add in ROIs that don't have specific Dx label region (need manual adjustments)
all_grois = sort([all_grois; {'TWM'}]);
all_rois = sort([all_rois; {'SMC';'dmPFC';'TWM';'FO'}]);

% Save new lists
out_fname = [root_dir 'PRJ_Stroop/data/atlases/full_roi_lists_' atlas_id '.mat'];
save(out_fname,'-v7.3','all_rois','all_grois');

end