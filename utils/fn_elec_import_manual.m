function fn_elec_import_manual(SBJ, proc_id, view_space, reg_type, atlas_id)
%% Import TSV file with manula adjustments of elec files
%   Added field:
%   elec.man_adj [0/1] - was this adjusted manually
%   elec.notes [str] - alternative ROIs usually

%% Load SBJ parameters and elec
% Check which root directory
[root_dir, ~] = fn_get_root_dir();
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);

if strcmp(reg_type,'v') || strcmp(reg_type,'s')
    reg_suffix = ['_' reg_type];
else
    reg_suffix = '';
end

elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' reg_suffix atlas_id '_compiled.mat'];
load(elec_fname);
save_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' reg_suffix atlas_id '_man.mat'];

%% Load manually adjusted TSV
tsv_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space '_' reg_suffix atlas_id '_man.tsv'];
tsv_file = fopen(tsv_fname, 'r');
headers = textscan(tsv_file, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 1);%, 'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
% Export order:
%   label, atlas_lab, atlas_prob, atlas_lab2, atlas_qryrng, ...
%   tissue, tissue_prob (4x), gm_weight, hemi, gROI, ROI, roi_flag,
%   man_adj, anatomy notes
man_adj = textscan(tsv_file, '%s %s %f %s %d %s %f %f %f %f %f %s %s %s %d %d %s', 'HeaderLines', 1,...
    'Delimiter', '\t', 'MultipleDelimsAsOne', 0);
fclose(tsv_file);

%% Add data back to elec
new_fields = {'tissue','hemi','gROI','ROI','man_adj','anat_notes'};
for f = 1:numel(new_fields)
    man_adj_ix = find(strcmp([headers{:}], new_fields{f}));
    if strcmp(new_fields{f},'man_adj')
        elec.man_adj = man_adj{man_adj_ix};
    elseif strcmp(new_fields{f},'anat_notes')
        elec.anat_notes = man_adj{man_adj_ix};
    else
        elec.(new_fields{f}) = man_adj{man_adj_ix};
    end
end

%% Save new elec
save(save_fname,'-v7.3','elec');

end

