function fn_elec_export_csv(SBJ, proc_id, view_space, reg_type, atlas_id)
%% Export elec struct to CSV
%   Adds gROI, ROI

%% Load elec and get ROI labels
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

%% Print data to csv file
csv_fname = [elec_fname(1:end-4) '.csv'];
csv = fopen(csv_fname,'w');

% Export order:
%   label, atlas_lab, atlas_prob, atlas_lab2, atlas_prob2, atlas_qryrng, ...
%   tissue, tissue_prob (4x), gm_weight, hemi, gROI, ROI, roi_flag
for e = 1:numel(elec.label)
    if isempty(elec.atlas_lab2{e})
        atlas_lab2 = '';
    else
        atlas_lab2 = strjoin(elec.atlas_lab2{e},'+');
    end
    fprintf(csv,'%s,%s,%.3f,%s,%.3f,%d,%s,%.3f,%.3f,%.3f,%.3f,%.1f,%s,%s,%s,%d\n',...
        elec.label{e},elec.atlas_lab{e},elec.atlas_prob(e),atlas_lab2,elec.atlas_prob2{e}(1),elec.atlas_qryrng(e),...
        elec.tissue{e},elec.tissue_prob(e,:),elec.gm_weight(e),elec.hemi{e},...
        elec.gROI{e},elec.ROI{e},elec.roi_flag(e));
end

fclose(csv);

end

