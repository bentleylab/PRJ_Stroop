function fn_elec_export_csv(SBJ, proc_id, view_space, reg_type, atlas_id, reref)
%% Export elec struct to CSV for manual edits on google spreadsheet
%   If orig, adds gROI and ROI

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
if reref
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' view_space reg_suffix atlas_id '_compiled.mat'];
else
    elec_fname = [SBJ_vars.dirs.recon,SBJ,'_elec_' proc_id '_',view_space,reg_suffix,'_orig_',atlas_id,'.mat'];
end
load(elec_fname);

%% Add fields for orig elecs
if ~reref
    % ROI matching
    elec.gROI = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'gROI');
    elec.ROI  = fn_atlas2roi_labels(elec.atlas_lab,atlas_id,'ROI');
    
    % Fake reref fields
    elec.gm_weight = elec.tissue_prob(:,strcmp('GM',elec.tissue_labels));
    elec.roi_flag  = zeros(size(elec.label));
end

%% Print data to csv file
csv_fname = [elec_fname(1:end-4) '.csv'];
fprintf('\tExporting %s...\n',csv_fname);
csv = fopen(csv_fname,'w');

% Export order:
%   label, atlas_lab, atlas_prob, atlas_lab2, atlas_qryrng, ...
%   tissue, tissue_prob (4x), gm_weight, hemi, gROI, ROI, roi_flag
for e = 1:numel(elec.label)
    if ~isempty(elec.atlas_lab2{e})
        if ischar(elec.atlas_lab2{e})
            atlas_lab2 = elec.atlas_lab2{e};
        else
            atlas_lab2 = strjoin(elec.atlas_lab2{e},'+');
        end
    else
        atlas_lab2 = '';
    end
    fprintf(csv,'%s,%s,%.3f,%s,%d,%s,%.3f,%.3f,%.3f,%.3f,%.1f,%s,%s,%s,%d\n',...
        elec.label{e},elec.atlas_lab{e},elec.atlas_prob(e),atlas_lab2,elec.atlas_qryrng(e),...
        elec.tissue{e},elec.tissue_prob(e,:),elec.gm_weight(e),elec.hemi{e},...
        elec.gROI{e},elec.ROI{e},elec.roi_flag(e));
end

fclose(csv);

end

