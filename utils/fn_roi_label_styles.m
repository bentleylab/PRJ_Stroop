function [labels, colors] = fn_roi_label_styles(roi_id)
%% Converts the name of a set of ROIs into labels, plotting colors
% condition_name: [str] 'ROI', 'gROI', 'INS', 'LPFC', 'MPFC', 'OFC'
% colors from http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=3
%
% einfo_col = 2 for specific ROIs, 3 for general ROIs

[root_dir, ~] = fn_get_root_dir();

% if length(cond_lab) == 1
switch roi_id
    case 'ROI'
        load([root_dir 'PRJ_Stroop/data/atlases/full_roi_lists_Dx.mat']);
        labels = all_rois;
        % Exclude FWM, '', OUT
        labels(strmatch('FWM',labels,'exact')) = [];
        labels(strmatch('TWM',labels,'exact')) = [];
        labels(strmatch('OUT',labels,'exact')) = [];
        labels(strmatch('',labels,'exact')) = [];
    case 'Yeo7'
        labels = {'Vis','SM','DAttn','VAttn','Limb','FP','Def'};
    case 'gPFC'
        labels = {'LPFC','MPFC'};
    case 'PFC'
        labels = {'FPC','dlPFC','vlPFC','PM','FO','ACC','dmPFC','SMC','aMCC','pMCC'};
    case 'Main3'
        labels = {'LPFC','MPFC','INS'};
    case 'mgROI'
        labels = {'LPFC','MPFC','SM','INS','OFC'};
    case 'gROI'
        labels = {'LPFC','MPFC','INS','SM','OFC','TMP'};%'PAR','AMG','HPC','OCC'};
    case 'all_gROI'
        labels = {'LPFC','MPFC','INS','SM','OFC','PAR','TMP','AMG','HPC','OCC','FWM','TWM'};
    case 'lat'
        labels = {'LPFC','SM','TMP'};%'PAR','OCC'
    case 'deep'
        labels = {'INS'};%'AMG','HPC'
    case 'mnLPFC'
        labels = {'dlPFC','vlPFC','PM','aMCC','SMC'};
    case 'thryROI'
        labels = {'dlPFC','vlPFC','PM','aMCC','SMC','aINS','FO'};
    case 'SM'
        labels = {'M1','S1'};
    case 'PAR'
        labels = {'SPL','IPL','PCC','PRC'};
    case 'TMP'
        labels = {'PT','STG','STS','ITC','VTC'};
    case 'LPFC'
        labels = {'FPC','dlPFC','vlPFC','PM','FO'};
    case 'MPFC'
        labels = {'ACC','dmPFC','SMC','aMCC','pMCC'};
    case 'INS'
        labels = {'aINS','pINS'};
    case 'OFC'
        labels = {'mOFC','lOFC'};
    case 'MTL'
        labels = {'AMG','HPC'};
    case {'tissue', 'tissueC'}
        labels = {'GM','WM','CSF','OUT'};
    case 'all'
        load([root_dir 'PRJ_Stroop/data/atlases/full_roi_lists_Dx.mat']);
        labels = all_rois;
    otherwise
        error(strcat('Unknown roi_id: ',roi_id));
end

% Get colors
colors = cell(size(labels));
for roi_ix = 1:numel(labels)
    colors{roi_ix} = fn_roi2color(labels{roi_ix});
end

end
