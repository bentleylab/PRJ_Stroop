function [elec_lab] = fn_atlas_lookup(elec, atlas)
% generate a list of roi_labels for a given elec struct and atlas
% INPUTS:
%   elec [struct] - elec struct
%   atlas [struct] - atlas loaded from ft_read_atlas


% edit generate_electable.m % more extensive version

%% Find ROI Labels of electrodes (Create elec table)
% Initialize the new fields
elec_lab = elec;
elec_lab.atlas_name   = atlas.name;
elec_lab.atlas_label  = cell(size(elec.label));
elec_lab.atlas_count  = zeros(size(elec.label));
elec_lab.atlas_qryrng = cell(size(elec.label));
elec_lab.atlas_label2  = cell(size(elec.label));
elec_lab.atlas_count2  = cell(size(elec.label));
elec_lab.atlas_qryrng2 = cell(size(elec.label));

cfg = [];
cfg.atlas = atlas;
cfg.inputcoord = elec.coordsys;
cfg.output = 'single';
cfg.maxqueryrange = 5;
for e = 1:numel(elec.label);
    % Search for this elec
    cfg.roi = elec.chanpos(e,:);
    report = ft_volumelookup(cfg,atlas);
    
    % Report label for this elec
    [sorted_cnt cnt_idx] = sort(report.count,1,'descend');
    match_cnt = find(sorted_cnt);
    
    % Assign highest match as label
    if numel(match_cnt)>=1
        elec_lab.atlas_label{e}  = report.name{cnt_idx(1)};
        elec_lab.atlas_count{e}  = report.count(cnt_idx(1));
        elec_lab.atlas_qryrng{e} = report.usedqueryrange{cnt_idx(1)};
        % add additional labels if needed
        if numel(match_cnt)>1
            elec_lab.atlas_label2{e} = report.name(cnt_idx(2:numel(match_cnt)));
            elec_lab.atlas_count2{e} = report.count(cnt_idx(2:numel(match_cnt)));
            elec_lab.atlas_qryrng2{e} = report.usedqueryrange(cnt_idx(2:numel(match_cnt)));
        else
            elec_lab.atlas_label2{e} = '';
            elec_lab.atlas_count2{e} = NaN;
            elec_lab.atlas_qryrng2{e} = [];
        end
    else   % No matches found
        error(['No matches found for elec: ' elec.label{e}]);
    end
%     for j=1:length(match_cnt);
%         found{j} = report.name{cnt_idx(j)};%[num2str(report.count(ind(j))) ': ' report.name{ind(j)}];
%     end
%     roi_labels{e} = strjoin(found,',');
end

% nf_cnt = 0;
% for e = 1:numel(elec.label)
%     disp([num2str(e) ': ' elec.label{e} ' = ' roi_labels{e}]);
%     if ~isempty(strfind(roi_labels{e},'no_label_found'))
%         nf_cnt= nf_cnt+1;
%     end
% end
% disp(['nf_cnt = ' num2str(nf_cnt)]);
% % Import atlas of interest
% atlas = ft_read_atlas([ft_dir 'template/atlas/aal/ROI_MNI_V4.nii']);
% [~, indx] = max(labels.count);
% labels.name(indx)


end


