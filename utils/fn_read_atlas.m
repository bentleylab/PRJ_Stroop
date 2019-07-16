function [atlas] = fn_read_atlas(atlas_id)
%% Custom version of ft_read_atlas to work with Yeo atlas
%   modeled after lines 1722:1746 from ft_read_atlas

[root_dir, ~] = fn_get_root_dir();

%% Load Atlas
% Select atlas filename
switch atlas_id
    case 'Yeo7'
        atlas_fname = [root_dir 'PRJ_Stroop/data/atlases/Yeo/Yeo_JNeurophysiol11_MNI152/'...
            'Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'];
        atlas_table_fname = [root_dir 'PRJ_Stroop/data/atlases/Yeo/Yeo_JNeurophysiol11_MNI152/'...
            'Yeo2011_7Networks_ColorLUT.txt'];
        value = 1:7;
        label = {
            'Visual'%'7Networks_1'
            'Somatomotor'%'7Networks_2'
            'Dorsal Attention'%'7Networks_3'
            'Ventral Attention'%7Networks_4'
            'Limbic'%7Networks_5'
            'Frontoparietal'%7Networks_6'
            'Default'%7Networks_7'
            };
        colors = [
            120  18 134;
            70 130 180;
            0 118  14;
            196  58 250;
            220 248 164;
            230 148  34;
            205  62  78
            ];
    case 'Yeo17'
        atlas_fname = [root_dir 'PRJ_Stroop/data/atlases/Yeo/Yeo_JNeurophysiol11_MNI152/'...
            'Yeo2011_17Networks_MNI152_FreeSurferConformed1mm_LiberalMask_colin27.nii'];
        atlas_table_fname = [root_dir 'PRJ_Stroop/data/atlases/Yeo/Yeo_JNeurophysiol11_MNI152/'...
            'Yeo2011_7Networks_ColorLUT.txt'];
        value = 1:17;
        label = {
            'Visual A'%'17Networks_1'
            'Visual B'%'17Networks_2'
            'Somatomotor A'%'17Networks_3'
            'Somatomotor B'%'17Networks_4'
            'Temporal Parietal'%'17Networks_5'
            'Dorsal Attention A'%'17Networks_6'
            'Dorsal Attention B'%'17Networks_7'
            'Salience A'%'17Networks_8'
            'Salience B'%'17Networks_9'
            'Control A'%'17Networks_10'
            'Control B'%17Networks_11'
            'Control C'%17Networks_12'
            'Default A'%'17Networks_13'
            'Default B'%'17Networks_14'
            'Default C'%'17Networks_15'
            'Limbic A'%'17Networks_16'
            'Limbic B'%'17Networks_17'
            };
        colors = [
            120  18 134;
            255   0   0;
            70 130 180;
            42 204 164;
            74 155  60;
            0 118  14;
            196  58 250;
            255 152 213;
            220 248 164;
            122 135  50;
            119 140 176;
            230 148  34;
            135  50  74;
            12  48 255;
            0   0 130;
            255 255   0;
            205  62  78
            ];
end

% Read in the volume
atlas = ft_read_mri(atlas_fname);
dat   = atlas.anatomy;
atlas = rmfield(atlas, 'anatomy');

% ft's read_annotation 
% % Load Labels and Mapping
% [v, p, c] = read_annotation(atlas_table_fname);
% 
% label = c.struct_names;
% rgba  = c.table(:,1:4);
% rgb   = c.table(:,5); % compound value that is used for the indexing in vector p
% index = ((1:c.numEntries)-1)';


%% get the unique values to save time later on
uval = unique(dat(:));
sel  = find(ismember(value, uval));
fprintf('subselecting %d labels from the total list of %d\n', numel(sel), numel(label));
value = value(sel);
label = label(sel);

% remap the values in the data
rois = zeros(size(dat));
cnt   = 0;
for k = 1:numel(value)
    sel = dat==value(k);
    if sum(sel(:))
        cnt = cnt+1;
        fprintf('re-indexing label %s to a value of %d (was %d)\n', label{k}, cnt, value(k));
        rois(sel)      = cnt;
        roilabel{cnt,1} = label{k};
    end
end
atlas.roi      = rois;
atlas.roilabel = roilabel;

end
