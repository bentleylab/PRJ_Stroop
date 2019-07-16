function [RGB] = fn_atlas2color(atlas_id,rois)
%% Returns the RGB color code to plot a given ROI
% INPUTS:
%   label [cell array of strs] - name (or list of names) of the ROI label
%       right now only Yeo7 and Yeo17
% OUTPUTS:
%   RGB [3 floats] - RGB for this ROI
%

switch atlas_id
    case 'Yeo7'
        labels = {
            'Vis'%'7Networks_1'
            'SM'%'7Networks_2'
            'DAttn'%'7Networks_3'
            'VAttn'%'7Networks_4'
            'Limb'%'7Networks_5'
            'FP'%'7Networks_6'
            'Def'%'7Networks_7'
            'OUT'
            };
        colors = [
            120  18 134;
            70 130 180;
            0 118  14;
            196  58 250;
            220 248 164;
            230 148  34;
            205  62  78;
            0 0 0
            ];
    case 'Yeo17'
        labels = {
            'Vis_A'%'17Networks_1'
            'Vis_B'%'17Networks_2'
            'SM_A'%'17Networks_3'
            'SM_B'%'17Networks_4'
            'TmpPar'%'17Networks_5'
            'DAttn_A'%'17Networks_6'
            'DAttn_B'%'17Networks_7'
            'Sal_A'%'17Networks_8'
            'Sal_B'%'17Networks_9'
            'Ctrl_A'%'17Networks_10'
            'Ctrl_B'%17Networks_11'
            'Ctrl_C'%17Networks_12'
            'Def_A'%'17Networks_13'
            'Def_B'%'17Networks_14'
            'Def_C'%'17Networks_15'
            'Limb_A'%'17Networks_16'
            'Limb_B'%'17Networks_17'
            'no_label_found'
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
            0 0 0
            ];
end

RGB = zeros([numel(rois) 3]);
for r = 1:numel(rois)
    RGB(r,:) = colors(strcmp(rois{r},labels),:)./256;
end

end
