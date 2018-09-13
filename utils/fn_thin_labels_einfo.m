function [x_lab,y_lab] = fn_thin_labels_einfo(einfo,sort_vec,thin_axes)
% Removes repeated labels.
%      einfo = cell array with electrod info:
%      label- name of electrode
%      ROI- specific region
%      gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%      ROI2- specific region of second electrode
%      tissue- primary tissue type
%      GM weight- percentage of electrode pair in GM
%      Out- 0/1 flag for whether this is partially out of the brain

%   x thins by 1 (immediately repeated);
%   y thins so that only one label is shown for each patt_ROI label
%   assumes einfo{:,6} is xLabels and einfo{:,7} is yLabels
%   and that einfo{:,2} is ROI and einfo{:,3} is Patt


last_lab_x='placeHolder';
last_lab_y='placeHolder';
x_lab = cell([size(einfo,1) 1]);
y_lab = cell([size(einfo,1) 1]);
for e=1:size(einfo,1)
    cur_lab_x =char(einfo{e,sort_vec(1)});
    if thin_axes(1)~=0 && strcmp(cur_lab_x,last_lab_x) %mod(e,thin_axes(1))  % thins based on factor in thin_axes
        cur_lab_x=' ';
    else
        last_lab_x=char(einfo{e,sort_vec(1)});
    end
    x_lab{e}=cur_lab_x;
    
    cur_lab_y =char(einfo{e,sort_vec(2)});
    if thin_axes(2)~=0 && strcmp(cur_lab_y,last_lab_y)
        cur_lab_y=' ';
    else
        last_lab_y=char(einfo{e,sort_vec(2)});
    end
    y_lab{e}=cur_lab_y;
end
end
