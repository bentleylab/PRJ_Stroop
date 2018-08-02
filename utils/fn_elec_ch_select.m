function [selected_elec] = fn_elec_ch_select(elec,ch_lab)
%% My version of ft_selectdata for .elec structs
% INPUTS:
%   elec [struct] - fieldtrip elec structure
%   ch_lab [str array] - cell array of desired channel labels (strings)

% Index of desired channels
% good_labels = ft_channelselection(ch_lab, elec.label);
[~,good_idx] = ismember(ch_lab,elec.label);

% Filter to desired channels
selected_elec = elec;
fields = fieldnames(elec);
for f = 1:numel(fields)
    if size(eval(['elec.' fields{f}]),1) == numel(elec.label)
        eval(['selected_elec.' fields{f} ' = elec.' fields{f} '(good_idx(good_idx~=0),:);']);
    end
end

end