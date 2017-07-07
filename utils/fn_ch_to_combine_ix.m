function avg_ch_ix = fn_ch_to_combine_ix(data_id)
% Returns which channels in a given dataset should be averaged for a grand average effect
% Inputs:
%   data_id [str] - the dataset that will be averaged (SBJ_ROI_ref)
% Outputs:
%   avg_ch_ix [int array] - indices of channels in dataset to average]
%       returns empty array for all channels

switch data_id
    case 'IR21_RC_WM'
        avg_ch_ix = [1:3]; %RC1-3, exclude RC4 for looking different
    case 'IR35_LAC_WM'
        avg_ch_ix = [1:3]; %LAC1-3, exclude LAC4 for noise
    case 'IR39_RAC_WM'
        avg_ch_ix = [];
    otherwise
        error(strcat('Channels to average not defined for data_id: ',data_id));
end
    
end