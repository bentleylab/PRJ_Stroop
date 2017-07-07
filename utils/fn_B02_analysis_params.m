function [analysis_id, env_flag, filt_lim] = fn_B02_analysis_params(analysis, HG_type)
%% Returns processing parameters for a given analysis type
% analysis: [str] 'HG', 'theta', 'beta', 'alpha', 'ERP'
% HG_type: [str] 'multiband', 'wideband'

switch analysis
    case 'HG'
        env_flag = 1;
        env_id = '_env';
        if strcmp(HG_type,'wideband')
            analysis_id = 'HGw';
            filt_cf = [110];
            filt_bw = 80;
        elseif strcmp(HG_type,'multiband')
            analysis_id = 'HGm';
            filt_cf = [70:10:150];
            filt_bw = 10;
        else
            error(strcmp('Unknown HG filter type:',HG_type));
        end
        %     y_bound_chunk = 0.02;       % extend ylim by chunks of this size
    case 'theta'
        env_flag = 1;
        env_id = '_env';
        analysis_id = 'theta';
        filt_cf = [6];
        filt_bw = 4;
        %     y_bound_chunk = 0.5;
    case 'alpha'
        env_flag = 1;
        env_id = '_env';
        analysis_id = 'alpha';
        filt_cf = [10];
        filt_bw = 4;
        %     y_bound_chunk = 0.5;
    case 'beta'
        env_flag = 1;
        env_id = '_env';
        analysis_id = 'beta';
        filt_cf = [23];
        filt_bw = 8;
        %     y_bound_chunk = 0.5;
    case 'ERP'
        env_flag = 0;
        env_id = '';
        analysis_id = 'ERP';
        filt_lim = [];
        filt_cf = [];
        %     y_bound_chunk = 5;
    otherwise
        error(strcat('Unknown analysis selection: ',analysis));
end

% Convert center frequencies and bandwidths to high/low pass limits
for f_ix = 1:length(filt_cf)
    filt_lim{f_ix} = [filt_cf(f_ix)-filt_bw/2 filt_cf(f_ix)+filt_bw/2];
end

end
