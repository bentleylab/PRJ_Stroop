function [data, header] = fn_cut_analysis_time_blocks(data, header, analysis_time)
%% Pull out sections of data to analyze, concatenate them together
%   Assumes analysis_time is in seconds

if(~isempty(analysis_time{1}))
    % Find correct data indices
    all_analysis_points = [];
    for time_grp = 1:length(analysis_time)
        all_analysis_points = [all_analysis_points ...
            (analysis_time{time_grp}(1)*header.sample_rate):(analysis_time{time_grp}(2)*header.sample_rate)];
    end
    % Check if any data indices are before 1
    if sum(all_analysis_points<1)
        disp('===================================================');
        fprintf('WARNING: %i ecog points less than 0\n',sum(all_analysis_points<1));
        disp('===================================================');
    end
    all_analysis_points(all_analysis_points<1) = [];
    % Check if any data indices are after the end of the recording
    if sum(all_analysis_points>header.n_samples)
        disp('===================================================');
        fprintf('WARNING: %i ecog points longer than n_samples\n',sum(all_analysis_points>header.n_samples));
        disp('===================================================');
    end
    all_analysis_points(all_analysis_points>header.n_samples) = [];
    
    % Select the appropriate data sections
    data = data(:,all_analysis_points);
    header.n_samples = size(data,2);
    header.length_in_seconds = header.n_samples / header.sample_rate;
end
header.analysis_time = analysis_time;

end

