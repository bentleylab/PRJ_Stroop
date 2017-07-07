function orig_times = fn_convert_evnt_times_at2full(times,analysis_time,sample_rate)
%% Convert time markers from trimmed analysis_time to original time
%   times [int array] - time points (samples) in time series already cut by analysis_time

% Convert analysis_time to samples
for ix = 1:numel(analysis_time)
    analysis_samples{ix} = analysis_time{ix}.*sample_rate;
    analysis_seg_len{ix} = diff(analysis_samples{ix});
end

% Add segment length back to times
orig_times = nan(size(times));
seg_ix = 1;
seg_end = analysis_samples{seg_ix}(2);
seg_start = analysis_samples{seg_ix}(1);
for t_ix = 1:length(times)
    % Compute length of previous segments to offset and avoid counting those twice
    if seg_ix==1
        offset = 0;
    else
        offset = sum(analysis_seg_len{1:seg_ix-1});
    end
    
    if times(t_ix)+seg_start-offset > seg_end
        seg_ix = seg_ix+1;
        seg_end = analysis_samples{seg_ix}(2);
        seg_start = analysis_samples{seg_ix}(1);
    end
    
    % Add back the start of the first segment
    orig_times(t_ix) = times(t_ix)+seg_start-offset;
end

end