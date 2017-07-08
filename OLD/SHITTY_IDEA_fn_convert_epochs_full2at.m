function cut_epochs = fn_convert_epochs_full2at(epochs,analysis_time,sample_rate,keep_partial)
%% Convert epochs from original time to trimmed analysis_time
%   Will set time points to NaN that are outside of analysis_time
% INPUTS:
%   times [int array] - time points (samples) in original time series (not yet cut)
%   analysis_time [cell array of tuples] - contains [start stop] times in SECONDS of valid segments of data
%   sample_rate [int] - sampling rate of the data (to convert time to samples)
%   keep_partial [0/1] - flag to handle cases where epochs stradle segment borders
%       0: leave time points outside of valid segments as NaN
%       1: adjust epoch edges to cover same window but only including valid points

error('do not use this function, it is way too complicated and has bugs and errors!!!');
% Convert analysis_time to samples
valid_samples = [];
for ix = 1:numel(analysis_time)
    analysis_samples{ix} = analysis_time{ix}.*sample_rate;
    analysis_seg_len{ix} = diff(analysis_samples{ix});
    valid_samples = [valid_samples analysis_samples{ix}(1):analysis_samples{ix}(2)];
end

% Initialize variables
cut_epochs = nan(size(epochs));
seg_ix    = 1;
seg_start = analysis_samples{seg_ix}(1);
seg_end   = analysis_samples{seg_ix}(2);
offset    = 0;
for t_ix = 1:size(epochs,1)
    %     if ~((epochs(t_ix,2) < seg_start) && (seg_ix==1)) %ignore epochs before the first segment
    % Check for overlap with this segment
    valid = [any(epochs(t_ix,1)==valid_samples) any(epochs(t_ix,2)==valid_samples)]
    
    % If overlapping at all with valid segment, adjust epoch
    if any(valid)
        if all(valid)   % Both are good, adjust them
            cut_epochs(t_ix,:) = epochs(t_ix,:)-seg_start;
        elseif keep_partial % Only one is good, adjust edge to fit in valid times
            if ~valid(1)
                cut_epochs(t_ix,1) = seg_start;
                cut_epochs(t_ix,2) = epochs(t_ix,2)-seg_start;
            else
                assert(~valid(2));
                cut_epochs(t_ix,1) = epochs(t_ix,1)-seg_start;
                cut_epochs(t_ix,2) = seg_start+analysis_seg_len{seg_ix};
            end
        end
    else
        % If before segment starts, do nothing
        if epochs(t_ix,2)<seg_start
            continue;
        elseif epochs(t_ix,1)>seg_end
            seg_ix    = seg_ix+1;
            % If completely past segment, advance segment or end loop
            if seg_ix>numel(analysis_samples)
                break;
            else
                seg_start = analysis_samples{seg_ix}(1);
                seg_end   = analysis_samples{seg_ix}(2);
                offset    = analysis_seg_len{seg_ix-1};
            end
        end
    end
    
end

%     % Subtract gap before this segment
%     cut_epochs(t_ix,:) = epochs(t_ix,:)-seg_start;
%     if cut_epochs(t_ix,1)<0
%         if keep_partial
%             cut_epochs(t_ix,1) = 1;
%             
%             switch
%                 case
%                     cut
%                     
%                     
%                     % Roll over segment limits if passed end of last segment
%                     if epochs(t_ix,2) > seg_end
%                     end
%                     %     end
%             end
%             
%             % Set times before the onset of the first segment to NaN
%             cut_epochs(
