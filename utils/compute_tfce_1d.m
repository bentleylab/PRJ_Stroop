function TFCE_scores = compute_tfce_1d(data, dh, max_val, E, H)
    % Compute TFCE scores for 1D data (time series) with dynamic step size.
    %
    % Inputs:
    %   data: 1D array of t-values or other test statistic
    %   dh: step size for integration via Riemann Sum
    %   E: Extent exponent (recommended: 2/3 from Mensen 2013 NeuroImage)
    %   H: Height exponent (recommended: 2 from Mensen 2013 NeuroImage)
    %
    % Output:
    %   TFCE_scores: 1D array of TFCE scores
    
    % Initialize TFCE scores
    TFCE_scores = zeros(size(data));
    
    % Get vector of "heights"
    h_vals = 0:dh:max_val;

    % Iterate over all thresholds
    for h = h_vals
        % Find clusters above the current threshold for positive values
        pos_thresh = data >= h;  % Positive clusters
        neg_thresh = -data >= h; % Negative clusters
        
        % Add TFCE contributions for positive clusters
        TFCE_scores = TFCE_scores + process_clusters(pos_thresh, data, h, E, H, dh);
        
        % Add TFCE contributions for negative clusters
        TFCE_scores = TFCE_scores + -1 * process_clusters(neg_thresh, data, h, E, H, dh);
    end
end