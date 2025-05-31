function contribution = process_clusters(thresh_mask, data, h, E, H, dh)
    % Calculate TFCE contribution for a given threshold and cluster mask
    %
    % Inputs:
    %   thresh_mask: Logical mask of points exceeding "threshold"
    %   data: Original data array
    %   h: Current threshold
    %   E: Extent exponent (suggested value: 2/3; Mensen 2013, Neuroimage)
    %   H: Height exponent (suggested value: 2; Mensen 2013, Neuroimage)
    %   dh: Step size for integration
    %
    % Output:
    %   contribution: Array of TFCE contributions for the current threshold
    
    contribution = zeros(size(data));
    clusters = bwlabel(thresh_mask);
    unique_clusters = unique(clusters(clusters > 0));
    
    if size(unique_clusters,1) ~= 1
        unique_clusters = unique_clusters';
    end

    for cl = unique_clusters
        % Calculate cluster extent (number of points)
        cluster_extent = sum(clusters == cl);
        
        % Add TFCE contribution for the cluster
        contribution(clusters == cl) = (cluster_extent^E) * (h^H) * dh; %Riemann sum implementation of the integral
    end
end
