function csd_mat = fn_create_CSD_matrix(crsspctrm, labelcmb, MPFC, LPFC)
    % Get the dimensions of the input data
    [numCombinations, numFrequencies] = size(crsspctrm);
    numMPFC = length(MPFC);
    numLPFC = length(LPFC);

    % Initialize the output matrix
    csd_mat = zeros(numMPFC, numLPFC, numFrequencies);

    % Fill the result matrix using the label combinations
    for i = 1:numCombinations
        % Find the indices of the current MPFC and LPFC labels
        mpfcIndex = find(strcmp(MPFC, labelcmb{i, 1}));
        lpfcIndex = find(strcmp(LPFC, labelcmb{i, 2}));

        % Place the cross-spectrum data into the appropriate location
        csd_mat(mpfcIndex, lpfcIndex, :) = crsspctrm(i, :);
    end
end
