function combinations = fn_create_channel_combinations(MPFC, LPFC)
    % Get the number of elements in MPFC and LPFC
    numMPFC = length(MPFC);
    numLPFC = length(LPFC);

    % Calculate the total number of combinations
    totalCombinations = numMPFC * numLPFC;

    % Initialize the output cell array
    combinations = cell(totalCombinations, 2);

    % Fill the combinations cell array
    index = 1;
    for i = 1:numMPFC
        for j = 1:numLPFC
            combinations(index, 1) = MPFC(i);
            combinations(index, 2) = LPFC(j);
            index = index + 1;
        end
    end
end
