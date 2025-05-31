
% Load  Data
load('PreparedData.mat'); 

% Define Parameters for the Statistics and TFCE
n_permutations = 1000;
nsteps = 100;
E = 0.66;
H = 2;

% Define the time vector (adjust according to your data)
nTimePoints = length(newt);

modelFormula = 'powerValues ~ CurrentConflict + PreviousConflict + PC + CurrentConflict:PreviousConflict + CurrentConflict:PC';


StatsResults.dlPFC = {};
for i = 1:numel(preparedData.dlPFC)
    currElec = preparedData.dlPFC{i};
    subjID = currElec.Subject;
    
    % Retrieve the subject-specific design matrix
    subjectDesignMatrix = preparedData.subjectDesignMatrices.(subjID);
    
    %% Update the design matrix as needed:
    subjectDesignMatrix.PC = subjectDesignMatrix.PC - 0.33;
    
    % Create CurrentConflict column
    subjectDesignMatrix.CurrentConflict = repmat({'NoConflict'}, height(subjectDesignMatrix), 1);
    subjectDesignMatrix.CurrentConflict(strcmp(subjectDesignMatrix.CurrentType, 'inc')) = {'Conflict'};
    
    % Create PreviousConflict column
    subjectDesignMatrix.PreviousConflict = repmat({'NoConflict'}, height(subjectDesignMatrix), 1);
    subjectDesignMatrix.PreviousConflict(strcmp(subjectDesignMatrix.PreviousType, 'inc')) = {'Conflict'};
    
    subjectDesignMatrix.CurrentConflict = categorical(subjectDesignMatrix.CurrentConflict);
    subjectDesignMatrix.CurrentConflict = reordercats(subjectDesignMatrix.CurrentConflict, {'NoConflict', 'Conflict'});
    
    subjectDesignMatrix.PreviousConflict = categorical(subjectDesignMatrix.PreviousConflict);
    subjectDesignMatrix.PreviousConflict = reordercats(subjectDesignMatrix.PreviousConflict, {'NoConflict', 'Conflict'});
    
    % For the current electrode, use a copy of the updated design matrix.
    tempDesignMatrix = subjectDesignMatrix;
    elecPower = currElec.PowerData;  % [nTrials x nTimePoints]
    
    empConflictT = NaN(nTimePoints, 1);
    conflictT_perm = NaN(nTimePoints, n_permutations);
    
    % Loop over time points to run the model
    for t = 1:nTimePoints
        tempDesign = tempDesignMatrix;
        tempDesign.powerValues = elecPower(:, t);
        lme = fitlm(tempDesign, modelFormula);
        empConflictT(t) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict'));
    end
    
    % Permutation testing: shuffle the power values across trials for each time point
    for permi = 1:n_permutations
        for t = 1:nTimePoints
            tempDesign = tempDesignMatrix;
            tempDesign.powerValues = elecPower(randperm(size(elecPower,1)), t);
            lme = fitlm(tempDesign, modelFormula);
            conflictT_perm(t, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict'));
        end
    end
    
    globalMaxConflict = max([max(abs(empConflictT)), max(abs(conflictT_perm(:)))]);
    empConflictTFCE = compute_tfce_1d(empConflictT, globalMaxConflict/nsteps, globalMaxConflict, E, H);
    
    permConflict = zeros(1, n_permutations);
    for permi = 1:n_permutations
        permConflict(permi) = max(compute_tfce_1d(conflictT_perm(:, permi), globalMaxConflict/nsteps, globalMaxConflict, E, H));
    end
    
    thresholdConflict = prctile(permConflict, 95);
    onsetIdx = find(empConflictTFCE > thresholdConflict, 1, 'first');
    if ~isempty(onsetIdx)
        onsetTime = newt(onsetIdx);
    else
        onsetTime = NaN;
    end
    
    result = struct();
    result.Subject = subjID;
    result.Electrode = currElec.Electrode;
    result.empConflictT = empConflictT;
    result.empConflictTFCE = empConflictTFCE;
    result.thresholdConflict = thresholdConflict;
    result.onsetTime_Conflict = onsetTime;
    
    StatsResults.dlPFC{end+1} = result;
end


StatsResults.dmPFC = {};
for i = 1:numel(preparedData.dmPFC)
    currElec = preparedData.dmPFC{i};
    subjID = currElec.Subject;
    subjectDesignMatrix = preparedData.subjectDesignMatrices.(subjID);
    
    %% Update the design matrix as above:
    subjectDesignMatrix.PC = subjectDesignMatrix.PC - 0.33;
    subjectDesignMatrix.CurrentConflict = repmat({'NoConflict'}, height(subjectDesignMatrix), 1);
    subjectDesignMatrix.CurrentConflict(strcmp(subjectDesignMatrix.CurrentType, 'inc')) = {'Conflict'};
    subjectDesignMatrix.PreviousConflict = repmat({'NoConflict'}, height(subjectDesignMatrix), 1);
    subjectDesignMatrix.PreviousConflict(strcmp(subjectDesignMatrix.PreviousType, 'inc')) = {'Conflict'};
    subjectDesignMatrix.CurrentConflict = categorical(subjectDesignMatrix.CurrentConflict);
    subjectDesignMatrix.CurrentConflict = reordercats(subjectDesignMatrix.CurrentConflict, {'NoConflict', 'Conflict'});
    subjectDesignMatrix.PreviousConflict = categorical(subjectDesignMatrix.PreviousConflict);
    subjectDesignMatrix.PreviousConflict = reordercats(subjectDesignMatrix.PreviousConflict, {'NoConflict', 'Conflict'});
    
    tempDesignMatrix = subjectDesignMatrix;
    elecPower = currElec.PowerData;
    
    empConflictT = NaN(nTimePoints, 1);
    conflictT_perm = NaN(nTimePoints, n_permutations);
    
    for t = 1:nTimePoints
        tempDesign = tempDesignMatrix;
        tempDesign.powerValues = elecPower(:, t);
        lme = fitlm(tempDesign, modelFormula);
        empConflictT(t) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict'));
    end
    
    for permi = 1:n_permutations
        for t = 1:nTimePoints
            tempDesign = tempDesignMatrix;
            tempDesign.powerValues = elecPower(randperm(size(elecPower,1)), t);
            lme = fitlm(tempDesign, modelFormula);
            conflictT_perm(t, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict'));
        end
    end
    
    globalMaxConflict = max([max(abs(empConflictT)), max(abs(conflictT_perm(:)))]);
    empConflictTFCE = compute_tfce_1d(empConflictT, globalMaxConflict/nsteps, globalMaxConflict, E, H);
    
    permConflict = zeros(1, n_permutations);
    for permi = 1:n_permutations
        permConflict(permi) = max(compute_tfce_1d(conflictT_perm(:, permi), globalMaxConflict/nsteps, globalMaxConflict, E, H));
    end
    
    thresholdConflict = prctile(permConflict, 97.5);
    onsetIdx = find(abs(empConflictTFCE) > thresholdConflict, 1, 'first');
    if ~isempty(onsetIdx)
        onsetTime = newt(onsetIdx);
    else
        onsetTime = NaN;
    end
    
    result = struct();
    result.Subject = subjID;
    result.Electrode = currElec.Electrode;
    result.empConflictT = empConflictT;
    result.empConflictTFCE = empConflictTFCE;
    result.thresholdConflict = thresholdConflict;
    result.onsetTime_Conflict = onsetTime;
    
    StatsResults.dmPFC{end+1} = result;
end

%% 5. Save the Statistics Results
save('StatsResults.mat', 'StatsResults', 'newt', '-v7.3');
fprintf('Statistical analysis complete. Results saved in StatsResults.mat\n');