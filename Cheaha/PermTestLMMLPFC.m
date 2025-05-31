function PermTestLMMLPFC(ROI,freqBand,alignEvent)
restoredefaultpath;
root_dir = '/data/user/anaskhan/';
ft_dir=[root_dir 'fieldtrip/'];

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);
addpath(ft_dir);
ft_defaults;

%% Load Data and Preprocess
load(['/data/user/anaskhan/PRJ_Stroop/results/PermLMMInput/' ROI freqBand alignEvent '.mat'])

GroupDesignMatrixL.RT = round(GroupDesignMatrixL.RT * 1000);
GroupDesignMatrixL.PC = GroupDesignMatrixL.PC - 0.33;

% Create categorical columns
GroupDesignMatrixL.CurrentConflict = categorical(repmat({'NoConflict'}, height(GroupDesignMatrixL), 1));
GroupDesignMatrixL.CurrentConflict(strcmp(GroupDesignMatrixL.CurrentType, 'inc')) = {'Conflict'};
GroupDesignMatrixL.CurrentConflict = reordercats(GroupDesignMatrixL.CurrentConflict, {'NoConflict', 'Conflict'});

GroupDesignMatrixL.PreviousConflict = categorical(repmat({'NoConflict'}, height(GroupDesignMatrixL), 1));
GroupDesignMatrixL.PreviousConflict(strcmp(GroupDesignMatrixL.PreviousType, 'inc')) = {'Conflict'};
GroupDesignMatrixL.PreviousConflict = reordercats(GroupDesignMatrixL.PreviousConflict, {'NoConflict', 'Conflict'});


n_permutations = 1000;
nsteps = 100; % Number of steps for TFCE
E = 0.66;
H = 2;

numTimePoints = length(newt);

% Initialize storage variables
conflictT = NaN(numTimePoints, n_permutations);
conflictBlockT = NaN(numTimePoints, n_permutations);
sequenceT = NaN(numTimePoints, n_permutations);
previousT = NaN(numTimePoints, n_permutations);
blockT = NaN(numTimePoints, n_permutations);

empconflictT = NaN(numTimePoints, 1);
empconflictBlockT = NaN(numTimePoints, 1);
empsequenceT = NaN(numTimePoints, 1);
emppreviousT = NaN(numTimePoints, 1);
empblockT = NaN(numTimePoints, 1);

%% Compute Observed T-Statistics
for tidx = 1:numTimePoints
    powerValues = [];
    for subj = 1:numel(GroupPowerMatrixL)
        powerValues = [powerValues; GroupPowerMatrixL{subj}(:, tidx)];
    end
    
    tempDesignMatrix = GroupDesignMatrixL;
    tempDesignMatrix.powerValues = powerValues;
    
    lme = fitlme(tempDesignMatrix, ...
        'powerValues ~ CurrentConflict + PreviousConflict + PC + CurrentConflict:PreviousConflict + CurrentConflict:PC + (1 | Subject)', ...
        'FitMethod', 'REML');
    
    empconflictT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict'));
    empconflictBlockT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PC:CurrentConflict_Conflict'));
    empsequenceT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict:PreviousConflict_Conflict'));
    emppreviousT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PreviousConflict_Conflict'));
    empblockT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PC'));
end

%% Permutation Testing
for permi = 1:n_permutations
    for tidx = 1:numTimePoints
        powerValues = [];
        for subj = 1:numel(GroupPowerMatrixL)
            powerValues = [powerValues; circshift(GroupPowerMatrixL{subj}(:, tidx), randi(height(GroupPowerMatrixL{subj}(:, tidx))))];
        end
        
        tempDesignMatrix = GroupDesignMatrixL;
        tempDesignMatrix.powerValues = powerValues;
        
        lme = fitlme(tempDesignMatrix, ...
            'powerValues ~ CurrentConflict + PreviousConflict + PC + CurrentConflict:PreviousConflict + CurrentConflict:PC + (1 | Subject)', ...
            'FitMethod', 'REML');
        
        conflictT(tidx, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict'));
        conflictBlockT(tidx, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PC:CurrentConflict_Conflict'));
        sequenceT(tidx, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'CurrentConflict_Conflict:PreviousConflict_Conflict'));
        previousT(tidx, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PreviousConflict_Conflict'));
        blockT(tidx, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PC'));
    end
end

%% Determine Separate Global Maxima for Each Condition
globalMaxConflictBlock = max([max(abs(empconflictBlockT(:))), max(abs(conflictBlockT(:)))]); 
globalMaxConflict = max([max(abs(empconflictT(:))), max(abs(conflictT(:)))]); 
globalMaxSequence = max([max(abs(empsequenceT(:))), max(abs(sequenceT(:)))]); 
globalMaxPrevious = max([max(abs(emppreviousT(:))), max(abs(previousT(:)))]); 
globalMaxBlock = max([max(abs(empblockT(:))), max(abs(blockT(:)))]); 


%% Compute TFCE for Observed Data
empCBTFCE = compute_tfce_1d(empconflictBlockT, globalMaxConflictBlock / nsteps, globalMaxConflictBlock, E, H);
empConflictTFCE = compute_tfce_1d(empconflictT, globalMaxConflict / nsteps, globalMaxConflict, E, H);
empSeqTFCE = compute_tfce_1d(empsequenceT, globalMaxSequence / nsteps, globalMaxSequence, E, H);
empPrevTFCE = compute_tfce_1d(emppreviousT, globalMaxPrevious / nsteps, globalMaxPrevious, E, H);
empBlockTFCE = compute_tfce_1d(empblockT, globalMaxBlock / nsteps, globalMaxBlock, E, H);

%% Compute TFCE for Permutation Data
permConflictBlock = zeros(1, n_permutations);
permConflict = zeros(1, n_permutations);
permSequence = zeros(1, n_permutations);
permPrevious = zeros(1, n_permutations);
permBlock = zeros(1, n_permutations);

for permi = 1:n_permutations
    permConflictBlock(permi) = max(compute_tfce_1d(conflictBlockT(:, permi), globalMaxConflictBlock / nsteps, globalMaxConflictBlock, E, H));
    permConflict(permi) = max(compute_tfce_1d(conflictT(:, permi), globalMaxConflict / nsteps, globalMaxConflict, E, H));
    permSequence(permi) = max(compute_tfce_1d(sequenceT(:, permi), globalMaxSequence / nsteps, globalMaxSequence, E, H));
    permPrevious(permi) = max(compute_tfce_1d(previousT(:, permi), globalMaxPrevious / nsteps, globalMaxPrevious, E, H));
    permBlock(permi) = max(compute_tfce_1d(blockT(:, permi), globalMaxBlock / nsteps, globalMaxBlock, E, H));
end

%% Save Results
save(['/data/user/anaskhan/PRJ_Stroop/results/PermTestResults/' ROI freqBand alignEvent '.mat'], ...
    'permSequence', 'permConflict', 'permConflictBlock', 'permPrevious', 'permBlock', ...
    'empCBTFCE', 'empSeqTFCE', 'empPrevTFCE', 'empConflictTFCE', 'empBlockTFCE',...
    'empsequenceT', 'empconflictT', 'emppreviousT', 'empconflictBlockT', 'empblockT', 'newt', '-v7.3');
end