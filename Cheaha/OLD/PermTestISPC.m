function PermTestISPC(placeholder1,placeholder2,alignEvent)
restoredefaultpath;
root_dir = '/data/user/anaskhan/';
ft_dir=[root_dir 'fieldtrip/'];

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);
addpath(ft_dir);
ft_defaults;

%% Load Data and Preprocess
load(['/data/user/anaskhan/PRJ_Stroop/results/PermLMMInput/' alignEvent '.mat'])

GroupDesign.RT = round(GroupDesign.RT * 1000);
GroupDesign.PC = GroupDesign.PC - 0.33;

% Create categorical columns
GroupDesign.CurrentConflict = categorical(repmat({'NoConflict'}, height(GroupDesign), 1));
GroupDesign.CurrentConflict(strcmp(GroupDesign.CurrentType, 'inc')) = {'Conflict'};
GroupDesign.CurrentConflict = reordercats(GroupDesign.CurrentConflict, {'NoConflict', 'Conflict'});

GroupDesign.PreviousConflict = categorical(repmat({'NoConflict'}, height(GroupDesign), 1));
GroupDesign.PreviousConflict(strcmp(GroupDesign.PreviousType, 'inc')) = {'Conflict'};
GroupDesign.PreviousConflict = reordercats(GroupDesign.PreviousConflict, {'NoConflict', 'Conflict'});

%% Define Parameters

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
    ispcs = [];
    for subj = 1:numel(GroupISPC)
        ispcs = [ispcs; GroupISPC{subj}(:, tidx)];
    end
    
    tempDesignMatrix = GroupDesign;
    tempDesignMatrix.ispcs = ispcs;
    
    lme = fitlme(tempDesignMatrix, ...
        'ispcs ~ CurrentConflict + PreviousConflict + PC + CurrentConflict:PreviousConflict + CurrentConflict:PC + (1 | Subject)', ...
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
        ispcs = [];
        for subj = 1:numel(GroupISPC)
            ispcs = [ispcs; circshift(GroupISPC{subj}(:, tidx), randi(height(GroupISPC{subj}(:, tidx))))];
        end
        
        tempDesignMatrix = GroupDesign;
        tempDesignMatrix.ispcs = ispcs;
        
        lme = fitlme(tempDesignMatrix, ...
            'ispcs ~ CurrentConflict + PreviousConflict + PC + CurrentConflict:PreviousConflict + CurrentConflict:PC + (1 | Subject)', ...
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
save(['/data/user/anaskhan/PRJ_Stroop/results/PermTestResults/ThetaISPC' alignEvent '.mat'], ...
    'permSequence', 'permConflict', 'permConflictBlock', 'permPrevious', 'permBlock', ...
    'empCBTFCE', 'empSeqTFCE', 'empPrevTFCE', 'empConflictTFCE', 'empBlockTFCE',...
    'empsequenceT', 'empconflictT', 'emppreviousT', 'empconflictBlockT', 'empblockT', 'newt', '-v7.3');
end

save('/Users/anaskhan/Desktop/PRJ_Stroop/results/LMMStats/THetaISPCS.mat', 'permSequence', 'permConflict', 'permConflictBlock', 'permPrevious', 'permBlock', ...
    'empCBTFCE', 'empSeqTFCE', 'empPrevTFCE', 'empConflictTFCE', 'empBlockTFCE',...
    'empsequenceT', 'empconflictT', 'emppreviousT', 'empconflictBlockT', 'empblockT', 'newt', '-v7.3');