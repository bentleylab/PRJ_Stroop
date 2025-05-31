function PermTestLMMMPFCITI(ROI,freqBand,alignEvent)
restoredefaultpath;
root_dir = '/data/user/anaskhan/';
ft_dir=[root_dir 'fieldtrip/'];

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);
addpath(ft_dir);
ft_defaults;

%% Load Data and Preprocess
% load(['/data/user/anaskhan/PRJ_Stroop/results/PermLMMInput/dmPFCThetaITIS.mat'])
load(['/data/user/anaskhan/PRJ_Stroop/results/PermLMMInput/dmPFCHFAITIS.mat'])

GroupDesignMatrixM.RT = round(GroupDesignMatrixM.RT*1000);
GroupDesignMatrixM.PI = zeros(height(GroupDesignMatrixM),1);
GroupDesignMatrixM(strcmpi(GroupDesignMatrixM.BlockType, 'mcon'),:).PI = repmat(0.167,sum(strcmpi(GroupDesignMatrixM.BlockType, 'mcon')),1);
GroupDesignMatrixM(strcmpi(GroupDesignMatrixM.BlockType, 'minc'),:).PI = repmat(0.500,sum(strcmpi(GroupDesignMatrixM.BlockType, 'minc')),1);
GroupDesignMatrixM(strcmpi(GroupDesignMatrixM.BlockType, 'same'),:).PI = repmat(0.333,sum(strcmpi(GroupDesignMatrixM.BlockType, 'same')),1);

GroupDesignMatrixM.PI = GroupDesignMatrixM.PI - 0.33;
%% Define Parameters

n_permutations = 1000;
nsteps = 100; % Number of steps for TFCE
E = 0.66;
H = 2;

numTimePoints = length(newt);

% Initialize storage variables

blockT = NaN(numTimePoints, n_permutations);

empblockT = NaN(numTimePoints, 1);

%% Compute Observed T-Statistics
for tidx = 1:numTimePoints
    powerValues = [];
    for subj = 1:numel(GroupPowerMatrixM)
        powerValues = [powerValues; GroupPowerMatrixM{subj}(:, tidx)];
    end
    
    tempDesignMatrix = GroupDesignMatrixM;
    tempDesignMatrix.powerValues = powerValues;
    
    lme = fitlme(tempDesignMatrix, 'powerValues ~ PI + (1 | Subject)', 'FitMethod', 'REML');
    
    empblockT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PI'));
end

%% Permutation Testing
for permi = 1:n_permutations
    for tidx = 1:numTimePoints
        powerValues = [];
        for subj = 1:numel(GroupPowerMatrixM)
            powerValues = [powerValues; circshift(GroupPowerMatrixM{subj}(:, tidx), randi(height(GroupPowerMatrixM{subj}(:, tidx))))];
        end
        
        tempDesignMatrix = GroupDesignMatrixM;
        tempDesignMatrix.powerValues = powerValues;
        
        lme = fitlme(tempDesignMatrix, 'powerValues ~ PI + (1 | Subject)', 'FitMethod', 'REML');

        blockT(tidx, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PI'));
    end
end

%% Determine Separate Global Maxima for Each Condition
globalMaxBlock = max([max(abs(empblockT(:))), max(abs(blockT(:)))]); 

%% Compute TFCE for Observed Data
empBlockTFCE = compute_tfce_1d(empblockT, globalMaxBlock / nsteps, globalMaxBlock, E, H);

%% Compute TFCE for Permutation Data
permBlock = zeros(1, n_permutations);

for permi = 1:n_permutations
    permBlock(permi) = max(compute_tfce_1d(blockT(:, permi), globalMaxBlock / nsteps, globalMaxBlock, E, H));
end
% %% Save theta Results
% save(['/data/user/anaskhan/PRJ_Stroop/results/PermTestResults/dmPFCThetaITIS.mat'], ...
%     'permBlock', 'empBlockTFCE', 'empblockT', 'newt', '-v7.3');
%% Save hfa Results
save(['/data/user/anaskhan/PRJ_Stroop/results/PermTestResults/dmPFCHFAITIS.mat'], ...
    'permBlock', 'empBlockTFCE', 'empblockT', 'newt', '-v7.3');
end