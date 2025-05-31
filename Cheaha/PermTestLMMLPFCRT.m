function PermTestLMMLPFCRT(ROI,freqBand,condition)
restoredefaultpath;
root_dir = '/data/user/anaskhan/';
ft_dir=[root_dir 'fieldtrip/'];

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);
addpath(ft_dir);
ft_defaults;

%% Load Data and Preprocess
load(['/data/user/anaskhan/PRJ_Stroop/results/PermLMMInput/dlPFCThetaScCRT.mat'])

GroupDesignMatrixL.RT = round(GroupDesignMatrixL.RT * 1000);

%% Define Parameters

n_permutations = 1000;
nsteps = 100; % Number of steps for TFCE
E = 0.66;
H = 2;

numTimePoints = length(newt);

% Initialize storage variables
RTT = NaN(numTimePoints, n_permutations);

empRTT = NaN(numTimePoints, 1);


%% Compute Observed T-Statistics
for tidx = 1:numTimePoints
    powerValues = [];
    for subj = 1:numel(GroupPowerMatrixL)
        powerValues = [powerValues; GroupPowerMatrixL{subj}(:, tidx)];
    end
    
    tempDesignMatrix = GroupDesignMatrixL;
    tempDesignMatrix.powerValues = powerValues;
    
    % lme = fitlme(tempDesignMatrix, 'powerValues ~ RT + (1 | Subject)', 'FitMethod', 'REML');
    lme = fitlme(tempDesignMatrix, 'powerValues ~ PrevRT + (1 | Subject)', 'FitMethod', 'REML');

    % empRTT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'RT'));
    empRTT(tidx) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'PrevRT'));

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
        
        lme = fitlme(tempDesignMatrix, 'powerValues ~ RT + (1 | Subject)', 'FitMethod', 'REML');
        
        RTT(tidx, permi) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name, 'RT'));

    end
end

%% Determine Separate Global Maxima for Each Condition
globalMaxRT = max([max(abs(empRTT(:))), max(abs(RTT(:)))]); 

%% Compute TFCE for Observed Data
empRTTFCE = compute_tfce_1d(empRTT, globalMaxRT / nsteps, globalMaxRT, E, H);


%% Compute TFCE for Permutation Data
permRT = zeros(1, n_permutations);

for permi = 1:n_permutations
    permRT(permi) = max(compute_tfce_1d(RTT(:, permi), globalMaxRT / nsteps, globalMaxRT, E, H));
end

%% Save Results
save(['/data/user/anaskhan/PRJ_Stroop/results/PermTestResults/dlPFCThetaScCRT.mat'], ...
    'permRT', 'empRTTFCE', 'empRTT', 'newt', '-v7.3');
end