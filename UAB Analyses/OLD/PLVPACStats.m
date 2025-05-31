files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

for f = 1:length(files)
    load(files{f})
    sbj_conflict(:,f) = mean(conflict(lowFreqs >= 4 & lowFreqs <= 8,:));
    sbj_noconflict(:,f) = mean(noconflict(lowFreqs >= 4 & lowFreqs <= 8,:));
end

n_permutations = 1000;
nsteps = 100; 
E = 0.66;
H = 2;

aggregateData = [sbj_conflict,sbj_noconflict];

% Making mapping vector
mapping = [ones(size(aggregateData,2)/2,1); ones(size(aggregateData,2)/2,1)+1];

% Compute the difference for each subject at each time point
diffData = sbj_conflict - sbj_noconflict;  % [nTimepoints x nSubjects]
[nTimepoints, nSubjects] = size(diffData);

% Pre-allocate a matrix to store t statistics for each permutation
tStatsPerm = zeros(nTimepoints, n_permutations);

observedT = mean(diffData,2) ./ (std(diffData,0,2) / sqrt(nSubjects));

for permi = 1:n_permutations

    % Generate a random sign flip for each subject (+1 or -1)
    signFlip = (rand(nSubjects,1) < 0.5) * 2 - 1;  % [nSubjects x 1] vector
    
    % Flip the sign of the differences
    
    permutedDiff = diffData .* signFlip';
    
    % Compute the mean and standard deviation for each time point
    meanDiff = mean(permutedDiff, 2);      % [nTimepoints x 1]
    stdDiff  = std(permutedDiff, 0, 2);      % [nTimepoints x 1]; (0 uses n-1 normalization)
    
    % Compute the t-statistic for each time point:
    tStatsPerm(:, permi) = meanDiff ./ (stdDiff / sqrt(nSubjects));
end

maxVal = max([max(tStatsPerm,[],'all'),max(observedT)]);

observedTFCE = compute_tfce_1d(observedT, maxVal / nsteps, maxVal, E, H);

for permi = 1:n_permutations
    permTmax(permi) = max(compute_tfce_1d(tStatsPerm(:, permi), maxVal / nsteps, maxVal, E, H));
end

cluster_thresh = prctile(permTmax, 100 - (100*0.025));
