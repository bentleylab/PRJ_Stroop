function filt = zfullrecording_HFA(filt,filtref,niter)
%% Z-score HFA bins to whole recording using the bootstrap method

baseline = cell2mat(filtref.trial);
% Number of trials
nTrials = length(filt.trial);

%% Begin the random sampling (sample nTrials random data points per iteration)

% Allocate memory for bootstrap distributions

boot_dists = NaN(niter,size(baseline,1));

for ii = 1:niter
    bl_sample_idxs = randsample(length(baseline),nTrials,true);
    boot_dists(ii,:) = mean(baseline(:,bl_sample_idxs),2);
end

boot_mean = mean(boot_dists,1)';
boot_sd   = std(boot_dists)';

clear boot_dists

filt.ztrial = zeros(nTrials, size(baseline,1), length(filt.trial{1,1}));

for ii = 1:nTrials
    % Extract the current trial data
    currentTrialData = filt.trial{ii};
    
    % Normalize the data: (currentTrialData - boot_mean) / boot_sd
    normalizedData = (currentTrialData - boot_mean) ./ boot_sd;
    
    % Assign the normalized data to the 3D matrix, noting the order of dimensions
    filt.ztrial(ii, :, :) = normalizedData;
end
end