function filt = zbaseline_HFA(filt, bl_start, bl_end, niter)
%% Z-score HFA bins to baseline period using the bootstrap method

% Get trials into matrix format
tempdata = cat(3,filt.trial{:});

% Number of trials
nTrials = size(tempdata,3);

% Get baseline start and end points in indices
start_idx = dsearchn(filt.time{1}',bl_start);
end_idx = dsearchn(filt.time{1}',bl_end);

% Extract baseline periods
baselines = tempdata(:, start_idx:end_idx, :);

% Pool baselines across trials
pooled_baselines = reshape(baselines, size(baselines,1), []);

%% Begin the random sampling (sample nTrials random data points per iteration)

% Allocate memory for bootstrap distributions

boot_dists = NaN(niter,size(baselines,1));

for ii = 1:niter
    bl_sample_idxs = randsample(length(pooled_baselines),nTrials,true);
    boot_dists(ii,:) = mean(pooled_baselines(:,bl_sample_idxs),2);
end

boot_mean = mean(boot_dists,1)';
boot_sd   = std(boot_dists)';

clear boot_dists


for ii = 1:nTrials
    
    % Normalize the current trial
    filt.trial{ii} = (filt.trial{ii} - boot_mean) ./ boot_sd;
    
end
end
