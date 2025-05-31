function tfr = zbaseline_tfr(tfr,time, bl_start, bl_end, niter)
%% Z-score HFA bins to baseline period using the bootstrap method

% Number of trials
nTrials = size(tfr.powspctrm,1);

% Get baseline start and end points in indices
start_idx = dsearchn(time',bl_start);
end_idx = dsearchn(time',bl_end);

% Extract baseline periods
baselines = tfr.powspctrm(:,:,:,start_idx:end_idx);
baselines = permute(baselines,[2 3 1 4]);
% Pool baselines across trials
pooled_baselines = reshape(baselines, size(baselines,1), size(baselines,2),[]);

%% Begin the random sampling (sample nTrials random data points per iteration)

% Allocate memory for bootstrap distributions

boot_dists = NaN(niter,size(baselines,1),size(baselines,2));

for ii = 1:niter
    bl_sample_idxs = randsample(length(pooled_baselines),nTrials,true);
    boot_dists(ii,:,:) = mean(pooled_baselines(:,:,bl_sample_idxs),3);
end

boot_mean = mean(boot_dists);
boot_sd   = std(boot_dists);

clear boot_dists pooled_baselines baselines

boot_mean = repmat(boot_mean, [nTrials, 1, 1]);
boot_sd = repmat(boot_sd, [nTrials, 1, 1]);

tfr.powspctrm = (tfr.powspctrm - boot_mean) ./ boot_sd;
end