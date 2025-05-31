function tfr = zfullrecording_tfr(tfr,tfr_full,niter)
%% Z-score HFA bins to whole recording using the bootstrap method

% Get Fs
Fs = 1/(tfr_full.time(2)-tfr_full.time(1));

% Clip edges to account for NaNs
tfr_full.powspctrm = tfr_full.powspctrm(:,:,round(1.5*Fs):end-round(1.5*Fs));

% Number of trials
nTrials = size(tfr.powspctrm,1);

%% Begin the random sampling (sample nTrials random data points per iteration)

% Allocate memory for bootstrap distributions

boot_dists = NaN(niter,size(tfr_full.powspctrm,1),size(tfr_full.powspctrm,2));

for ii = 1:niter
    bl_sample_idxs = randsample(size(tfr_full.powspctrm,3),nTrials,true);
    boot_dists(ii,:,:) = mean(tfr_full.powspctrm(:,:,bl_sample_idxs),3);
end

boot_mean = mean(boot_dists,1);
boot_sd   = std(boot_dists);

if any(any(isnan(boot_mean)))
    error('NaNs at edges. Adjust clipping length.')
end

clear boot_dists


if any(size(boot_mean,[2 3]) ~= size(tfr.powspctrm,[2 3]))
    error('Dimensions of bootstrapped structures and tfr do not match.')
end

tfr.powspctrm = (tfr.powspctrm - boot_mean)./boot_sd;


end