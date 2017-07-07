function [pval,PLV_diff_hist_sorted] = fn_PLV_time_surr_epoch_diff_rand_PLV_ts(real_plv_diff,plv_ts1,plv_ts2,nboots)
% Calculates the pvalue given the PLV value and phase time series.
% rand_PLV_ts version:
%   1) takes PLV time series from the two windows of interest,
%   2) shuffles the time points between them into random halves
%   3) averages those randomized time series,
%   4) adds difference between these random PLVs to a distribution
%   5) compares real PLV for distribution to calculate a p value
% NOTE: This is 2-tailed, so I make both the real value and the null
% distribution absolute values to avoid +/- differences, just is it different overall

rng('shuffle'); % seed randi with time

% Check if windows are the same length to avoid bias towards PLV of longer windows
if length(plv_ts1) ~= length(plv_ts2)
    error('ERROR: Two plv time series are not the same length.');
end

%Create histogram of PLV values
plv_diff_hist = nan([1 nboots]);

% Combine plv time series
plv_ts_comb = horzcat(plv_ts1,plv_ts2);
if length(plv_ts_comb)<=length(plv_ts1)
    error('ERROR: horzcat(plv_ts1,plv_ts2) did not work; try transposing time series!');
end

for boot = 1:nboots
    % Randomly split the combined time series into two halves
    shuffle_ind1 = randperm(length(plv_ts_comb),length(plv_ts1));
    shuffle_ind2 = setdiff(1:length(plv_ts_comb),shuffle_ind1);
    
    % Calculate PLV and the difference for this iteration
    plv1 = mean(abs(plv_ts_comb(shuffle_ind1)));
    plv2 = mean(abs(plv_ts_comb(shuffle_ind2)));
    plv_diff_hist(boot) = plv1-plv2;
    
    clear shuffle_ind1 shuffle_ind2
end
if sum(isnan(plv_diff_hist))>0
    error('ERROR: Still have NaNs in the plv_diff_hist, check for problems...');
end

% Return PLV p-value
PLV_diff_hist_sorted = sort(plv_diff_hist);
PLV_diff_hist_sorted_abs = abs(PLV_diff_hist_sorted);
PLV_logical = logical(abs(real_plv_diff) > PLV_diff_hist_sorted_abs);
pval = (length(PLV_diff_hist_sorted_abs)-(length(find(PLV_logical))))/nboots;

end


