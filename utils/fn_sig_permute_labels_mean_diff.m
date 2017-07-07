function [pval,diff_hist_sorted] = fn_sig_permute_labels_mean_diff(dist1,dist2,nboots)
% Calculates the two-sided pvalue for difference in means between two distributions.
% Shuffles the labels of data across the two distributions, then averages
% those randomized distributions; repeat this process nboots times
%   1) calculate actual difference in means between two distributions
%   2) shuffles labels of values in two distributions
%   3) averages those randomized distributions
%   4) adds difference between these randomized means to a distribution
%   5) compares real difference to sorted, absolute values of permuted distribution to obtain p value
% NOTE: This is 2-tailed, so I make both the real value and the null
% distribution absolute values to avoid +/- differences, just is it different overall

rng('shuffle'); % seed randi with time

% Check if windows are the same length to avoid bias towards longer windows
%if length(dist1) ~= length(dist2)
%    disp('WARNING: Two distributions are not the same size; using size of smaller one.');
%end
min_sample_len = min(length(dist1),length(dist2));
real_diff = mean(dist1)-mean(dist2);

% Initialize histogram of mean differences
diff_hist = nan([1 nboots]);

% Combine time series
if size(dist1,1)==1
    ts_comb = horzcat(dist1,dist2);
else
    ts_comb = vertcat(dist1,dist2);
end
if length(ts_comb)<=length(dist1)
    error('ERROR: horzcat(dist1,dist2) did not work; try transposing time series!');
end

for boot = 1:nboots
    % Randomly split the combined time series into two halves
    shuffle_ind1 = randperm(length(ts_comb),min_sample_len);
    shuffle_ind2 = randsample(setdiff(1:length(ts_comb),shuffle_ind1),min_sample_len);
    
    % Calculate means and their difference for this iteration
    diff_hist(boot) = mean(ts_comb(shuffle_ind1))-mean(ts_comb(shuffle_ind2));
    
    clear shuffle_ind1 shuffle_ind2
end
if sum(isnan(diff_hist))>0
    error('ERROR: Still have NaNs in the diff_hist, check for problems...');
end

% Return p-value
diff_hist_sorted = sort(diff_hist);
diff_hist_sorted_abs = abs(diff_hist_sorted);
diff_logical = logical(abs(real_diff) > diff_hist_sorted_abs);
pval = (length(diff_hist_sorted_abs)-(length(find(diff_logical))))/nboots;

end


