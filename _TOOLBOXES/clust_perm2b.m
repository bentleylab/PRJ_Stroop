function [pval, t_orig, clust_info, est_alpha]=clust_perm2b(real_diffs, null_diffs, muhats, sigmahats, chan_hood,fwer,tail,thresh_p)
%
% clust_perm2-Independent samples cluster-based permutation test using the "cluster
%             mass" statistic and a null hypothesis of a mean of zero.  This
%             function can handle multiple electrodes and time points/frequencies.  
%             This test was originally proposed for MRI data by Bullmore et al.
%             (1999) and for EEG/MEG analysis by Maris & Oostenveld (2007).
%
% Required Inputs:
%  real_diffs, null_diffs, muhats, sigmahats
%  chan_hood - 2D symmetric binary matrix that indicates which channels are
%              considered neighbors of other channels. E.g., if
%              chan_hood(2,10)=1, then Channel 2 and Channel 10 are
%              nieghbors. You can produce a chan_hood matrix using the
%              function spatial_neighbors.m.
%  fwer            - Desired family-wise error rate (i.e., alpha level) {default=.05}
%  tail            - [1 | 0 | -1] If tail=1, the alternative hypothesis is that the
%                    mean of the data is greater than 0 (upper tailed test).  If tail=0,
%                    the alternative hypothesis is that the mean of the data is different
%                    than 0 (two tailed test).  If tail=-1, the alternative hypothesis
%                    is that the mean of the data is less than 0 (lower tailed test).
%                    {default: 0}
%  thresh_p        - The test-wise p-value threshold for cluster inclusion. If
%                    a channel/time-point has a t-score that corresponds to an
%                    uncorrected p-value greater than thresh_p, it is assigned
%                    a p-value of 1 and not considered for clustering. Note
%                    that thresh_p automatically takes into account the tail of
%                    the test (e.g., you will get positive and negative t-score
%                    thresholds for a two-tailed test).
% Outputs:
%  pval       - p-value at each time point and electrode (corrected for
%                multiple comparisons via the permutation test)
%  t_orig     - t-score at each time point and electrode
%  clust_info - A struct variable containing information about the
%               clusters found.  Depending on the tail of the test it will
%               be composed of all or half of the following fields:
%                 pos_clust_pval: p-values of the positive clusters
%                 pos_clust_mass: t-score mass of the positive clusters
%                 pos_clust_ids:  channel x time point matrix that
%                   indicated which positive cluster each channel/time point
%                   pair belongs to. E.g., if pos_clust_ids(1,2)=4, then
%                   the second time point of the first channel is a
%                   member of the fourth cluster. 0 indicates that the
%                   channel/time point is not a member of any positive
%                   cluster.
%                 neg_clust_pval: p-values of the negative clusters
%                 neg_clust_mass: t-score mass of the negative clusters
%                 neg_clust_ids:  channel x time point matrix that
%                   indicated which negative cluster each channel/time point
%                   pair belongs to. E.g., if neg_clust_ids(1,2)=4, then
%                   the second time point of the first channel is a
%                   member of the fourth cluster. 0 indicates that the
%                   channel/time point is not a member of any negative
%                   cluster.
%  est_alpha  - The estimated family-wise alpha level of the test.  With 
%               permutation tests, a finite number of p-values are possible.
%               This function tries to use an alpha level that is as close 
%               as possible to the desired alpha level.  However, if the 
%               sample size is small, a very limited number of p-values are 
%               possible and the desired family-wise alpha level may be 
%               impossible to acheive.
%
% Author:
% David Groppe
% May, 2011
% Kutaslab, San Diego
%
% References:
% Bullmore, E. T., Suckling, J., Overmeyer, S., Rabe-Hesketh, S., Taylor, 
% E., & Brammer, M. J. (1999). Global, voxel, and cluster tests, by theory 
% and permutation, for a difference between two groups of structural MR 
% images of the brain. IEEE Transactions on Medical Imaging, 18(1), 32-42. 
% doi:10.1109/42.750253
%
% Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of 
% EEG- and MEG-data. Journal of Neuroscience Methods, 164(1), 177-190. 
% doi:10.1016/j.jneumeth.2007.03.024
%
% Manly, B.F.J. (1997) Randomization, bootstrap, and Monte Carlo methods in
% biology. 2nd ed. Chapmn and Hall, London.


n_perm = size(null_diffs,1);

n_chan = size(real_diffs,1);
n_tpt = size(real_diffs,2);
if(n_chan > n_tpt)
  n_chan = size(real_diffs,2);
  n_tpt = size(real_diffs,1);
end

% Factors that are used to compute t-scores.  Saves time to compute them
% now rather than to compute them anew for each permutation.
if tail
    %one tailed test
    thresh_t=norminv(thresh_p,0,1); %note unless thresh_p is greater than .5, thresh_t will be negative
else
    %two tailed test
    thresh_t=norminv(thresh_p/2,0,1);
end

if tail==0
    mx_clust_massNEG=zeros(1,n_perm);
    mx_clust_massPOS=zeros(1,n_perm);
else
    mx_clust_mass=zeros(1,n_perm);
end


print_iters = ceil(n_perm/40);
for perm=1:n_perm

    %compute t-scores
    all_t = (squeeze(null_diffs(perm,:,:))-muhats)./sigmahats;
    
    %form t-scores into clusters
    if tail==0,
        %two-tailed test
        
        %positive clusters
        [clust_ids, n_clust]=find_clusters(all_t,-thresh_t,chan_hood,1); %note, thresh_t should be negative by default
        mx_clust_massPOS(perm)=find_mx_mass(clust_ids,all_t,n_clust,1);
        
        %negative clusters
        [clust_ids, n_clust]=find_clusters(all_t,thresh_t,chan_hood,-1);
        mx_clust_massNEG(perm)=find_mx_mass(clust_ids,all_t,n_clust,-1);
        
    elseif tail>0
        %upper tailed test
        [clust_ids, n_clust]=find_clusters(all_t,-thresh_t,chan_hood,1); %note, thresh_t should be negative by default
        mx_clust_mass(perm)=find_mx_mass(clust_ids,all_t,n_clust,1);
    else
        %lower tailed test
        [clust_ids, n_clust]=find_clusters(all_t,thresh_t,chan_hood,-1); %note, thresh_t should be negative by default
        mx_clust_mass(perm)=find_mx_mass(clust_ids,all_t,n_clust,-1);
    end
    
    if(mod(perm,print_iters) == 0)
        fprintf('.');
    end
end
fprintf('\n');

%Compute critical t's
if tail==0,    
    %two-tailed, test statistic is biggest absolute value of all t's
    mx_clust_mass=abs([mx_clust_massPOS mx_clust_massNEG]);
    tmx_ptile(2)=prctile(mx_clust_mass,100-100*fwer);
    tmx_ptile(1)=-tmx_ptile(2);
    est_alpha=mean(mx_clust_mass>=tmx_ptile(2));
elseif tail==1,
    %upper tailed
    tmx_ptile=prctile(mx_clust_mass,100-100*fwer);
    est_alpha=mean(mx_clust_mass>=tmx_ptile);
else
    %tail=-1, lower tailed
    tmx_ptile=prctile(mx_clust_mass,fwer*100);
    est_alpha=mean(mx_clust_mass<=tmx_ptile);
end

%Compute t-scores of actual observations
t_orig = (real_diffs-muhats)./sigmahats;

%compute p-values
pval=ones(n_chan,n_tpt);
if tail==0,
    %two-tailed test
    pval=pval*2; %default p-value for channel/time point pairs is 2
    
    %positive clusters
    [clust_ids, n_clust]=find_clusters(t_orig,-thresh_t,chan_hood,1); %note thresh_t is negative by default
    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_massPOS>=clust_mass)*2; %multiply by two to correct for doing an upper AND lower tailed test
        pval(use_ids)=clust_p;
        clust_info.pos_clust_pval(a)=clust_p;
        clust_info.pos_clust_mass(a)=clust_mass;
    end
    
    %negative clusters
    [clust_ids, n_clust]=find_clusters(t_orig,thresh_t,chan_hood,-1); %note thresh_t is negative by default
    clust_info.neg_clust_pval=ones(1,n_clust);
    clust_info.neg_clust_mass=zeros(1,n_clust);
    clust_info.neg_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_massNEG<=clust_mass)*2; %multiply by two to correct for doing an upper AND lower tailed test
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
elseif tail>0
    %positive clusters
    [clust_ids, n_clust]=find_clusters(t_orig,-thresh_t,chan_hood,1); %note thresh_t is negative by default
    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_mass>=clust_mass); 
        pval(use_ids)=clust_p;
        clust_info.pos_clust_pval(a)=clust_p;
        clust_info.pos_clust_mass(a)=clust_mass;
    end
else
    %negative clusters
    [clust_ids, n_clust]=find_clusters(t_orig,thresh_t,chan_hood,-1); %note thresh_t is negative by default
    clust_info.neg_clust_pval=ones(1,n_clust);
    clust_info.neg_clust_mass=zeros(1,n_clust);
    clust_info.neg_clust_ids=clust_ids;
    for a=1:n_clust,
        use_ids=find(clust_ids==a);
        clust_mass=sum(t_orig(use_ids));
        clust_p=mean(mx_clust_mass<=clust_mass);
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
end


%%% End of Main Function %%%


function mx_clust_mass=find_mx_mass(clust_ids,data_t,n_clust,tail)

mx_clust_mass=0;
if tail<0
    %looking for most negative cluster mass
    for z=1:n_clust,
        use_ids=(clust_ids==z);
        use_mass=sum(data_t(use_ids));
        if use_mass<mx_clust_mass,
            mx_clust_mass=use_mass;
        end
    end
elseif tail>0,
    %looking for most positive cluster mass
    for z=1:n_clust,
        use_ids=(clust_ids==z);
        use_mass=sum(data_t(use_ids));
        if use_mass>mx_clust_mass,
            mx_clust_mass=use_mass;
        end
    end
end

