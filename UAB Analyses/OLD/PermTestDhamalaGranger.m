files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
%%
for subj = 1:length(files)
    load(files{subj})
    net_incLtM = granger_inc(:,:,2,1)-granger_inc(:,:,1,2);
    net_incLtM = net_incLtM(:,dsearchn(time',-0.2):dsearchn(time',1.25));
    
    net_conLtM = granger_con(:,:,2,1)-granger_con(:,:,1,2);
    net_conLtM = net_conLtM(:,dsearchn(time',-0.2):dsearchn(time',1.25));

    con(:,:,subj) = net_conLtM;
    inc(:,:,subj) = net_incLtM;
end
time = time(dsearchn(time',-0.2):dsearchn(time',1.25));
% Number of subjects
nSubjects = 8;

% Generate all possible binary combinations
permutations = dec2bin(1:(2^nSubjects-1)) - '0';

allSamples = cat(4,con,inc); % nFrequencies x nTime x nSubjects x nConditions
allSamples = squeeze(mean(allSamples(frex >= 4 & frex <= 8,:,:,:),1)); % Average over theta frequencies

trueDiff = mean(allSamples(:,:,2) - allSamples(:,:,1),2);
alpha = 0.05;
npermutations = size(permutations, 1);
perm_dist = NaN(npermutations,length(trueDiff));

% Apply permutations
for permi = 1:npermutations
    perm = permutations(permi, :);
    
    % Shuffle condition labels within subset of subjects
    for j = 1:nSubjects
        if perm(j) == 1
            sbj_diffs(:,j) = allSamples(:,j,1) - allSamples(:,j,2);
        else
            sbj_diffs(:,j) = allSamples(:,j,2) - allSamples(:,j,1);
        end
    end
    
    % Compute mean over subject differences
    perm_dist(permi,:) = mean(sbj_diffs,2);
end

% Get mean and standard deviation of chance distribution
mean_h0 = squeeze(mean(perm_dist))';
std_h0 = squeeze(std(perm_dist))';

% z-score empirical values using mean and std of H0 dist
zmap = (trueDiff-mean_h0) ./ std_h0;

% convert to Z value
zcrit = abs(norminv(alpha/2));

% threshold at critical p-value by setting subthreshold values to 0
zmap(abs(zmap)<zcrit) = 0;
zmap = abs(zmap);

%%% cluster-based correction

% initialize matrices for cluster-based correctiion
max_cluster_masses = zeros(1,npermutations);

% loop through permutations
for permi = 1:npermutations

    % take each perm map, and transform to Z
    threshimg = squeeze(perm_dist(permi,:))';
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % threshold iamge at p-value
    threshimg(abs(threshimg)<zcrit) = 0;
    threshimg = abs(threshimg);
    % find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0

        for ii = 1:islands.NumObjects

            % sum z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));

        end
    
        % store size of biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% find cluster threshold based using 95th percentile of cluster masses
% dist

cluster_thresh = prctile(max_cluster_masses,100-(100*0.05));

% set subthreshold clusters in real thresholded zmap to 0
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects

    %if real cluster masses are too small remove them by setting to 0
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end
%% T-F full map
for subj = 1:length(files)
    load(files{subj})
    net_incLtM = granger_inc(:,:,2,1)-granger_inc(:,:,1,2);
    net_incLtM = net_incLtM(:,nearest(time,-0.2):nearest(time,1.25));
    
    net_conLtM = granger_con(:,:,2,1)-granger_con(:,:,1,2);
    net_conLtM = net_conLtM(:,nearest(time,-0.2):nearest(time,1.25));

    con(:,:,subj) = squeeze(net_conLtM);
    inc(:,:,subj) = squeeze(net_incLtM);
end

time = time(nearest(time,-0.2):nearest(time,1.25));
nSubjects = 8;
nFrequencies = size(con, 1);
nTimePoints = length(time);

% Generate all possible binary combinations
permutations = dec2bin(1:(2^nSubjects-1)) - '0';

% Combine condition data into one array
allSamples = cat(4, con, inc); % nFrequencies x nTimePoints x nSubjects x 2Conditions

% Compute the true difference (inc - con) at each time-frequency point
trueDiff = mean(allSamples(:,:,:,2) - allSamples(:,:,:,1), 3);

alpha = 0.05;
npermutations = length(permutations);
perm_dist = NaN(npermutations, nFrequencies, nTimePoints);

% Apply permutations across frequencies and time points
for permi = 1:npermutations
    perm = permutations(permi, :);
    
    % Shuffle condition labels within subjects
    for j = 1:nSubjects
        if perm(j) == 1
            sbj_diffs(:,:,j) = allSamples(:,:,j,1) - allSamples(:,:,j,2);
        else
            sbj_diffs(:,:,j) = allSamples(:,:,j,2) - allSamples(:,:,j,1);
        end
    end
    
    % Compute the mean over subjects for the permutation
    perm_dist(permi,:,:) = mean(sbj_diffs, 3);
end

% Compute the mean and standard deviation of the permutation distribution
mean_h0 = squeeze(mean(perm_dist,1));
std_h0 = squeeze(std(perm_dist,[],1));

% Z-score the true differences using the null distribution
zmap = (trueDiff - mean_h0) ./ std_h0;

% Set the critical z-value for alpha = 0.05
zcrit = abs(norminv(alpha/2));

% Threshold the zmap at the critical z-value
zmap(abs(zmap) < zcrit) = 0;

%%% Cluster-based correction

% Initialize array for max cluster masses from permutations
max_cluster_masses = zeros(1,npermutations);

% Loop through permutations for cluster-based correction
for permi = 1:npermutations

    % Z-score the permutation map
    threshimg = squeeze(perm_dist(permi,:,:));
    threshimg = (threshimg - mean_h0) ./ std_h0;

    % Threshold at the critical z-value
    threshimg(abs(threshimg) < zcrit) = 0;
    
    % Find clusters
    islands = bwconncomp(threshimg);

    if numel(islands.PixelIdxList) > 0
        tempclustmasses = zeros(1,islands.NumObjects);
        for ii = 1:islands.NumObjects
            % Sum the z-statistic mass of clusters
            tempclustmasses(ii) = sum(threshimg(islands.PixelIdxList{ii}));
        end
        
        % Store the size of the biggest cluster
        max_cluster_masses(permi) = max(tempclustmasses);
    end
end

% Find cluster size threshold at 95th percentile of the null distribution
cluster_thresh = prctile(max_cluster_masses, 95);

% Apply the threshold to the real zmap
islands = bwconncomp(zmap);

for ii = 1:islands.NumObjects
    % If real cluster mass is smaller than threshold, remove it
    if sum(zmap(islands.PixelIdxList{ii})) < cluster_thresh
        zmap(islands.PixelIdxList{ii}) = 0;
    else
        zmap(islands.PixelIdxList{ii}) = 1;
    end
end
