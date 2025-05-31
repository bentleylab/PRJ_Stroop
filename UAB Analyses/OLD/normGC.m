files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

SBJs = cellfun(@(x) x(1:end-4), files, 'UniformOutput', false);
root_dir = '/Users/anaskhan/Desktop/';
%% z-score GC values
for sbj = 1:numel(SBJs)
    load(files{sbj});

    LtM = cellfun(@(x) squeeze(sum(x(:,5:9,:),2)), LtM, 'UniformOutput', false);
    MtL = cellfun(@(x) squeeze(sum(x(:,5:9,:),2)), MtL, 'UniformOutput', false);

    tempGC = cellfun(@(x, y) log10(x) - log10(y), LtM, MtL, 'UniformOutput', false);
    
    tempGC = cat(3,tempGC{:});
    tempGC = permute(tempGC,[3 1 2]);

    baselines = reshape(tempGC(:,1:75,:),size(tempGC,1),[]);
   
    niter = 1000;
    boot_dists = NaN(niter,size(baselines,1));

    for ii = 1:niter
        bl_sample_idxs = randsample(length(baselines),size(tempGC,3),true);
        boot_dists(ii,:) = mean(baselines(:,bl_sample_idxs),2);
    end
    
    boot_mean = mean(boot_dists,1)';
    boot_sd   = std(boot_dists)';
    
    clear boot_dists

    zLtM = (tempGC - boot_mean)./boot_sd;
    save(files{sbj},'zLtM','-append')

end
%%

inc = [];
con = [];
for sbj = 1:numel(SBJs)
    load(files{sbj},'zLtM','grangerTime');
    load(['/Users/anaskhan/Desktop/PRJ_Stroop/results/ThetaISPC/' files{sbj}],'coherenceDesign');

    % LtM = cellfun(@(x) squeeze(sum(x(:,5:9,:),2)), LtM, 'UniformOutput', false);
    % MtL = cellfun(@(x) squeeze(sum(x(:,5:9,:),2)), MtL, 'UniformOutput', false);
    % 
    % tempGC = cellfun(@(x, y) x - y, LtM, MtL, 'UniformOutput', false);
    % 
    % tempGC = cat(3,tempGC{:});
    % tempGC = permute(tempGC,[2 3 1]);

    trials = coherenceDesign([coherenceDesign.ElectrodePair] == 1,:);

    % con(:,sbj) = mean(mean(tempGC(strcmpi(trials.CurrentType,'con'),:,:)),2);
    % inc(:,sbj) = mean(mean(tempGC(strcmpi(trials.CurrentType,'inc'),:,:)),2);

    zLtM = permute(zLtM,[3 1 2]);
    con(:,sbj) = mean(mean(zLtM(strcmpi(trials.CurrentType,'con'),:,:)),2);
    inc(:,sbj) = mean(mean(zLtM(strcmpi(trials.CurrentType,'inc'),:,:)),2);
end

figure(1)
boundedline(grangerTime,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',grangerTime,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r')
legend({'iC','cC'})

figure(2)
boundedline(frex,mean(iC,2),std(iC,[],2)./sqrt(size(iC,2)),'r',frex,mean(cC,2),std(cC,[],2)./sqrt(size(cC,2)),'b')
legend({'iC','cC'})

%% Stats

GroupDesign = [];
GroupGC = [];
for sbj = 1:numel(SBJs)
    load(files{sbj})
    load(['/Users/anaskhan/Desktop/PRJ_Stroop/results/ThetaISPC/' files{sbj}],'coherenceDesign');

    LtM = cellfun(@(x) squeeze(sum(x(:,5:9,:),2)), LtM, 'UniformOutput', false);
    MtL = cellfun(@(x) squeeze(sum(x(:,5:9,:),2)), MtL, 'UniformOutput', false);

    tempGC = cellfun(@(x, y) log10(x) - log10(y), LtM, MtL, 'UniformOutput', false);
    
    tempGC = cat(3,tempGC{:});
    tempGC = permute(tempGC,[2 3 1]);

    % rowsToremove = strcmpi(coherenceDesign.PreviousType,'None');
    shortDesign = coherenceDesign(coherenceDesign.ElectrodePair == 1,:);
    % rowsToremove_gc= strcmpi(shortDesign.PreviousType,'None');

    % coherenceDesign(rowsToremove,:) = [];
    % tempGC(rowsToremove_gc,:,:) = [];
    % shortDesign(rowsToremove_gc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.CurrentType,'neu');% | strcmpi(coherenceDesign.PreviousType,'neu');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.CurrentType,'neu');% | strcmpi(shortDesign.PreviousType,'neu');
    tempGC(rowsToremove_gc,:,:) = [];
    shortDesign(rowsToremove_gc,:) = [];

    % rowsToremove = strcmpi(coherenceDesign.CurrentType,'inc');
    % coherenceDesign(rowsToremove,:) = [];
    % 
    % rowsToremove_gc = strcmpi(shortDesign.CurrentType,'inc');
    % tempGC(:,:,:,rowsToremove_gc) = [];
    % shortDesign(rowsToremove_gc,:) = [];


    % Focus on only pI trials
    % rowsToremove = strcmpi(coherenceDesign.PreviousType,'con');
    % coherenceDesign(rowsToremove,:) = [];
    % 
    % rowsToremove_gc = strcmpi(shortDesign.PreviousType,'con');
    % tempGC(:,:,:,rowsToremove_gc) = [];
    % shortDesign(rowsToremove_gc,:) = [];

    coherenceDesign.Subject = repmat(sbj,size(coherenceDesign,1),1);
    GroupDesign = cat(1,GroupDesign,coherenceDesign);
    GroupGC{sbj} = tempGC;
end

GroupDesign.RT = round(GroupDesign.RT*1000);
%%
GroupDesign.CurrentType = categorical(GroupDesign.CurrentType);
GroupDesign.CurrentType = reordercats(GroupDesign.CurrentType, {'con', 'inc'});
%%
GroupDesign.PreviousType = categorical(GroupDesign.PreviousType);
GroupDesign.PreviousType = reordercats(GroupDesign.PreviousType, {'con', 'inc'});
%%

[~,~,numTimePoints] = size(tempGC);

% Conflict Effect
conflictP = NaN(numTimePoints, 1);
conflictB = NaN(numTimePoints, 1);
conflictT = NaN(numTimePoints, 1);


% Previous effect
% previousP = NaN(numFrequencies, 1);
% previousB = NaN(numFrequencies, 1);

% Previous effect
% adjustP = NaN(numFrequencies, 1);
% adjustB = NaN(numFrequencies, 1);
% adjustT = NaN(numFrequencies, 1);

for time = 1:numTimePoints

    % Initialize a column for power values
    gcs = [];
    
    for subj = 1:numel(GroupGC)
        % Number of electrodes pairs for current subject
        n_pairs = size(GroupGC{subj},2);
        for pair = 1:n_pairs
            
            % Extract the power value for the current electrode, frequency, and time
            gc = squeeze(GroupGC{subj}(:,pair,time));
            gcs = [gcs; gc];

        end
    end

    
    % Add these power values to the design matrix
    tempDesignMatrix = GroupDesign;
    tempDesignMatrix.GC = gcs;
    
    % Fit the LME model
    % lme = fitlme(tempDesignMatrix, 'GC ~ PreviousType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

    % Conflict Model
    lme = fitlme(tempDesignMatrix, 'GC ~ CurrentType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

    % Adjustment Model
    % lme = fitlme(tempDesignMatrix, 'GC ~ RT + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

    % Extract coefficients and their p-values
    aov = anova(lme,'DFMethod','Satterthwaite');
    
    % Store the p-values and F-stats

    % adjustP(freq) = aov.pValue(strcmpi(aov.Term,'RT'));
    % adjustB(freq) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'RT'));
    % adjustT(freq) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'RT'));

    conflictP(time) = aov.pValue(strcmpi(aov.Term,'CurrentType'));
    conflictB(time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
    conflictT(time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));

    
end
%%
threshP = conflictP < 0.05;
islands = bwconncomp(threshP);

sig_clusters = find(cellfun(@(x) length(x),islands.PixelIdxList) > 25);

for ii = 1:length(sig_clusters)
    sig_idxs{ii} = islands.PixelIdxList{sig_clusters(ii)};
end

plot(frex,adjustT,'-k','LineWidth',1)

% Calculate the y-value for the horizontal line
min_adjustB = min(adjustB);
y_line = min_adjustB * 0.9;

% Find the indices where h == 1
indices = h == 1;
hold on
% Plot the horizontal line only at those indices
plot(frex(indices), y_line * ones(size(frex(indices))), '-r');
hold off
%% Conflict
[h,~,~,~] = fdr_bh(conflictP);

plot(frex,conflictT,'-k','LineWidth',1)

% Calculate the y-value for the horizontal line
min_adjustB = min(adjustB);
y_line = min_adjustB * 0.9;

% Find the indices where h == 1
indices = h == 1;
hold on
% Plot the horizontal line only at those indices
plot(frex(indices), y_line * ones(size(frex(indices))), '-r');
hold off

%% Response Locked granger (baseline period)

GroupDesign = [];
GroupGC = [];
for sbj = 1:numel(SBJs)
    load(files{sbj})
    tempGC = log10(l_to_m_GC) - log10(m_to_l_GC);
    clear coherenceDesign
    SBJ = SBJs{sbj};
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    [coherenceDesign,~] = fn_create_LMM_design(trial_info, size(tempGC,2), size(tempGC,1),'coherence');

    rowsToremove = strcmpi(coherenceDesign.PostType,'None');
    shortDesign = coherenceDesign(coherenceDesign.ElectrodePair == 1,:);
    rowsToremove_gc= strcmpi(shortDesign.PostType,'None');

    coherenceDesign(rowsToremove,:) = [];
    tempGC(:,:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.PostType,'neu') | strcmpi(coherenceDesign.CurrentType,'neu');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.PostType,'neu') | strcmpi(shortDesign.CurrentType,'neu');
    tempGC(:,:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.CurrentType,'con');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.CurrentType,'con');
    tempGC(:,:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];


    % Focus on only iC trials
    rowsToremove = strcmpi(coherenceDesign.PostType,'con');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.PostType,'con');
    tempGC(:,:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];

    coherenceDesign.Subject = repmat(sbj,size(coherenceDesign,1),1);
    GroupDesign = cat(1,GroupDesign,coherenceDesign);
    GroupGC{sbj} = tempGC;
end

GroupDesign.PostRT = round(GroupDesign.PostRT*1000);

