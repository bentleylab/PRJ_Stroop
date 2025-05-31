files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

SBJs = cellfun(@(x) x(1:end-4), files, 'UniformOutput', false);
root_dir = '/Users/anaskhan/Desktop/';
%%
GroupDesign = [];
GroupGC = [];
for sbj = 1:numel(SBJs)
    load(files{sbj})
    load(['/Users/anaskhan/Desktop/PRJ_Stroop/results/ThetaISPC/' files{sbj}],'coherenceDesign')

    tempGC = cellfun(@(x, y) log10(x) - log10(y), LtM, MtL, 'UniformOutput', false);
    
    tempGC = cat(3,tempGC{:});
    tempGC = permute(tempGC,[3 1 2]);    

    rowsToremove = strcmpi(coherenceDesign.PreviousType,'None');
    shortDesign = coherenceDesign(coherenceDesign.ElectrodePair == 1,:);
    rowsToremove_gc= strcmpi(shortDesign.PreviousType,'None');

    coherenceDesign(rowsToremove,:) = [];
    tempGC(:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.CurrentType,'neu') | strcmpi(coherenceDesign.PreviousType,'neu');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.CurrentType,'neu') | strcmpi(shortDesign.PreviousType,'neu');
    tempGC(:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.CurrentType,'inc');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.CurrentType,'inc');
    tempGC(:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];


    % Focus on only pI trials
    rowsToremove = strcmpi(coherenceDesign.PreviousType,'con');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.PreviousType,'con');
    tempGC(:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];

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

[~,numFrequencies,~] = size(tempGC);

% Conflict Effect
% conflictP = NaN(numFrequencies, 1);
% conflictB = NaN(numFrequencies, 1);
% conflictT = NaN(numFrequencies, 1);


% Previous effect
% previousP = NaN(numFrequencies, 1);
% previousB = NaN(numFrequencies, 1);

% Previous effect
adjustP = NaN(numFrequencies, 1);
adjustB = NaN(numFrequencies, 1);
adjustT = NaN(numFrequencies, 1);

for freq = 1:numFrequencies

    % Initialize a column for power values
    gcs = [];
    
    for subj = 1:numel(GroupGC)
        % Number of electrode pairs
        n_pairs = size(GroupGC{subj},1);
        for elec = 1:n_pairs
            % Extract the power value for the current electrode, frequency, and time
            gc = squeeze(GroupGC{subj}(elec, freq, :));
            gcs = [gcs; gc];
        end
    end

    
    % Add these power values to the design matrix
    tempDesignMatrix = GroupDesign;
    tempDesignMatrix.GC = gcs;
    
    % Fit the LME model
    % lme = fitlme(tempDesignMatrix, 'GC ~ PreviousType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

    % Conflict Model
    % lme = fitlme(tempDesignMatrix, 'GC ~ CurrentType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

    % Adjustment Model
    lme = fitlme(tempDesignMatrix, 'GC ~ RT + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

    % Extract coefficients and their p-values
    aov = anova(lme,'DFMethod','Satterthwaite');
    
    % Store the p-values and F-stats

    adjustP(freq) = aov.pValue(strcmpi(aov.Term,'RT'));
    adjustB(freq) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'RT'));
    adjustT(freq) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'RT'));

    % conflictP(freq) = aov.pValue(strcmpi(aov.Term,'CurrentType'));
    % conflictB(freq) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
    % conflictT(freq) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));

    
end