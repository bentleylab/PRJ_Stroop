files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

SBJs = cellfun(@(x) x(1:end-4), files, 'UniformOutput', false);
root_dir = '/Users/anaskhan/Desktop/';
%%
for sbj = 1:numel(SBJs)
    load(files{sbj})
    iI(:,sbj) = squeeze(mean(mean(mean(l_to_m_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'inc') & strcmpi(coherenceDesign.PreviousType,'inc')) - m_to_l_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'inc') & strcmpi(coherenceDesign.PreviousType,'inc')),4),2),1));
    cI(:,sbj) = squeeze(mean(mean(mean(l_to_m_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'inc') & strcmpi(coherenceDesign.PreviousType,'con')) - m_to_l_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'inc') & strcmpi(coherenceDesign.PreviousType,'con')),4),2),1));
    
    iC(:,sbj) = squeeze(mean(mean(mean(l_to_m_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'con') & strcmpi(coherenceDesign.PreviousType,'inc')) - m_to_l_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'con') & strcmpi(coherenceDesign.PreviousType,'inc')),4),2),1));
    cC(:,sbj) = squeeze(mean(mean(mean(l_to_m_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'con') & strcmpi(coherenceDesign.PreviousType,'con')) - m_to_l_GC(:,:,:,strcmpi(coherenceDesign.CurrentType,'con') & strcmpi(coherenceDesign.PreviousType,'con')),4),2),1));
    
end

boundedline(frex,mean(cI,2),std(cI,[],2)./sqrt(size(cI,2)),'r',frex,mean(iI,2),std(iI,[],2)./sqrt(size(iI,2)),'b')
legend({'cI','iI'})

figure(2)
boundedline(frex,mean(iC,2),std(iC,[],2)./sqrt(size(iC,2)),'r',frex,mean(cC,2),std(cC,[],2)./sqrt(size(cC,2)),'b')
legend({'iC','cC'})

%% Stats

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

    rowsToremove = strcmpi(coherenceDesign.PreviousType,'None');
    shortDesign = coherenceDesign(coherenceDesign.ElectrodePair == 1,:);
    rowsToremove_gc= strcmpi(shortDesign.PreviousType,'None');

    coherenceDesign(rowsToremove,:) = [];
    tempGC(:,:,:,rowsToremove_gc) = [];
    shortDesign(rowsToremove_gc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.CurrentType,'neu');% | strcmpi(coherenceDesign.PreviousType,'neu');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_gc = strcmpi(shortDesign.CurrentType,'neu');% | strcmpi(shortDesign.PreviousType,'neu');
    tempGC(:,:,:,rowsToremove_gc) = [];
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

[~,~,numFrequencies,~] = size(tempGC);

% Conflict Effect
conflictP = NaN(numFrequencies, 1);
conflictB = NaN(numFrequencies, 1);
conflictT = NaN(numFrequencies, 1);


% Previous effect
% previousP = NaN(numFrequencies, 1);
% previousB = NaN(numFrequencies, 1);

% Previous effect
% adjustP = NaN(numFrequencies, 1);
% adjustB = NaN(numFrequencies, 1);
% adjustT = NaN(numFrequencies, 1);

for freq = 1:numFrequencies

    % Initialize a column for power values
    gcs = [];
    
    for subj = 1:numel(GroupGC)
        % Number of electrodes for the current subject
        n_MPFC_elecs = size(GroupGC{subj}, 1);
        n_LPFC_elecs = size(GroupGC{subj}, 2);
        for elecM = 1:n_MPFC_elecs
            for elecL = 1:n_LPFC_elecs
                % Extract the power value for the current electrode, frequency, and time
                gc = squeeze(GroupGC{subj}(elecM,elecL, freq, :));
                gcs = [gcs; gc];
            end
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

    conflictP(freq) = aov.pValue(strcmpi(aov.Term,'CurrentType'));
    conflictB(freq) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
    conflictT(freq) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));

    
end
%%
[h,~,~,~] = fdr_bh(adjustP);

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
%% Band-specific

% Initialize a column for power values
gcs = [];

for subj = 1:numel(GroupGC)
    % Number of electrodes for the current subject
    n_MPFC_elecs = size(GroupGC{subj}, 1);
    n_LPFC_elecs = size(GroupGC{subj}, 2);
    for elecM = 1:n_MPFC_elecs
        for elecL = 1:n_LPFC_elecs
            % Extract the power value for the current electrode, frequency, and time
            gc = squeeze(mean(GroupGC{subj}(elecM,elecL, 5:9, :),3));
            gcs = [gcs; gc];
        end
    end
end


% Add these power values to the design matrix
tempDesignMatrix = GroupDesign;
tempDesignMatrix.GC = gcs;

% Fit the LME model
% lme = fitlme(tempDesignMatrix, 'GC ~ PreviousType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

% Conflict model
lme = fitlme(tempDesignMatrix, 'GC ~ CurrentType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');


% Adjustment Model
% lme = fitlme(tempDesignMatrix, 'GC ~ RT + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

% Extract coefficients and their p-values
aov = anova(lme,'DFMethod','Satterthwaite');

% Store the p-values and F-stats

% adjustP = aov.pValue(strcmpi(aov.Term,'RT'));
% adjustB = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'RT'));
% adjustT = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'RT'));


conflictP = aov.pValue(strcmpi(aov.Term,'CurrentType'));
conflictB = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
conflictT = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));

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

%% Band-specific

% Initialize a column for power values
gcs = [];

for subj = 1:numel(GroupGC)
    % Number of electrodes for the current subject
    n_MPFC_elecs = size(GroupGC{subj}, 1);
    n_LPFC_elecs = size(GroupGC{subj}, 2);
    for elecM = 1:n_MPFC_elecs
        for elecL = 1:n_LPFC_elecs
            % Extract the power value for the current electrode, frequency, and time
            gc = squeeze(mean(GroupGC{subj}(elecM,elecL, 5:13, :),3));
            gcs = [gcs; gc];
        end
    end
end


% Add these power values to the design matrix
tempDesignMatrix = GroupDesign;
tempDesignMatrix.GC = gcs;

% Fit the LME model
% lme = fitlme(tempDesignMatrix, 'GC ~ PreviousType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

% Adjustment Model
lme = fitlme(tempDesignMatrix, 'GC ~ PostRT + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

% Extract coefficients and their p-values
aov = anova(lme,'DFMethod','Satterthwaite');

% Store the p-values and F-stats

adjustP = aov.pValue(strcmpi(aov.Term,'PostRT'));
adjustB = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'PostRT'));
adjustT = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'PostRT'));

%% Plot delta predicted values

predictedDelta = predict(lme);
RTquintiles = quantile([tempDesignMatrix.RT],[0.2 0.4 0.6 0.8 1.0]);
RTs = tempDesignMatrix.RT;


% predictedRT = predict(lme);
% Deltaquintiles = quantile([tempDesignMatrix.GC],[0.2 0.4 0.6 0.8 1.0]);
% Deltas = tempDesignMatrix.GC;
% Initialize an array to hold the quintile indices
quintileIndices = zeros(size(RTs));
% quintileIndices = zeros(size(Deltas));

% % Assign each RT to a quintile
for i = 1:length(RTs)
    if RTs(i) <= RTquintiles(1)
        quintileIndices(i) = 1;
    elseif RTs(i) <= RTquintiles(2)
        quintileIndices(i) = 2;
    elseif RTs(i) <= RTquintiles(3)
        quintileIndices(i) = 3;
    elseif RTs(i) <= RTquintiles(4)
        quintileIndices(i) = 4;
    else
        quintileIndices(i) = 5;
    end
end


% Assign each DeltaGC to a quintile
% for i = 1:length(Deltas)
%     if Deltas(i) <= Deltaquintiles(1)
%         quintileIndices(i) = 1;
%     elseif Deltas(i) <= Deltaquintiles(2)
%         quintileIndices(i) = 2;
%     elseif Deltas(i) <= Deltaquintiles(3)
%         quintileIndices(i) = 3;
%     elseif Deltas(i) <= Deltaquintiles(4)
%         quintileIndices(i) = 4;
%     else
%         quintileIndices(i) = 5;
%     end
% end
% Combine predictedDelta and quintileIndices into a table for easier plotting
% Combine predictedDelta and quintileIndices into a table for easier plotting
dataTable = table(predictedDelta, quintileIndices);

% Create a figure
figure;

% Create the boxplot with customized properties
h = boxplot(dataTable.predictedDelta, dataTable.quintileIndices, ...
    'Labels', {'Q1', 'Q2', 'Q3', 'Q4', 'Q5'}, ...
    'Symbol', '', ... % Remove individual data points
    'Whisker', 0); % Remove whiskers

% Set box properties
boxes = findobj(gca, 'Tag', 'Box');
for j = 1:length(boxes)
    patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), 'r', ...
        'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Translucent red box without edge color
end

% Set median line properties
medians = findobj(gca, 'Tag', 'Median');
for j = 1:length(medians)
    medians(j).LineWidth = 2; % Thick median line
end

% Set extreme points
% Find the maximum and minimum points to plot them as extreme points
for j = 1:5
    statVals(j).maxGC = max(predictedDelta(quintileIndices == j));
    statVals(j).minGC = min(predictedDelta(quintileIndices == j));
    statVals(j).medGC = median(predictedDelta(quintileIndices == j));
    statVals(j).IQR = [quantile(predictedDelta(quintileIndices == j),0.25),quantile(predictedDelta(quintileIndices == j),0.75)];
end



% Customize the axes
xlabel('RT Quintile');
ylabel('Predicted Delta');
title('Box Plots of Predicted Delta for Each RT Quintile');

% figure;
% boxplot(log10(dataTable.predictedRT), dataTable.quintileIndices, ...
%     'Labels', {'Q1', 'Q2', 'Q3', 'Q4', 'Q5'});
% xlabel('Delta GC Quintile');
% ylabel('Log(Predicted RT)');
% title('Box Plots of Predicted Delta for Each RT Quintile');

