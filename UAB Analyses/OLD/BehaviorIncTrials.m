files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

meanRT_cI = [];
meanRT_iI = [];

for subj = 1:numel(files)
    load(files{subj}, "lpfcDesign");

    lpfcDesign(lpfcDesign.Electrode ~= 1,:) = [];
    % Remove rows from design matrix that have previous type "None"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'None');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that have previous or current type "neu"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'neu') | strcmpi(lpfcDesign.CurrentType, 'neu');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that are current type "con"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'con');
    lpfcDesign(rowsToRemove, :) = [];
    
    % Calculate mean RT for congruent (cI) trials
    cI_index = strcmpi(lpfcDesign.PreviousType, 'con');
    meanRT_cI(subj) = round(mean(lpfcDesign.RT(cI_index)*1000));

    % Calculate mean RT for incongruent (iI) trials
    iI_index = strcmpi(lpfcDesign.PreviousType, 'inc');
    meanRT_iI(subj) = round(mean(lpfcDesign.RT(iI_index)*1000));
end

% Calculate mean and SEM for each condition
mean_cI = mean(meanRT_cI);
mean_iI = mean(meanRT_iI);
sem_cI = std(meanRT_cI) / sqrt(length(meanRT_cI));
sem_iI = std(meanRT_iI) / sqrt(length(meanRT_iI));

% Create a figure for the plot
figure;

% Plot the means as bars
b = bar([1, 2], [mean_cI, mean_iI], 'FaceColor', 'flat', 'BarWidth', 0.5); % Adjust BarWidth for skinnier bars
b.CData(1,:) = [0.75, 0.75, 1];% Blue color for cI
b.CData(2,:) = [1, 0.75, 0.75]; % Red color for iI

b.EdgeColor = 'black'; % Black edge color
b.LineWidth = 1.5; % Make edges thicker

hold on;

% Add error bars for SEM
errorbar([1, 2], [mean_cI, mean_iI], [sem_cI, sem_iI], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Overlay individual data points and connect them
for subj = 1:length(meanRT_cI)
    plot([1.1, 1.9], [meanRT_cI(subj), meanRT_iI(subj)], '-k', 'LineWidth', 1); % Connect with black lines
end

% Enhance plot
set(gca, 'XTick', [1 2], 'XTickLabel', {'cI', 'iI'});
ylabel('Mean RT (ms)');
title('Group Level Gratton Effect: Incongruent');
ylim([min([meanRT_cI, meanRT_iI]) * 0.95, max([meanRT_cI, meanRT_iI]) * 1.05]); % Adjust y-limits slightly
yticks([600 1000 1400])
%%
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];

files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
files = cellfun(@(x) strrep(x, '.mat', ''), files, 'UniformOutput', false);

meanRT_mcon_cI = [];
meanRT_mcon_iI = [];
meanRT_minc_cI = [];
meanRT_minc_iI = [];

for subj = 1:numel(files)

    SBJ = files{subj};

    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

    lpfcDesign = fn_create_PERM_design(trial_info);

    % Remove rows from design matrix that have previous type "None"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'None');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that have previous or current type "neu"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'neu') | strcmpi(lpfcDesign.CurrentType, 'neu');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that are current type "con"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'con');
    lpfcDesign(rowsToRemove, :) = [];

    % Filter for 'mcon' BlockType
    mcon_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "mcon"), :);

    % Extract cI and iI pairs for 'mcon'
    mcon_cI = mcon_trials(strcmpi(mcon_trials.PreviousType, "con"), :);
    mcon_iI = mcon_trials(strcmpi(mcon_trials.PreviousType, "inc"), :);
    
    % Filter for 'minc' BlockType
    minc_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "minc"), :);
    % Extract cI and iI pairs for 'minc'
    minc_cI = minc_trials(strcmpi(minc_trials.PreviousType, "con"), :);
    minc_iI = minc_trials(strcmpi(minc_trials.PreviousType, "inc"), :);

    if isempty(mcon_iI)
        continue
    end
    meanRT_mcon_cI(subj) = round(mean(mcon_cI.RT * 1000));
    meanRT_mcon_iI(subj) = round(mean(mcon_iI.RT * 1000));
    meanRT_minc_cI(subj) = round(mean(minc_cI.RT * 1000));
    meanRT_minc_iI(subj) = round(mean(minc_iI.RT * 1000));
end

clearvars -except meanRT_minc_iI meanRT_minc_cI meanRT_mcon_iI meanRT_mcon_cI

meanRT_mcon_cI([2 9]) = [];
meanRT_mcon_iI([2 9]) = [];
meanRT_minc_cI([2 9]) = [];
meanRT_minc_iI([2 9]) = [];

% Calculate mean and SEM for mcon
mean_mcon_cI = mean(meanRT_mcon_cI);
mean_mcon_iI = mean(meanRT_mcon_iI);
sem_mcon_cI = std(meanRT_mcon_cI) / sqrt(length(meanRT_mcon_cI));
sem_mcon_iI = std(meanRT_mcon_iI) / sqrt(length(meanRT_mcon_iI));

% Calculate mean and SEM for minc
mean_minc_cI = mean(meanRT_minc_cI);
mean_minc_iI = mean(meanRT_minc_iI);
sem_minc_cI = std(meanRT_minc_cI) / sqrt(length(meanRT_minc_cI));
sem_minc_iI = std(meanRT_minc_iI) / sqrt(length(meanRT_minc_iI));

% Plot for mcon
figure;
b_mcon = bar([1, 2], [mean_mcon_cI, mean_mcon_iI], 'FaceColor', 'flat', 'BarWidth', 0.5);
b_mcon.CData(1,:) = [1, 0.75, 0.75]; % Pink for cI
b_mcon.CData(2,:) = [0.75, 0.75, 1]; % Blue for iI
b_mcon.EdgeColor = 'black';
b_mcon.LineWidth = 1.5;
hold on;
errorbar([1, 2], [mean_mcon_cI, mean_mcon_iI], [sem_mcon_cI, sem_mcon_iI], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
for subj = 1:length(meanRT_mcon_cI)
    plot([1.1, 1.9], [meanRT_mcon_cI(subj), meanRT_mcon_iI(subj)], '-k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'cI', 'iI'});
ylabel('Mean Reaction Time (ms)');
title('Mostly Congruent: cI vs iI');
ylim([min([meanRT_mcon_cI, meanRT_mcon_iI]) * 0.95, max([meanRT_mcon_cI, meanRT_mcon_iI]) * 1.05]);

% Plot for minc
figure;
b_minc = bar([1, 2], [mean_minc_cI, mean_minc_iI], 'FaceColor', 'flat', 'BarWidth', 0.5);
b_minc.CData(1,:) = [1, 0.75, 0.75]; % Pink for cI
b_minc.CData(2,:) = [0.75, 0.75, 1]; % Blue for iI
b_minc.EdgeColor = 'black';
b_minc.LineWidth = 1.5;
hold on;
errorbar([1, 2], [mean_minc_cI, mean_minc_iI], [sem_minc_cI, sem_minc_iI], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
for subj = 1:length(meanRT_minc_cI)
    plot([1.1, 1.9], [meanRT_minc_cI(subj), meanRT_minc_iI(subj)], '-k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'cI', 'iI'});
ylabel('Mean Reaction Time (ms)');
title('Mostly Incongruent: cI vs iI');
ylim([min([meanRT_minc_cI, meanRT_minc_iI]) * 0.95, max([meanRT_minc_cI, meanRT_minc_iI]) * 1.05]);

%% "Same" BlockType

root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];

files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
files = cellfun(@(x) strrep(x, '.mat', ''), files, 'UniformOutput', false);

meanRT_same_cI = [];
meanRT_same_iI = [];

for subj = 1:numel(files)
    SBJ = files{subj};

    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

    lpfcDesign = fn_create_PERM_design(trial_info);

    % Remove rows from design matrix that have previous type "None"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'None');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that have previous or current type "neu"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'neu') | strcmpi(lpfcDesign.CurrentType, 'neu');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that are current type "con"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'con');
    lpfcDesign(rowsToRemove, :) = [];

    % Filter for 'same' BlockType
    same_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "same"), :);

    % Extract cI and iI pairs for 'same'
    same_cI = same_trials(strcmpi(same_trials.PreviousType, "con"), :);
    same_iI = same_trials(strcmpi(same_trials.PreviousType, "inc"), :);

    % Calculate the mean reaction time for cI and iI in 'same'
    meanRT_same_cI(subj) = round(mean(same_cI.RT * 1000));
    meanRT_same_iI(subj) = round(mean(same_iI.RT * 1000));
end

clearvars -except meanRT_same_iI meanRT_same_cI

% Calculate mean and SEM for 'same'
mean_same_cI = mean(meanRT_same_cI);
mean_same_iI = mean(meanRT_same_iI);
sem_same_cI = std(meanRT_same_cI) / sqrt(length(meanRT_same_cI));
sem_same_iI = std(meanRT_same_iI) / sqrt(length(meanRT_same_iI));

% Plot for 'same'
figure;
b_same = bar([1, 2], [mean_same_cI, mean_same_iI], 'FaceColor', 'flat', 'BarWidth', 0.5);
b_same.CData(1,:) = [1, 0.75, 0.75]; % Pink for cI
b_same.CData(2,:) = [0.75, 0.75, 1]; % Blue for iI
b_same.EdgeColor = 'black';
b_same.LineWidth = 1.5;
hold on;
errorbar([1, 2], [mean_same_cI, mean_same_iI], [sem_same_cI, sem_same_iI], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
for subj = 1:length(meanRT_same_cI)
    plot([1.1, 1.9], [meanRT_same_cI(subj), meanRT_same_iI(subj)], '-k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'cI', 'iI'});
ylabel('Mean Reaction Time (ms)');
title('Equal Congruent: cI vs iI');
ylim([min([meanRT_same_cI, meanRT_same_iI]) * 0.95, max([meanRT_same_cI, meanRT_same_iI]) * 1.05]);
%% LMM for RTs
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

GroupDesign = [];

for subj = 1:numel(files)
    load(files{subj}, "lpfcDesign");

    lpfcDesign(lpfcDesign.Electrode ~= 1,:) = [];
    % Remove rows from design matrix that have previous type "None"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'None');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that have previous or current type "neu"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'neu') | strcmpi(lpfcDesign.CurrentType, 'neu');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that are current type "con"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'con');
    lpfcDesign(rowsToRemove, :) = [];
    
    lpfcDesign.RT = round(lpfcDesign.RT * 1000);

    GroupDesign = [GroupDesign; lpfcDesign];
end

clearvars -except GroupDesign

GroupDesign.PC = GroupDesign.PC - mean(GroupDesign.PC);
GroupDesign.PreviousType = categorical(GroupDesign.PreviousType);
GroupDesign.PreviousType = reordercats(GroupDesign.PreviousType, {'con', 'inc'});

lme = fitlme(GroupDesign, 'RT ~ PreviousType + PC + (1 | Subject)', 'FitMethod', 'REML');
altlme = fitlme(GroupDesign, 'RT ~ PreviousType + PC + PreviousType*PC + (1 | Subject)', 'FitMethod', 'REML');
results = compare(lme,altlme)
aov = anova(altlme,'DFMethod','Satterthwaite');

