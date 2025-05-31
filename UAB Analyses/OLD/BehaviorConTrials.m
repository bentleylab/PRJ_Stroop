
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

meanRT_cC = [];
meanRT_iC = [];

for subj = 1:numel(files)
    load(files{subj}, "lpfcDesign");

    lpfcDesign(lpfcDesign.Electrode ~= 1,:) = [];
    % Remove rows from design matrix that have previous type "None"
    rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'None');
    lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that have previous or current type "neu"
    % rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'neu') | strcmpi(lpfcDesign.CurrentType, 'neu');
    % lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that are current type "inc"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'inc');
    lpfcDesign(rowsToRemove, :) = [];
    
    % Calculate mean RT for congruent (cC) trials
    cC_index = strcmpi(lpfcDesign.PreviousType, 'con') | strcmpi(lpfcDesign.PreviousType, 'neu');
    meanRT_cC(subj) = round(mean(lpfcDesign.RT(cC_index)*1000));

    % Calculate mean RT for incongruent (iC) trials
    iC_index = strcmpi(lpfcDesign.PreviousType, 'inc');
    meanRT_iC(subj) = round(mean(lpfcDesign.RT(iC_index)*1000));
end

% Calculate mean and SEM for each condition
mean_cC = mean(meanRT_cC);
mean_iC = mean(meanRT_iC);
sem_cC = std(meanRT_cC) / sqrt(length(meanRT_cC));
sem_iC = std(meanRT_iC) / sqrt(length(meanRT_iC));

% Create a figure for the plot
figure;

% Plot the means as bars
b = bar([1, 2], [mean_cC, mean_iC], 'FaceColor', 'flat', 'BarWidth', 0.5); % Adjust BarWidth for skinnier bars
b.CData(1,:) = [0.75, 0.75, 1];% Blue color for cC
b.CData(2,:) = [1, 0.75, 0.75]; % Red color for iC
b.EdgeColor = 'black'; % Black edge color
b.LineWidth = 1.5; % Make edges thicker

hold on;

% Add error bars for SEM
errorbar([1, 2], [mean_cC, mean_iC], [sem_cC, sem_iC], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Overlay individual data points and connect them
for subj = 1:length(meanRT_cC)
    plot([1.1, 1.9], [meanRT_cC(subj), meanRT_iC(subj)], '-k', 'LineWidth', 1); % Connect with black lines
end

% Add horizontal black line at the top from center of "C" bar to center of "I" bar
yline_position = max([meanRT_cC, meanRT_iC]) * 1.05; % Position slightly above the highest data point
plot([1, 2], [yline_position, yline_position], '-k', 'LineWidth', 1.5); % Horizontal line

% Add black star slightly above the center of the line
star_x = 1.5; % Center between "C" and "I"
star_y = yline_position + (0.01 * yline_position); % Slightly above the horizontal line
text(star_x, star_y, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontSize', 18); % Black star

% Enhance plot
set(gca, 'XTick', [1 2], 'XTickLabel', {'Previous NoConflict', 'Previous Conflict'});
ylabel('Mean RT (ms)');
title('Group Level Gratton Effect: Current NoConflict');
ylim([min([meanRT_cC, meanRT_iC]) * 0.95, star_y * 1.05]); % Adjust y-limits slightly
yticks([600 750 900 1050 1200])
%%
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];

files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
files = cellfun(@(x) strrep(x, '.mat', ''), files, 'UniformOutput', false);

meanRT_mcon_cC = [];
meanRT_mcon_iC = [];
meanRT_minc_cC = [];
meanRT_minc_iC = [];

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

    % Remove rows from design matrix that are current type "inc"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'inc');
    lpfcDesign(rowsToRemove, :) = [];

    % Filter for 'mcon' BlockType
    mcon_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "mcon"), :);

    % Extract cI and iI pairs for 'mcon'
    mcon_cC = mcon_trials(strcmpi(mcon_trials.PreviousType, "con"), :);
    mcon_iC = mcon_trials(strcmpi(mcon_trials.PreviousType, "inc"), :);
    
    % Filter for 'minc' BlockType
    minc_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "minc"), :);
    % Extract cI and iI pairs for 'minc'
    minc_cC = minc_trials(strcmpi(minc_trials.PreviousType, "con"), :);
    minc_iC = minc_trials(strcmpi(minc_trials.PreviousType, "inc"), :);

    if isempty(minc_cC)
        continue
    end
    meanRT_mcon_cC(subj) = round(mean(mcon_cC.RT * 1000));
    meanRT_mcon_iC(subj) = round(mean(mcon_iC.RT * 1000));
    meanRT_minc_cC(subj) = round(mean(minc_cC.RT * 1000));
    meanRT_minc_iC(subj) = round(mean(minc_iC.RT * 1000));
end

clearvars -except meanRT_minc_iC meanRT_minc_cC meanRT_mcon_iC meanRT_mcon_cC

meanRT_mcon_cC([5 7]) = [];
meanRT_mcon_iC([5 7]) = [];
meanRT_minc_cC([5 7]) = [];
meanRT_minc_iC([5 7]) = [];

% Calculate mean and SEM for mcon
mean_mcon_cC = mean(meanRT_mcon_cC);
mean_mcon_iC = mean(meanRT_mcon_iC);
sem_mcon_cC = std(meanRT_mcon_cC) / sqrt(length(meanRT_mcon_cC));
sem_mcon_iC = std(meanRT_mcon_iC) / sqrt(length(meanRT_mcon_iC));

% Calculate mean and SEM for minc
mean_minc_cC = mean(meanRT_minc_cC);
mean_minc_iC = mean(meanRT_minc_iC);
sem_minc_cC = std(meanRT_minc_cC) / sqrt(length(meanRT_minc_cC));
sem_minc_iC = std(meanRT_minc_iC) / sqrt(length(meanRT_minc_iC));

% Plot for mcon
figure;
b_mcon = bar([1, 2], [mean_mcon_cC, mean_mcon_iC], 'FaceColor', 'flat', 'BarWidth', 0.5);
b_mcon.CData(1,:) = [1, 0.75, 0.75]; % Pink for cC
b_mcon.CData(2,:) = [0.75, 0.75, 1]; % Blue for iC
b_mcon.EdgeColor = 'black';
b_mcon.LineWidth = 1.5;
hold on;
errorbar([1, 2], [mean_mcon_cC, mean_mcon_iC], [sem_mcon_cC, sem_mcon_iC], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
for subj = 1:length(meanRT_mcon_cC)
    plot([1.1, 1.9], [meanRT_mcon_cC(subj), meanRT_mcon_iC(subj)], '-k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'cC', 'iC'});
ylabel('Mean Reaction Time (ms)');
title('Mostly Congruent: cC vs iC');
ylim([min([meanRT_mcon_cC, meanRT_mcon_iC]) * 0.95, max([meanRT_mcon_cC, meanRT_mcon_iC]) * 1.05]);

% Plot for minc
figure;
b_minc = bar([1, 2], [mean_minc_cC, mean_minc_iC], 'FaceColor', 'flat', 'BarWidth', 0.5);
b_minc.CData(1,:) = [1, 0.75, 0.75]; % Pink for cC
b_minc.CData(2,:) = [0.75, 0.75, 1]; % Blue for iC
b_minc.EdgeColor = 'black';
b_minc.LineWidth = 1.5;
hold on;
errorbar([1, 2], [mean_minc_cC, mean_minc_iC], [sem_minc_cC, sem_minc_iC], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
for subj = 1:length(meanRT_minc_cC)
    plot([1.1, 1.9], [meanRT_minc_cC(subj), meanRT_minc_iC(subj)], '-k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'cC', 'iC'});
ylabel('Mean Reaction Time (ms)');
title('Mostly Incongruent: cC vs iC');
ylim([min([meanRT_minc_cC, meanRT_minc_iC]) * 0.95, max([meanRT_minc_cC, meanRT_minc_iC]) * 1.05]);

%% "Same" BlockType

root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];

files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
files = cellfun(@(x) strrep(x, '.mat', ''), files, 'UniformOutput', false);

meanRT_same_cC = [];
meanRT_same_iC = [];

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

    % Remove rows from design matrix that are current type "inc"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'inc');
    lpfcDesign(rowsToRemove, :) = [];

    % Filter for 'same' BlockType
    same_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "same"), :);

    % Extract cC and iC pairs for 'same'
    same_cC = same_trials(strcmpi(same_trials.PreviousType, "con"), :);
    same_iC = same_trials(strcmpi(same_trials.PreviousType, "inc"), :);

    % Calculate the mean reaction time for cC and iC in 'same'
    meanRT_same_cC(subj) = round(mean(same_cC.RT * 1000));
    meanRT_same_iC(subj) = round(mean(same_iC.RT * 1000));
end

clearvars -except meanRT_same_iC meanRT_same_cC

% Calculate mean and SEM for 'same'
mean_same_cC = mean(meanRT_same_cC);
mean_same_iC = mean(meanRT_same_iC);
sem_same_cC = std(meanRT_same_cC) / sqrt(length(meanRT_same_cC));
sem_same_iC = std(meanRT_same_iC) / sqrt(length(meanRT_same_iC));

% Plot for 'same'
figure;
b_same = bar([1, 2], [mean_same_cC, mean_same_iC], 'FaceColor', 'flat', 'BarWidth', 0.5);
b_same.CData(1,:) = [1, 0.75, 0.75]; % Pink for cC
b_same.CData(2,:) = [0.75, 0.75, 1]; % Blue for iC
b_same.EdgeColor = 'black';
b_same.LineWidth = 1.5;
hold on;
errorbar([1, 2], [mean_same_cC, mean_same_iC], [sem_same_cC, sem_same_iC], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
for subj = 1:length(meanRT_same_cC)
    plot([1.1, 1.9], [meanRT_same_cC(subj), meanRT_same_iC(subj)], '-k', 'LineWidth', 1);
end
set(gca, 'XTick', [1 2], 'XTickLabel', {'cC', 'iC'});
ylabel('Mean Reaction Time (ms)');
title('Equal Congruent: cC vs iC');
ylim([min([meanRT_same_cC, meanRT_same_iC]) * 0.95, max([meanRT_same_cC, meanRT_same_iC]) * 1.05]);
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

    % Remove rows from design matrix that are current type "inc"
    rowsToRemove = strcmpi(lpfcDesign.CurrentType, 'inc');
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

