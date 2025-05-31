%% Simple Conflict effect overall
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
%%
meanRT_C = [];
meanRT_I = [];

for subj = 1:numel(files)
    load(files{subj}, "lpfcDesign");

    lpfcDesign(lpfcDesign.Electrode ~= 1,:) = [];
    % Remove rows from design matrix that have previous type "None"
    % rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'None');
    % lpfcDesign(rowsToRemove, :) = [];

    % Calculate mean RT for congruent (C) trials
    C_index = strcmpi(lpfcDesign.CurrentType, 'con') | strcmpi(lpfcDesign.CurrentType, 'neu');
    meanRT_C(subj) = round(mean(lpfcDesign.RT(C_index))*1000);

    % Calculate mean RT for incongruent (I) trials
    I_index = strcmpi(lpfcDesign.CurrentType, 'inc');
    meanRT_I(subj) = round(mean(lpfcDesign.RT(I_index))*1000);
end

% Calculate mean and SEM for each condition
mean_C = mean(meanRT_C);
mean_I = mean(meanRT_I);
sem_C = std(meanRT_C) / sqrt(length(meanRT_C));
sem_I = std(meanRT_I) / sqrt(length(meanRT_I));

% Create a figure for the plot
figure;

% Plot the means as bars
b = bar([1, 2], [mean_C, mean_I], 'FaceColor', 'flat', 'BarWidth', 0.5); % Adjust BarWidth for skinnier bars
b.CData(1,:) = [0.75, 0.75, 1]; % Blue congruent
b.CData(2,:) = [1, 0.75, 0.75]; % Red Incongruent
b.EdgeColor = 'black'; % Black edge color
b.LineWidth = 1.5; % Make edges thicker

hold on;

% Add error bars for SEM
errorbar([1, 2], [mean_C, mean_I], [sem_C, sem_I], 'k', 'LineStyle', 'none', 'LineWidth', 1.5);

% Overlay individual data points and connect them
for subj = 1:length(meanRT_C)
    plot([1.1, 1.9], [meanRT_C(subj), meanRT_I(subj)], '-k', 'LineWidth', 1); % Connect with black lines
end

% Add horizontal black line at the top from center of "C" bar to center of "I" bar
yline_position = max([meanRT_C, meanRT_I]) * 1.05; % Position slightly above the highest data point
plot([1, 2], [yline_position, yline_position], '-k', 'LineWidth', 1.5); % Horizontal line

% Add black star slightly above the center of the line
star_x = 1.5; % Center between "C" and "I"
star_y = yline_position + (0.01 * yline_position); % Slightly above the horizontal line
text(star_x, star_y, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontSize', 18); % Black star

% Tidy up plot
set(gca, 'XTick', [1 2], 'XTickLabel', {'NoConflict', 'Conflict'});
ylabel('Mean RT (ms)');
title('Group Level Stroop Effect');
ylim([min([meanRT_C, meanRT_I]) * 0.95, star_y * 1.05]); % Adjust y-limits slightly
%%

meanRT_CMC = [];
meanRT_IMC = [];
meanRT_CMI = [];
meanRT_IMI = [];
meanRT_CSA = [];
meanRT_ISA = [];

for subj = 1:numel(files)

    load(files{subj})

    % Remove rows from design matrix that have previous type "None"
    % rowsToRemove = strcmpi(lpfcDesign.PreviousType, 'None');
    % lpfcDesign(rowsToRemove, :) = [];

    % Remove rows from design matrix that have previous or current type "neu"
    % rowsToRemove =  strcmpi(lpfcDesign.CurrentType, 'neu') | strcmpi(lpfcDesign.PreviousType, 'neu');
    % lpfcDesign(rowsToRemove, :) = [];

    % Filter for 'mcon' BlockType
    mcon_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "mcon"), :);

    % Extract CMC and IMC pairs for 'mcon'
    CMC = mcon_trials(strcmpi(mcon_trials.CurrentType, "con") | strcmpi(mcon_trials.CurrentType, "neu"), :);
    IMC = mcon_trials(strcmpi(mcon_trials.CurrentType, "inc"), :);
    
    % Filter for 'minc' BlockType
    minc_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "minc"), :);

    % Extract CMI and IMI pairs for 'minc'
    CMI = minc_trials(strcmpi(minc_trials.CurrentType, "con") | strcmpi(minc_trials.CurrentType, "neu"), :);
    IMI = minc_trials(strcmpi(minc_trials.CurrentType, "inc"), :);

    % Filter for 'same' BlockType
    same_trials = lpfcDesign(strcmpi(lpfcDesign.BlockType, "same"), :);

    % Extract CMI and IMI pairs for 'minc'
    CSA = same_trials(strcmpi(same_trials.CurrentType, "con") | strcmpi(same_trials.CurrentType, "neu"), :);
    ISA = same_trials(strcmpi(same_trials.CurrentType, "inc"), :);

    meanRT_CMC(subj) = round(mean(CMC.RT * 1000));
    meanRT_IMC(subj) = round(mean(IMC.RT * 1000));
    meanRT_CMI(subj) = round(mean(CMI.RT * 1000));
    meanRT_IMI(subj) = round(mean(IMI.RT * 1000));
    meanRT_CSA(subj) = round(mean(CSA.RT * 1000));
    meanRT_ISA(subj) = round(mean(ISA.RT * 1000));
end

clearvars -except meanRT_CMC meanRT_IMC meanRT_CMI meanRT_IMI meanRT_ISA meanRT_CSA

% Calculate mean and SEM for mcon
mean_CMC = mean(meanRT_CMC);
mean_IMC = mean(meanRT_IMC);
sem_CMC = std(meanRT_CMC) / sqrt(length(meanRT_CMC));
sem_IMC = std(meanRT_IMC) / sqrt(length(meanRT_IMC));

% Calculate mean and SEM for minc
mean_CMI = mean(meanRT_CMI);
mean_IMI = mean(meanRT_IMI);
sem_CMI = std(meanRT_CMI) / sqrt(length(meanRT_CMI));
sem_IMI = std(meanRT_IMI) / sqrt(length(meanRT_IMI));

% Calculate mean and SEM for same
mean_CSA = mean(meanRT_CSA);
mean_ISA = mean(meanRT_ISA);
sem_CSA = std(meanRT_CSA) / sqrt(length(meanRT_CSA));
sem_ISA = std(meanRT_ISA) / sqrt(length(meanRT_ISA));

% Plot 
nsubjects = length(meanRT_ISA);

figure;
hold on
for i = 1:numel(meanRT_ISA)
    plot([1 2], [meanRT_IMI(i)-meanRT_CMI(i) meanRT_ISA(i)-meanRT_CSA(i)],'-k','Marker','o','MarkerSize',4,'MarkerFaceColor','k','LineWidth',0.5)
    plot([2 3], [meanRT_ISA(i)-meanRT_CSA(i) meanRT_IMC(i)-meanRT_CMC(i)],'-k','Marker','o','MarkerSize',4,'MarkerFaceColor','k','LineWidth',0.5)  
end

xlim([0 4])
ylim([-20 500])
yticks([0 125 250 375 500])
ylabel('Conflict-NoConflict RT (ms)')
xlabel('Proportion Congruent')
set(gca, 'XTick', [1 2 3], 'XTickLabel', {'16%', '33%', '50%'});

% Add horizontal black line at the top from center of "C" bar to center of "I" bar
yline_position = 475; % Position slightly above the highest data point
plot([1, 3], [yline_position, yline_position], '-k', 'LineWidth', 1.5); % Horizontal line

% Add black star slightly above the center of the line
star_x = 2; % Center between "C" and "I"
star_y = yline_position + (0.01 * yline_position); % Slightly above the horizontal line
text(star_x, star_y, '*', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', 'k', 'FontSize', 18); % Black star


plot([1 2],[nanmean(meanRT_IMI-meanRT_CMI) nanmean(meanRT_ISA-meanRT_CSA)],'r','LineWidth',1.5)
errorbar([1, 2], [nanmean(meanRT_IMI-meanRT_CMI) nanmean(meanRT_ISA-meanRT_CSA)], [nanstd(meanRT_IMI-meanRT_CMI)./sqrt(nsubjects), nanstd(meanRT_ISA-meanRT_CSA)./sqrt(nsubjects)], 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

plot([2 3],[mean(meanRT_ISA-meanRT_CSA) mean(meanRT_IMC-meanRT_CMC)],'r','LineWidth',1.5)
errorbar([2, 3], [mean(meanRT_ISA-meanRT_CSA) mean(meanRT_IMC-meanRT_CMC)], [std(meanRT_ISA-meanRT_CSA)./sqrt(nsubjects), std(meanRT_IMC-meanRT_CMC)./sqrt(nsubjects)], 'r', 'LineStyle', 'none', 'LineWidth', 1.5);

title('Stroop Effect by Block Type')

%% Unity Plot
conflictMC = meanRT_IMC - meanRT_CMC;
conflictMI = meanRT_IMI - meanRT_CMI;
scatter(conflictMI,conflictMC,190,'Marker','o','MarkerFaceColor',[0.4940 0.1840 0.5560],'MarkerEdgeColor',[0.4940 0.1840 0.5560])
xlim([0 500])
ylim([0 500])
hold on
plot([0 500],[0 500],'--k','LineWidth',1.5)
hold off
ylabel('50% PC')
xlabel('16% PC')
xticks([0 250 500])
yticks([0 250 500])
title('Stroop Effect Modulated by Conflict Expectation')
