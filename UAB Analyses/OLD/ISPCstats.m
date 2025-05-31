% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
%%
% Initialize a cell arrays to store power data for each subject
GroupISPC = cell(numel(files), 1);
GroupDesign = [];

for subj = 1:numel(files)

    load(files{subj},"coherenceDesign","MPFC","LPFC","zISPC","frex","newt")

    rowsToremove = strcmpi(coherenceDesign.PreviousType,'None');
    shortDesign = coherenceDesign(coherenceDesign.ElectrodePair == 1,:);
    rowsToremove_ispc = strcmpi(shortDesign.PreviousType,'None');

    coherenceDesign(rowsToremove,:) = [];
    zISPC(rowsToremove_ispc,:,:,:,:) = [];
    shortDesign(rowsToremove_ispc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.CurrentType,'neu') | strcmpi(coherenceDesign.PreviousType,'neu');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_ispc = strcmpi(shortDesign.CurrentType,'neu') | strcmpi(shortDesign.PreviousType,'neu');
    zISPC(rowsToremove_ispc,:,:,:,:) = [];
    shortDesign(rowsToremove_ispc,:) = [];

    rowsToremove = strcmpi(coherenceDesign.CurrentType,'inc');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_ispc = strcmpi(shortDesign.CurrentType,'inc');
    zISPC(rowsToremove_ispc,:,:,:,:) = [];
    shortDesign(rowsToremove_ispc,:) = [];

    % zISPC = mean(zISPC(:,:,:,frex >= 2 & frex <= 8,:),4);

    % Focus on only pI trials
    rowsToremove = strcmpi(coherenceDesign.PreviousType,'con');
    coherenceDesign(rowsToremove,:) = [];

    rowsToremove_ispc = strcmpi(shortDesign.PreviousType,'con');
    zISPC(rowsToremove_ispc,:,:,:,:) = [];
    shortDesign(rowsToremove_ispc,:) = [];

    GroupDesign = cat(1,GroupDesign,coherenceDesign);
    GroupISPC{subj} = zISPC;
end

GroupDesign.PC = GroupDesign.PC - mean(GroupDesign.PC);
GroupDesign.RT = round(GroupDesign.RT*1000);
%%
GroupDesign.CurrentType = categorical(GroupDesign.CurrentType);
GroupDesign.CurrentType = reordercats(GroupDesign.CurrentType, {'con', 'inc'});
%%
GroupDesign.PreviousType = categorical(GroupDesign.PreviousType);
GroupDesign.PreviousType = reordercats(GroupDesign.PreviousType, {'con', 'inc'});
%%

[~,~,~,numFrequencies,numTimePoints] = size(zISPC);

% Initialize matrices to store p-values and F-statistics for each frequency-time point

conflictP = NaN(numFrequencies,numTimePoints,1);
conflictB = NaN(numFrequencies,numTimePoints,1);
conflictT = NaN(numFrequencies,numTimePoints,1);

% Conflict
% conflictP = NaN(numTimePoints,1);
% conflictB = NaN(numTimePoints,1);

% % Proportion congruent
% conrateP = NaN(numFrequencies, numTimePoints);
% conrateB = NaN(numFrequencies, numTimePoints);
% 
% % Previous effect
% previousP = NaN(numFrequencies, numTimePoints);
% previousB = NaN(numFrequencies, numTimePoints);
% 
% % Sequence effect (conflict effect modulated by previous trial inc vs con)
% sequenceP = NaN(numFrequencies, numTimePoints);
% sequenceB = NaN(numFrequencies, numTimePoints);
% 
% % Conflict effect modulated by proportion congruent
% conXrateP = NaN(numFrequencies, numTimePoints);
% conXrateB = NaN(numFrequencies, numTimePoints);
% 
% % Sequence modulated by proportion congruent
% seqXrateP = NaN(numFrequencies, numTimePoints);
% seqXrateB = NaN(numFrequencies, numTimePoints);

for freq = 1:numFrequencies
    for time = 1:numTimePoints
        % Initialize a column for power values
        ispcs = [];
        
        for subj = 1:numel(GroupISPC)
            % Number of electrodes for the current subject
            n_MPFC_elecs = size(GroupISPC{subj}, 2);
            n_LPFC_elecs = size(GroupISPC{subj}, 3);
            for elecM = 1:n_MPFC_elecs
                for elecL = 1:n_LPFC_elecs
                    % Extract the power value for the current electrode, frequency, and time
                    ispc = GroupISPC{subj}(:,elecM,elecL, freq, time);
                    ispcs = [ispcs; ispc];
                end
            end
        end
    
        
            % Add these power values to the design matrix
            tempDesignMatrix = GroupDesign;
            tempDesignMatrix.ISPC = ispcs;
            
            % Fit the LME model
            % lme = fitlme(tempDesignMatrix, 'ISPC ~ CurrentType + PreviousType + PC + CurrentType*PreviousType + CurrentType*PC + CurrentType:PreviousType:PC + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');
            % lme = fitlme(tempDesignMatrix, 'ISPC ~ CurrentType + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');

            % Adjustment Model
            lme = fitlme(tempDesignMatrix, 'ISPC ~ RT + (1 | Subject) + (1 | Subject:ElectrodePair)', 'FitMethod', 'REML');
    
            % Extract coefficients and their p-values
            aov = anova(lme,'DFMethod','Satterthwaite');
            
            % Store the p-values and F-stats
    
            % conflictP(freq,time) = aov.pValue(strcmpi(aov.Term,'CurrentType'));
            % conflictB(freq,time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
            % conflictT(freq,time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));

            adjustP(freq,time) = aov.pValue(strcmpi(aov.Term,'RT'));
            adjustB(freq,time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'RT'));
            adjustT(freq,time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'RT'));

            % conflictP(time) = aov.pValue(strcmpi(aov.Term,'CurrentType'));
            % conflictB(time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
            % 
            % conrateP(freq, time) = aov.pValue(strcmpi(aov.Term,'PC'));
            % conrateB(freq, time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'PC'));
            % 
            % previousP(freq, time) = aov.pValue(strcmpi(aov.Term,'PreviousType'));
            % previousB(freq, time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'PreviousType_inc'));
    
            % sequenceP(freq, time) = aov.pValue(strcmpi(aov.Term,'CurrentType:PreviousType'));
            % sequenceB(freq, time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc:PreviousType_inc'));
            % 
            % conXrateP(freq, time) = aov.pValue(strcmpi(aov.Term,'CurrentType:PC'));
            % conXrateB(freq, time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc:PC'));
            % 
            % seqXrateP(freq, time) = aov.pValue(strcmpi(aov.Term,'CurrentType:PC:PreviousType'));
            % seqXrateB(freq, time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc:PC:PreviousType_inc'));
    end
end
%%
[~,fidx] = arrayfun(@(x) min(abs(x-frex)), [2 4 8 12 30 40]);
frexticks = frex(fidx);
clear fidx
%% iC Adjustment Effect

threshP = adjustP < 0.001;
islands = bwconncomp(threshP);

for ii = 1:islands.NumObjects
    [row_frex,col_time] = ind2sub([size(adjustP)],islands.PixelIdxList{ii});
    if ~(max(row_frex) - min(row_frex) >= 3 && max(col_time) - min(col_time) >= 13)
        threshP(islands.PixelIdxList{ii}) = 0;
    else
        threshP(islands.PixelIdxList{ii}) = 1;
    end
end

contourf(newt, frex, adjustT, 50, 'LineColor', 'none')
hold on
contour(newt, frex, threshP, 1, 'linecolor', 'k', 'linewidth', 1.0);
set(gca, 'YScale', 'log')
yticks(frexticks)
yticklabels([2 4 8 12 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.8*max(adjustB,[],'all'),3), round(0.8*max(adjustB,[],'all'),3)]);
lims = clim;
set(cb, 'YTick', [lims(1) 0 lims(2)])
cb.Label.String = ['\bf', '\beta', '\rm']; 
title('dlPFC Response Aligned: iC Power vs RT')
fontsize(gcf, 16, 'points')

%% Conflict effect
threshP = conflictP < 0.001;
islands = bwconncomp(threshP);

for ii = 1:islands.NumObjects
    [row_frex,col_time] = ind2sub([size(conflictP)],islands.PixelIdxList{ii});
    if ~(max(row_frex) - min(row_frex) >= 3 && max(col_time) - min(col_time) >= 13)
        threshP(islands.PixelIdxList{ii}) = 0;
    else
        threshP(islands.PixelIdxList{ii}) = 1;
    end
end

contourf(newt,frex,conflictT,50,'LineColor','none')
hold on
contour(newt,frex,threshP,1,'linecolor','k','linewidth',1.0);
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 30 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.9*max(conflictT,[],'all'),1) round(0.9*max(conflictT,[],'all'),1)]);
lims = clim;
set(cb,'YTick',[lims(1) 0 lims(2)])
cb.Label.String = ['\bf' 'Standardized Regression Coeff (t)' '\rm']; 
title('dlPFC: Conflict Main Effect')
fontsize(gcf,16,'points')
%% ISPC by current trial type

for i = 1:numel(files)
    
    load(files{i})
    coherenceDesign([coherenceDesign.ElectrodePair] ~= 1,:) = [];

    Con(:,:,i) = mean(mean(mean(zISPC(strcmpi(coherenceDesign.CurrentType,'con'),:,:,:,:)),2),3);
   
    Inc(:,:,i) = mean(mean(mean(zISPC(strcmpi(coherenceDesign.CurrentType,'inc'),:,:,:,:)),2),3);
end

contourf(newt,frex,mean(Con,3),50,'EdgeColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 30 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-1 1])

figure(2)
contourf(newt,frex,mean(Inc,3),50,'EdgeColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 30 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-1 1])

figure(3)
contourf(newt,frex,mean(Inc-Con,3),50,'EdgeColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-1 1])

%% ISPC grand average

for i = 1:numel(files)
    
    load(files{i})
    coherenceDesign([coherenceDesign.ElectrodePair] ~= 1,:) = [];

    Alltrials(:,:,i) = mean(mean(mean(zISPC(:,:,:,:,:)),2),3);
   
end

contourf(newt,frex,mean(Alltrials,3),50,'EdgeColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 30 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-1 1])
lims = clim;
set(cb,'YTick',[lims(1) 0 lims(2)])
cb.Label.String = ['\bf' 'ISPC (z)' '\rm']; 
fontsize(gcf,16,'points')
