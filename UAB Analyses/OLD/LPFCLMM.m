% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
%%
% Initialize a cell arrays to store power data for each subject
GroupPowerMatrixL = cell(numel(files), 1);
GroupDesignMatrixL = [];

for subj = 1:numel(files)

    load(files{subj})
    
    load(files{subj}, 'lpfcDesign','tfr','LPFC','frex','newt')

    tfr.powspctrm = tfr.powspctrm(:,ismember(tfr.label,LPFC),1:36,:);

    frex = frex(1:36);

    % rowsToremove = strcmpi(lpfcDesign.PreviousType,'None');
    shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);
    % rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

    % lpfcDesign(rowsToremove,:) = [];
    % tfr.powspctrm(rowsToremove_power,:,:,:) = [];
    % shortDesign(rowsToremove_power,:) = [];

    rowsToremove = strcmpi(lpfcDesign.CurrentType,'neu');% | strcmpi(lpfcDesign.PreviousType,'neu') | strcmpi(lpfcDesign.CurrentType,'inc');
    lpfcDesign(rowsToremove,:) = [];

    rowsToremove_power = strcmpi(shortDesign.CurrentType,'neu');% | strcmpi(shortDesign.PreviousType,'neu') | strcmpi(shortDesign.CurrentType,'inc');
    tfr.powspctrm(rowsToremove_power,:,:,:) = [];
    shortDesign(rowsToremove_power,:) = [];

    % Focus on only pI trials
    % rowsToremove = strcmpi(lpfcDesign.PreviousType,'con');
    % lpfcDesign(rowsToremove,:) = [];
    % 
    % rowsToremove_power = strcmpi(shortDesign.PreviousType,'con');
    % tfr.powspctrm(rowsToremove_power,:,:,:) = [];
    % shortDesign(rowsToremove_power,:) = [];

    GroupDesignMatrixL = cat(1,GroupDesignMatrixL,lpfcDesign);
    GroupPowerMatrixL{subj} = tfr.powspctrm;
end

GroupDesignMatrixL.RT = round(GroupDesignMatrixL.RT*1000);
GroupDesignMatrixL.PC = GroupDesignMatrixL.PC - mean(GroupDesignMatrixL.PC);
%%
GroupDesignMatrixL.CurrentType = categorical(GroupDesignMatrixL.CurrentType);
GroupDesignMatrixL.CurrentType = reordercats(GroupDesignMatrixL.CurrentType, {'con', 'inc'});
%%
GroupDesignMatrixL.PreviousType = categorical(GroupDesignMatrixL.PreviousType);
GroupDesignMatrixL.PreviousType = reordercats(GroupDesignMatrixL.PreviousType, {'con', 'inc'});
%%

[~,~,numFrequencies,numTimePoints] = size(tfr.powspctrm);

% Initialize matrices to store p-values and F-statistics for each frequency-time point

% Conflict
% conflictP = NaN(numFrequencies, numTimePoints);
% conflictB = NaN(numFrequencies, numTimePoints);
% conflictT = NaN(numFrequencies, numTimePoints);

adjustP = NaN(numFrequencies, numTimePoints);
adjustB = NaN(numFrequencies, numTimePoints);
adjustT = NaN(numFrequencies, numTimePoints);

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
        powerValues = [];
        
        for subj = 1:numel(GroupPowerMatrixL)
            % Number of electrodes for the current subject
            numElectrodes = size(GroupPowerMatrixL{subj}, 2);
            
            for elec = 1:numElectrodes
                % Extract the power value for the current electrode, frequency, and time
                powerValue = GroupPowerMatrixL{subj}(:,elec, freq, time);
                powerValues = [powerValues; powerValue];
            end
        end
        
        % Add these power values to the design matrix
        tempDesignMatrix = GroupDesignMatrixL;
        tempDesignMatrix.powerValues = powerValues;
        
        % Fit the LME model
        % lme = fitlme(tempDesignMatrix, 'powerValues ~ CurrentType + PreviousType + PC + CurrentType:PreviousType + CurrentType:PC + PreviousType:PC+ CurrentType:PreviousType:PC + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');
        % lme = fitlme(tempDesignMatrix, 'powerValues ~ CurrentType + PreviousType + CurrentType:PreviousType + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');
        % lme = fitlme(tempDesignMatrix, 'powerValues ~ CurrentType + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');

        % Adjustment Model
        lme = fitlme(tempDesignMatrix, 'powerValues ~ RT + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');

        % Extract coefficients and their p-values
        aov = anova(lme,'DFMethod','Satterthwaite');
        
        % Store the p-values and F-stats

        adjustP(freq,time) = aov.pValue(strcmpi(aov.Term,'RT'));
        adjustB(freq,time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'RT'));
        adjustT(freq,time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'RT'));

        % conflictP(freq, time) = aov.pValue(strcmpi(aov.Term,'CurrentType'));
        % conflictB(freq, time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
        % conflictT(freq, time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));

        % conrateP(freq, time) = aov.pValue(strcmpi(aov.Term,'PC'));
        % conrateB(freq, time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'PC'));

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

[~,fidx] = arrayfun(@(x) min(abs(x-frex)), [2 4 8 12 20 30 40]);
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
yticklabels([2 4 8 12 20 30 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.8*max(adjustT,[],'all'),3), round(0.8*max(adjustT,[],'all'),3)]);
lims = clim;
set(cb, 'YTick', [lims(1) 0 lims(2)])
cb.Label.String = ['\bf', 'Std Regression Coeff (t)', '\rm']; 
title('dlPFC Cue Aligned: Post-conflict Adjustment')
fontsize(gcf, 16, 'points')
%% Conflict main effect       
% [h,~,~,~] = fdr_bh(conflictP);
%% Keep p-value if apart of cluster that lasts >= 100 ms
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
yticklabels([2 4 8 12 20 30 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.9*max(conflictT,[],'all'),1) round(0.9*max(conflictT,[],'all'),1)]);
lims = clim;
set(cb,'YTick',[lims(1) 0 lims(2)])
cb.Label.String = ['\bf' 'Standardized Regression Coeff (t)' '\rm']; 
title('dlPFC: Conflict Main Effect')
fontsize(gcf,16,'points')

%% Proportion Congruent Effect
[h,~,~,~] = fdr_bh(conrateP);

contourf(newt, frex, conrateB, 50, 'LineColor', 'none')
hold on
contour(newt, frex, h, 1, 'linecolor', 'k', 'linewidth', 1.0);
set(gca, 'YScale', 'log')
yticks(frexticks)
yticklabels([2, 4, 8, 12, 20, 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.9*max(conrateB,[],'all'),1), round(0.9*max(conrateB,[],'all'),1)]);
lims = clim;
set(cb, 'YTick', [lims(1) 0 lims(2)])
cb.Label.String = ['\bf', '\beta', '\rm']; 
title('LPFC: PC Main Effect')
fontsize(gcf, 16, 'points')
%% Previous Effect
[h,~,~,~] = fdr_bh(previousP);

contourf(newt, frex, previousB, 50, 'LineColor', 'none')
hold on
contour(newt, frex, h, 1, 'linecolor', 'k', 'linewidth', 1.0);
set(gca, 'YScale', 'log')
yticks(frexticks)
yticklabels([2, 4, 8, 12, 20, 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.9*max(previousB,[],'all'),1), round(0.9*max(previousB,[],'all'),1)]);
lims = clim;
set(cb, 'YTick', [lims(1) 0 lims(2)])
cb.Label.String = ['\bf', '\beta', '\rm']; 
title('LPFC: Previous Main Effect')
fontsize(gcf, 16, 'points')
%% Sequence Effect
[h,~,~,~] = fdr_bh(sequenceP);

contourf(newt, frex, sequenceB, 50, 'LineColor', 'none')
hold on
contour(newt, frex, h, 1, 'linecolor', 'k', 'linewidth', 1.0);
set(gca, 'YScale', 'log')
yticks(frexticks)
yticklabels([2, 4, 8, 12, 20, 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.9*max(sequenceB,[],'all'),1), round(0.9*max(sequenceB,[],'all'),1)]);
lims = clim;
set(cb, 'YTick', [lims(1) 0 lims(2)])
cb.Label.String = ['\bf', '\beta', '\rm']; 
title('LPFC: Conflict x Previous Interaction')
fontsize(gcf, 16, 'points')

%% Conflict Effect Modulated by Proportion Congruent
[h,~,~,~] = fdr_bh(conXrateP);

contourf(newt, frex, conXrateB, 50, 'LineColor', 'none')
hold on
contour(newt, frex, h, 1, 'linecolor', 'k', 'linewidth', 1.0);
set(gca, 'YScale', 'log')
yticks(frexticks)
yticklabels([2, 4, 8, 12, 20, 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.9*max(conXrateB,[],'all'),1), round(0.9*max(conXrateB,[],'all'),1)]);
lims = clim;
set(cb, 'YTick', [lims(1) 0 lims(2)])
cb.Label.String = ['\bf', '\beta', '\rm']; 
title('LPFC: Conflict x PC')
fontsize(gcf, 16, 'points')
%% Sequence Modulated by Proportion Congruent
[h,~,~,~] = fdr_bh(seqXrateP);

contourf(newt, frex, seqXrateB, 50, 'LineColor', 'none')
hold on
contour(newt, frex, h, 1, 'linecolor', 'k', 'linewidth', 1.0);
set(gca, 'YScale', 'log')
yticks(frexticks)
yticklabels([2, 4, 8, 12, 20, 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([round(-0.9*max(seqXrateB,[],'all'),1), round(0.9*max(seqXrateB,[],'all'),1)]);
lims = clim;
set(cb, 'YTick', [lims(1) 0 lims(2)])
cb.Label.String = ['\bf', '\beta', '\rm']; 
title('LPFC: Sequence x PC')
fontsize(gcf, 16, 'points')
%% HFA Stats

% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];
proc_id = 'main_ft';

% Initialize a cell arrays to store power data for each subject
GroupPowerMatrixL = cell(numel(files), 1);
GroupDesignMatrixL = [];
%
for subj = 1:numel(files)

    load(files{subj})
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);

    dlpfcChans = elec.label(ismember(elec.ROI,'dlPFC'));
    nChans = numel(dlpfcChans);

    dlpfc_idx = ismember(hfa.label,dlpfcChans);
    

    [lpfcDesign,~] = fn_create_LMM_design(trial_info,nChans,1,'power');

    rowsToremove = strcmpi(lpfcDesign.PreviousType,'None');
    shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);
    rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

    errors = logical(trial_info.error);
    hfa.ztrial(errors,:,:) = [];

    lpfcDesign(rowsToremove,:) = [];
    hfa.ztrial(rowsToremove_power,:,:) = [];
    shortDesign(rowsToremove_power,:) = [];

    rowsToremove = strcmpi(lpfcDesign.CurrentType,'neu') | strcmpi(lpfcDesign.PreviousType,'neu') | strcmpi(lpfcDesign.CurrentType,'inc');
    lpfcDesign(rowsToremove,:) = [];

    rowsToremove_power = strcmpi(shortDesign.CurrentType,'neu') | strcmpi(shortDesign.PreviousType,'neu') | strcmpi(shortDesign.CurrentType,'inc');
    hfa.ztrial(rowsToremove_power,:,:) = [];
    shortDesign(rowsToremove_power,:) = [];

    hfa.ztrial(:,~dlpfc_idx,:) = [];

    
    % Focus on only pI trials
    rowsToremove = strcmpi(lpfcDesign.PreviousType,'con');
    lpfcDesign(rowsToremove,:) = [];

    rowsToremove_power = strcmpi(shortDesign.PreviousType,'con');
    hfa.ztrial(rowsToremove_power,:,:) = [];
    shortDesign(rowsToremove_power,:) = [];

    lpfcDesign.Subject = repmat(subj,size(lpfcDesign,1),1);
    GroupDesignMatrixL = cat(1,GroupDesignMatrixL,lpfcDesign);
    GroupPowerMatrixL{subj} = hfa.ztrial;
end
%%
GroupDesignMatrixL.RT = round(GroupDesignMatrixL.RT*1000);
GroupDesignMatrixL.PC = GroupDesignMatrixL.PC - mean(GroupDesignMatrixL.PC);
%%
GroupDesignMatrixL.CurrentType = categorical(GroupDesignMatrixL.CurrentType);
GroupDesignMatrixL.CurrentType = reordercats(GroupDesignMatrixL.CurrentType, {'con', 'inc'});
%%
GroupDesignMatrixL.PreviousType = categorical(GroupDesignMatrixL.PreviousType);
GroupDesignMatrixL.PreviousType = reordercats(GroupDesignMatrixL.PreviousType, {'con', 'inc'});
%%
[~,~,numTimePoints] = size(GroupPowerMatrixL{1});

% conflictP = NaN(numTimePoints,1);
% conflictB = NaN(numTimePoints,1);
% conflictT = NaN(numTimePoints,1);

adjustP = NaN(numTimePoints,1);
adjustB = NaN(numTimePoints,1);
adjustT = NaN(numTimePoints,1);

% previousP = NaN(numTimePoints,1);
% previousB = NaN(numTimePoints,1);
% previousT = NaN(numTimePoints,1);

for time = 1:numTimePoints
    % Initialize a column for power values
    powerValues = [];
    
    for subj = 1:numel(GroupPowerMatrixL)
        % Number of electrodes for the current subject
        numElectrodes = size(GroupPowerMatrixL{subj}, 2);
        
        for elec = 1:numElectrodes
            % Extract the power value for the current electrode, frequency, and time
            powerValue = GroupPowerMatrixL{subj}(:,elec, time);
            powerValues = [powerValues; powerValue];
        end
    end
    
    % Add these power values to the design matrix
    tempDesignMatrix = GroupDesignMatrixL;
    tempDesignMatrix.powerValues = powerValues;
    
    % Fit the LME model
    % lme = fitlme(tempDesignMatrix, 'powerValues ~ CurrentType + PreviousType + PC + CurrentType:PreviousType + CurrentType:PC + PreviousType:PC+ CurrentType:PreviousType:PC + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');
    % lme = fitlme(tempDesignMatrix, 'powerValues ~ CurrentType + PreviousType + CurrentType:PreviousType + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');
    % lme = fitlme(tempDesignMatrix, 'powerValues ~ CurrentType + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');

    % Previous type effect model
    % lme = fitlme(tempDesignMatrix, 'powerValues ~ PreviousType + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');


    % Adjustment Model
    lme = fitlme(tempDesignMatrix, 'powerValues ~ RT + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');

    % Extract coefficients and their p-values
    aov = anova(lme,'DFMethod','Satterthwaite');
    
    % Store the p-values and F-stats

    adjustP(time) = aov.pValue(strcmpi(aov.Term,'RT'));
    adjustB(time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'RT'));
    adjustT(time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'RT'));

    % conflictP(time) = aov.pValue(strcmpi(aov.Term,'CurrentType'));
    % conflictB(time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));
    % conflictT(time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'CurrentType_inc'));

    % previousP(time) = aov.pValue(strcmpi(aov.Term,'PreviousType'));
    % previousB(time) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'PreviousType_inc'));
    % previousT(time) = lme.Coefficients.tStat(strcmpi(lme.Coefficients.Name,'PreviousType_inc'));

   
end
%% Plot conflict stats
%[h,~,~,~] = fdr_bh(conflictP);
%% Keep p-value if apart of cluster that lasts >= 100 ms

% Shrink length of trial based on RT (S)
% conflictT = conflictT(31:151);
% conflictB = conflictB(31:151);
% conflictP = conflictP(31:151);
% hfa.newtime = hfa.newtime(31:151);

% Shrink length of trial based on RT (R)
conflictT = conflictT(26:176);
conflictB = conflictB(26:176);
conflictP = conflictP(26:176);
hfa.newtime = hfa.newtime(26:176);

threshP = conflictP < 0.05;
islands = bwconncomp(threshP);

sig_clusters = find(cellfun(@(x) length(x),islands.PixelIdxList) > 10);

for ii = 1:length(sig_clusters)
    sig_idxs{ii} = islands.PixelIdxList{sig_clusters(ii)};
end
%%
trialsidx = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);

con = [];
inc = [];

for i = 1:numel(GroupPowerMatrixL)
    temptrials = trialsidx([trialsidx.Subject] == i,:);
    con(:,i) = squeeze(mean(mean(GroupPowerMatrixL{i}([temptrials.CurrentType] == 'con',:,:)),2));
    inc(:,i) = squeeze(mean(mean(GroupPowerMatrixL{i}([temptrials.CurrentType] == 'inc',:,:)),2));
end

[hl, ~] = boundedline(hfa.newtime,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',hfa.newtime,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-0.2 2])
ylabel('HFA (z)')
xlabel('Time (s)')
sigtimes = NaN(length(con),1);
sigtimes(sig_idxs) = -0.9;
hold on
plot(hfa.newtime(31:end),sigtimes,'-k','LineWidth',2.5)
hold off

set(gca,'FontSize',16)

% Different way of plotting: plot the t-stats

plot(hfa.newtime,conflictT,'r')
hold on
for ii = 1:length(sig_idxs)
    plot(hfa.newtime(sig_idxs{ii}),conflictT(sig_idxs{ii}),'r','LineWidth',3)
end
hold off

title('dlPFC HFA Conflict Main Effect')
ylabel('Std Regression Coeff (t)')
xlabel('Time (s)')
RTs = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);
meanRT = mean(RTs.RT)/1000;
% xline(meanRT,'--k')
xline(0,'--k')
xlim([-0.75 0.75])
xticks([-0.75 -0.5 -0.25 0 0.25 0.5 0.75])
ylim([-3 4.5])
set(gca,'FontSize',16)

%% Plot previous effect on congruent trials
% h,~,~,~] = fdr_bh(conflictP);

% Shrink length of trial based on RT (S)
% previousT = previousT(31:151);
% previousB = previousB(31:151);
% previousP = previousP(31:151);
% hfa.newtime = hfa.newtime(31:151);

% Shrink length of trial based on RT (R)
previousT = previousT(26:176);
previousB = previousB(26:176);
previousP = previousP(26:176);
hfa.newtime = hfa.newtime(26:176);

threshP = previousP < 0.05;
islands = bwconncomp(threshP);

sig_clusters = find(cellfun(@(x) length(x),islands.PixelIdxList) > 10);

for ii = 1:length(sig_clusters)
    sig_idxs{ii} = islands.PixelIdxList{sig_clusters(ii)};
end

cC = [];
iC = [];

sbjs = unique(trialsidx.Subject);
for i = 1:numel(GroupPowerMatrixL)
    temptrials = trialsidx([trialsidx.Subject] == sbjs(i),:);
    cC(:,i) = squeeze(mean(mean(GroupPowerMatrixL{i}([temptrials.PreviousType] == 'con',:,:)),2));
    iC(:,i) = squeeze(mean(mean(GroupPowerMatrixL{i}([temptrials.PreviousType] == 'inc',:,:)),2));
end

[hl, ~] = boundedline(hfa.newtime,mean(cC,2),std(cC,[],2)./sqrt(size(cC,2)),'b',hfa.newtime,mean(iC,2),std(iC,[],2)./sqrt(size(iC,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
xlim([-0.2 2])
ylabel('HFA (z)')
xlabel('Time (s)')

plot(hfa.newtime,previousT,'r')
hold on
for ii = 1:length(sig_idxs)
    plot(hfa.newtime(sig_idxs{ii}),previousT(sig_idxs{ii}),'r','LineWidth',3)
end
hold off

title('dlPFC HFA Previous Main Effect')
ylabel('Std Regression Coeff (t)')
xlabel('Time (s)')
RTs = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);
meanRT = mean(RTs.RT)/1000;
% xline(meanRT,'--k')
xline(0,'--k')
xlim([-0.75 0.75])
xticks([-0.75 -0.5 -0.25 0 0.25 0.5 0.75])
set(gca,'FontSize',16)
%% Adjustment iC

% Shrink length of trial based on RT (S)
adjustT = adjustT(31:151);
adjustB = adjustB(31:151);
adjustP = adjustP(31:151);
hfa.newtime = hfa.newtime(31:151);

% Shrink length of trial based on RT (R)
% adjustT = adjustT(26:176);
% adjustB = adjustB(26:176);
% adjustP = adjustP(26:176);
% hfa.newtime = hfa.newtime(26:176);

threshP = adjustP < 0.05;
islands = bwconncomp(threshP);

sig_clusters = find(cellfun(@(x) length(x),islands.PixelIdxList) > 10);

for ii = 1:length(sig_clusters)
    sig_idxs{ii} = islands.PixelIdxList{sig_clusters(ii)};
end

% Different way of plotting: plot the t-stats
plot(hfa.newtime,adjustT,'r')
hold on
for ii = 1:length(sig_idxs)
    plot(hfa.newtime(sig_idxs{ii}),adjustT(sig_idxs{ii}),'r','LineWidth',3)
end
hold off

title('dlPFC HFA Post-conflict Adjustment')
ylabel('Std Regression Coeff (t)')
xlabel('Time (s)')
RTs = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);
meanRT = mean(RTs.RT)/1000;
xline(meanRT,'--k')
% xline(0,'--k')
xlim([-0.75 0.75])
xticks([-0.75 -0.5 -0.25 0 0.25 0.5 0.75])
set(gca,'FontSize',16)