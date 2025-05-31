restoredefaultpath;
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];

addpath(ft_dir);
ft_defaults;

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);

%%
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

SBJs = cellfun(@(x) x(1:end-4), files, 'UniformOutput', false);

allRTs = [];
for sbj = 1:numel(SBJs)
    SBJ = SBJs{sbj};
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    trials = fn_create_LMM_design(trial_info,1,1,'power');
    trials(strcmpi(trials.PreviousType,'None'),:) = [];
    trials.Subject = repmat(sbj,size(trials,1),1);
    allRTs = [allRTs;trials];
end

%%
% Create CurrentConflict column
allRTs.CurrentConflict = repmat({'NoConflict'}, height(allRTs), 1); % Default to 'NoConflict'
allRTs.CurrentConflict(strcmp(allRTs.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc'

% Create PreviousConflict column
allRTs.PreviousConflict = repmat({'NoConflict'}, height(allRTs), 1); % Default to 'NoConflict'
allRTs.PreviousConflict(strcmp(allRTs.PreviousType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where PreviousType is 'inc'


%%
allRTs.CurrentType = cellfun(@(x) [upper(x(1)), lower(x(2:end))], allRTs.CurrentType, 'UniformOutput', false);

allRTs.RT = round(allRTs.RT*1000);

allRTs.PC = allRTs.PC - 0.33;
%%
allRTs.CurrentConflict = categorical(allRTs.CurrentConflict);
allRTs.CurrentConflict = reordercats(allRTs.CurrentConflict, {'NoConflict', 'Conflict'});
%%
allRTs.PreviousConflict = categorical(allRTs.PreviousConflict);
allRTs.PreviousConflict = reordercats(allRTs.PreviousConflict, {'NoConflict', 'Conflict'});

allRTs.logRT = log(allRTs.RT);

lme = fitlme(allRTs,'logRT ~ CurrentConflict + PreviousConflict + PC + PreviousConflict:PC + CurrentConflict:PreviousConflict + CurrentConflict:PC + CurrentConflict:PreviousConflict:PC + (1 | Subject)','FitMethod', 'REML')

%% Export to R

% Write the table to a CSV file
writetable(allRTs, 'StroopRTData.csv');

%% RT Histograms
histogram(allRTs.RT(strcmpi(allRTs.CurrentConflict,'NoConflict')),'FaceColor','b','FaceAlpha',0.5,'Normalization','probability')
hold on
histogram(allRTs.RT(strcmpi(allRTs.CurrentConflict,'Conflict')),'FaceColor','r','FaceAlpha',0.5,'Normalization','probability')
ylabel('Probability')
xlabel('RT (ms)')
xticks([500 1000 1500 2000])
title('RT Distributions by Condition')
