% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];
proc_id = 'main_ft';

% Initialize a cell arrays to store power data for each subject
GroupPowerMatrixL = cell(numel(files), 1);
GroupDesignMatrixL = [];
cnt = 1;
for subj = 1:numel(files)

    load(files{subj}, 'lpfcDesign','tfr','LPFC','frex','newt')

    % load(['/Users/anaskhan/Desktop/PRJ_Stroop/results/Power/ThetaStim/' files{subj}])
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])
    % sigLPFC = intersect(LPFC,sigChans);

    % Conflict encoding
    dlpfcChans = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(dlpfcChans)
        continue
    else
        beta_power = squeeze(mean(tfr.powspctrm(:,ismember(tfr.label,dlpfcChans),frex >= 12 & frex <= 30,:),3));

        lpfcDesign([lpfcDesign.Electrode] > numel(dlpfcChans),:) = [];
    
        rowsToremove = strcmpi(lpfcDesign.PreviousType,'None');
        shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');
    
        lpfcDesign(rowsToremove,:) = [];
        beta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power ,:) = [];
    
        % Follow-up on only Inc trials
        rowsToremove = ~strcmpi(lpfcDesign.CurrentType,'inc');
        lpfcDesign(rowsToremove,:) = [];

        rowsToremove_power = ~strcmpi(shortDesign.CurrentType,'inc');
        beta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];

        % Only next trial congruent
        rowsToremove = ~strcmpi(lpfcDesign.PostType,'con');
        lpfcDesign(rowsToremove,:) = [];

        rowsToremove_power = ~strcmpi(shortDesign.PostType,'con');
        beta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];
    
        GroupDesignMatrixL = cat(1,GroupDesignMatrixL,lpfcDesign);
        GroupPowerMatrixL{cnt} = mean(beta_power(:,:,newt >= 0 & newt <= 0.55),3);
        cnt = cnt + 1;
    end
end

GroupPowerMatrixL = GroupPowerMatrixL(1:cnt-1);
%%
GroupDesignMatrixL.PostRT = round(GroupDesignMatrixL.PostRT*1000);
%%
GroupDesignMatrixL.CurrentType = categorical(GroupDesignMatrixL.CurrentType);
GroupDesignMatrixL.CurrentType = reordercats(GroupDesignMatrixL.CurrentType, {'con', 'inc'});
%%
GroupDesignMatrixL.PreviousType = categorical(GroupDesignMatrixL.PreviousType);
GroupDesignMatrixL.PreviousType = reordercats(GroupDesignMatrixL.PreviousType, {'con', 'inc'});
%% LMM

% Initialize a column for power values
powerValues = [];

for subj = 1:numel(GroupPowerMatrixL)

    % Number of electrodes for the current subject
    numElectrodes = size(GroupPowerMatrixL{subj}, 2);

    for elec = 1:numElectrodes
        % Extract beta power
        powerValue = GroupPowerMatrixL{subj}(:,elec);
       
        powerValues = [powerValues; powerValue];
    end
end

% Add these power values to the design matrix
tempDesignMatrix = GroupDesignMatrixL;
tempDesignMatrix.powerValues = powerValues;

% Fit the LME model
lme = fitlme(tempDesignMatrix, 'PostRT ~ powerValues + (1 | Subject) + (1 | Subject:Electrode)', 'FitMethod', 'REML');


% Extract coefficients and their p-values
aov = anova(lme,'DFMethod','Satterthwaite');

