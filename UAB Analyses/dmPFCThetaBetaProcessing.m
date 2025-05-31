restoredefaultpath;
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];
proc_id = 'main_ft';

addpath(ft_dir);
ft_defaults;

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);
%%
% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

% Initialize a cell arrays to store power data for each subject
GroupPowerMatrixM = cell(numel(files), 1);
GroupDesignMatrixM = [];

total_active = 0;
cnt = 1;

for subj = 1:numel(files)

    load(files{subj}, 'mpfcDesign','tfr','MPFC','frex','newt')
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    daccChans = intersect(elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC')),w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(daccChans)
        continue
    else

        nChans = numel(daccChans);

        total_active = total_active + nChans;

        % theta_power = squeeze(mean(tfr.powspctrm(:,ismember(tfr.label,daccChans),frex >= 4 & frex <= 8,:),3));
        theta_power = squeeze(mean(tfr.powspctrm(:,ismember(tfr.label,daccChans),frex >= 12 & frex <= 30,:),3)); %actually beta, didn't want to rename

        mpfcDesign([mpfcDesign.Electrode] > numel(daccChans),:) = [];
    
        rowsToremove = strcmpi(mpfcDesign.PreviousType,'None');
        shortDesign = mpfcDesign(mpfcDesign.Electrode == 1,:);
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

        mpfcDesign(rowsToremove,:) = [];
        theta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power ,:) = [];

        % rowsToremove = ~(strcmpi(mpfcDesign.PreviousType,'neu') | strcmpi(mpfcDesign.PreviousType,'con')) | strcmpi(mpfcDesign.CurrentType,'inc');
        % mpfcDesign(rowsToremove,:) = [];
        % 
        % rowsToremove_power = ~(strcmpi(shortDesign.PreviousType,'neu') | strcmpi(shortDesign.PreviousType,'con')) | strcmpi(shortDesign.CurrentType,'inc');
        % theta_power(rowsToremove_power,:,:) = [];
        % shortDesign(rowsToremove_power,:) = [];

        % For ITI Stim
        % rowsToremove = ~strcmpi(mpfcDesign.PreviousType,'neu');
        % mpfcDesign(rowsToremove,:) = [];
        % 
        % rowsToremove_power = ~strcmpi(shortDesign.PreviousType,'neu');
        % theta_power(rowsToremove_power,:,:) = [];
        % shortDesign(rowsToremove_power,:) = [];

        % Only 16% Con & Current Inc
        % rowsToremove = ~strcmpi(mpfcDesign.BlockType,'minc') | ~strcmpi(mpfcDesign.CurrentType,'inc');
        % mpfcDesign(rowsToremove,:) = [];
        % 
        % rowsToremove_power = ~strcmpi(shortDesign.BlockType,'minc') | ~strcmpi(shortDesign.CurrentType,'inc');
        % theta_power(rowsToremove_power,:,:) = [];
        % shortDesign(rowsToremove_power,:) = [];

        shortDesign.Subject = repmat(subj,height(shortDesign),1);

        GroupDesignMatrixM = cat(1,GroupDesignMatrixM,shortDesign);

        if numel(size(theta_power)) == 2
            GroupPowerMatrixM{cnt} = theta_power;
           
        else
            GroupPowerMatrixM{cnt} = squeeze(mean(theta_power,2));

        end
        cnt = cnt + 1;
    end
end

GroupPowerMatrixM = GroupPowerMatrixM(1:cnt-1);

GroupDesignMatrixM.RT = round(GroupDesignMatrixM.RT*1000);
GroupDesignMatrixM.PC = GroupDesignMatrixM.PC - 0.33;
% Create CurrentConflict column
GroupDesignMatrixM.CurrentConflict = repmat({'NoConflict'}, height(GroupDesignMatrixM), 1); % Default to 'NoConflict'
GroupDesignMatrixM.CurrentConflict(strcmp(GroupDesignMatrixM.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc'

% Create PreviousConflict column
GroupDesignMatrixM.PreviousConflict = repmat({'NoConflict'}, height(GroupDesignMatrixM), 1); % Default to 'NoConflict'
GroupDesignMatrixM.PreviousConflict(strcmp(GroupDesignMatrixM.PreviousType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where PreviousType is 'inc'


GroupDesignMatrixM.CurrentConflict = categorical(GroupDesignMatrixM.CurrentConflict);
GroupDesignMatrixM.CurrentConflict = reordercats(GroupDesignMatrixM.CurrentConflict, {'NoConflict', 'Conflict'});

GroupDesignMatrixM.PreviousConflict = categorical(GroupDesignMatrixM.PreviousConflict);
GroupDesignMatrixM.PreviousConflict = reordercats(GroupDesignMatrixM.PreviousConflict, {'NoConflict', 'Conflict'});

GroupDesignMatrixM.BlockType = categorical(GroupDesignMatrixM.BlockType);
GroupDesignMatrixM.BlockType = reordercats(GroupDesignMatrixM.BlockType, {'minc', 'same', 'mcon'});

RTs = GroupDesignMatrixM([GroupDesignMatrixM.Electrode] == 1,:);
meanRT = mean(RTs.RT)/1000;

%% TFR Plot

% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

% Initialize a cell arrays to store power data for each subject
GroupPowerMatrixM = cell(numel(files), 1);
GroupDesignMatrixM = [];
cnt = 1;
for subj = 1:numel(files)

    load(files{subj}, 'mpfcDesign','tfr','MPFC','frex','newt')
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    daccChans = intersect(elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC')),w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(daccChans)
        continue
    else
        power = squeeze(mean(tfr.powspctrm(:,ismember(tfr.label,daccChans),:,:),2));

        mpfcDesign([mpfcDesign.Electrode] > numel(daccChans),:) = [];
    
        shortDesign = mpfcDesign(mpfcDesign.Electrode == 1,:);

        shortDesign.Subject = repmat(subj,height(shortDesign),1);

        GroupDesignMatrixM = cat(1,GroupDesignMatrixM,shortDesign);

        GroupPowerMatrixM{cnt} = power;
       
        cnt = cnt + 1;
    end
end

GroupPowerMatrixM = GroupPowerMatrixM(1:cnt-1);

RTs = GroupDesignMatrixM([GroupDesignMatrixM.Electrode] == 1,:);
meanRT = mean(RTs.RT);

for i = 1:length(GroupPowerMatrixM)
    comp_power(:,:,i) = squeeze(mean(GroupPowerMatrixM{i}));
end

contourf(newt,frex,mean(comp_power,3),50,'EdgeColor','none')
set(gca,'YScale','log')
yticks([2 4 8 12 30 75 150])
% xticks([-0.5 0 0.5]) % R-locked
xticks([0 0.6 1.2]) % S-locked

cb = colorbar;
clim([-0.6 0.6])
set(cb,'yTick',[-0.6 0 0.6])
cb.Label.String = ['\bf' 'Power (z)' '\rm'];
% xline(0,'--k','LineWidth',1.25) % R-locked
xline(meanRT,'--k','LineWidth',1.25) % S-locked

title('dmPFC: Cue Aligned Grand Average')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% Median Split by RT
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

GroupPowerMatrixM_low = cell(numel(files), 1);
GroupPowerMatrixM_high = cell(numel(files), 1);
GroupDesignMatrixM_low = [];
GroupDesignMatrixM_high = [];

cnt = 1;
for subj = 1:numel(files)

    load(files{subj}, 'mpfcDesign','tfr','MPFC','frex','newt')
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    daccChans = intersect(...
        elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC')), ...
        w2{1,1}.label(w2{1,1}.sig_chans) ...
    );

    if isempty(daccChans)
        continue
    else
        % Extract theta power (or whichever band you prefer)
        theta_power = squeeze(mean(tfr.powspctrm(:, ...
                                 ismember(tfr.label,daccChans), ...
                                 frex >= 4 & frex <= 8, :),3));

        mpfcDesign([mpfcDesign.Electrode] > numel(daccChans), :) = [];
    
        rowsToremove = strcmpi(mpfcDesign.PreviousType,'None');
        shortDesign = mpfcDesign(mpfcDesign.Electrode == 1,:);
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

        mpfcDesign(rowsToremove,:) = [];
        theta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];


        rowsToremove = (strcmpi(mpfcDesign.PreviousType,'neu') | ...
                         strcmpi(mpfcDesign.PreviousType,'con')) | ...
                          strcmpi(mpfcDesign.CurrentType,'inc');
        mpfcDesign(rowsToremove,:) = [];

        rowsToremove_power = (strcmpi(shortDesign.PreviousType,'neu') | ...
                               strcmpi(shortDesign.PreviousType,'con')) | ...
                                strcmpi(shortDesign.CurrentType,'inc');
        theta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];

        shortDesign.Subject = repmat(subj,height(shortDesign),1);

        medianRT = median(shortDesign.RT);
        lowInds  = shortDesign.RT <= medianRT;
        highInds = shortDesign.RT >  medianRT;

        % Subset shortDesign
        shortDesign_low = shortDesign(lowInds, :);
        shortDesign_high = shortDesign(highInds, :);

        % Subset theta_power
        theta_power_low = theta_power(lowInds, :, :);
        theta_power_high = theta_power(highInds, :, :);

        % Concatenate into group-level arrays
        GroupDesignMatrixM_low = cat(1, GroupDesignMatrixM_low, shortDesign_low);
        GroupDesignMatrixM_high = cat(1, GroupDesignMatrixM_high, shortDesign_high);

        % Depending on how your power data are shaped, do the same logic
        if numel(size(theta_power)) == 2
            GroupPowerMatrixM_low{cnt} = theta_power_low;
            GroupPowerMatrixM_high{cnt} = theta_power_high;
        else
            GroupPowerMatrixM_low{cnt}  = squeeze(mean(theta_power_low, 2));
            GroupPowerMatrixM_high{cnt} = squeeze(mean(theta_power_high,2));
        end
   
    cnt = cnt + 1;
    end
end

GroupPowerMatrixM_high = GroupPowerMatrixM_high(1:cnt-1);
GroupPowerMatrixM_low = GroupPowerMatrixM_low(1:cnt-1);

