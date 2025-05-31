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
GroupPowerMatrixL = cell(numel(files), 1);
GroupDesignMatrixL = [];
cnt = 1;
for subj = 1:numel(files)

    load(files{subj}, 'lpfcDesign','tfr','LPFC','frex','newt')
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Task encoding
    dlpfcChans = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(dlpfcChans)
        continue
    else
        % theta_power = squeeze(mean(tfr.powspctrm(:,ismember(tfr.label,dlpfcChans),frex >= 4 & frex <= 8,:),3));
        theta_power = squeeze(mean(tfr.powspctrm(:,ismember(tfr.label,dlpfcChans),frex >= 12 & frex <= 30,:),3)); %actually beta, didn't want to rename

        lpfcDesign([lpfcDesign.Electrode] > numel(dlpfcChans),:) = [];
    
        rowsToremove = strcmpi(lpfcDesign.PreviousType,'None');
        shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

        lpfcDesign(rowsToremove,:) = [];
        theta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power ,:) = [];

        % rowsToremove = ~(strcmpi(lpfcDesign.PreviousType,'neu') | strcmpi(lpfcDesign.PreviousType,'con')) | strcmpi(lpfcDesign.CurrentType,'inc');
        % lpfcDesign(rowsToremove,:) = [];
        % 
        % rowsToremove_power = ~(strcmpi(shortDesign.PreviousType,'neu') | strcmpi(shortDesign.PreviousType,'con')) | strcmpi(shortDesign.CurrentType,'inc');
        % theta_power(rowsToremove_power,:,:) = [];
        % shortDesign(rowsToremove_power,:) = [];

        shortDesign.Subject = repmat(subj,height(shortDesign),1);

        GroupDesignMatrixL = cat(1,GroupDesignMatrixL,shortDesign);

        if numel(size(theta_power)) == 2
            
            GroupPowerMatrixL{cnt} = theta_power;
        else
           
            GroupPowerMatrixL{cnt} = squeeze(mean(theta_power,2));
        end
        cnt = cnt + 1;
    end
end

GroupPowerMatrixL = GroupPowerMatrixL(1:cnt-1);
%
GroupDesignMatrixL.RT = round(GroupDesignMatrixL.RT*1000);
GroupDesignMatrixL.PC = GroupDesignMatrixL.PC - 0.33;
% Create CurrentConflict column
GroupDesignMatrixL.CurrentConflict = repmat({'NoConflict'}, height(GroupDesignMatrixL), 1); % Default to 'NoConflict'
GroupDesignMatrixL.CurrentConflict(strcmp(GroupDesignMatrixL.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc'

% Create PreviousConflict column
GroupDesignMatrixL.PreviousConflict = repmat({'NoConflict'}, height(GroupDesignMatrixL), 1); % Default to 'NoConflict'
GroupDesignMatrixL.PreviousConflict(strcmp(GroupDesignMatrixL.PreviousType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where PreviousType is 'inc'

%
GroupDesignMatrixL.CurrentConflict = categorical(GroupDesignMatrixL.CurrentConflict);
GroupDesignMatrixL.CurrentConflict = reordercats(GroupDesignMatrixL.CurrentConflict, {'NoConflict', 'Conflict'});

GroupDesignMatrixL.PreviousConflict = categorical(GroupDesignMatrixL.PreviousConflict);
GroupDesignMatrixL.PreviousConflict = reordercats(GroupDesignMatrixL.PreviousConflict, {'NoConflict', 'Conflict'});

GroupDesignMatrixL.BlockType = categorical(GroupDesignMatrixL.BlockType);
GroupDesignMatrixL.BlockType = reordercats(GroupDesignMatrixL.BlockType, {'minc', 'same', 'mcon'});
%

RTs = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);
meanRT = mean(RTs.RT)/1000;

%% TFR Plot

% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

% Initialize a cell arrays to store power data for each subject
GroupPowerMatrixL = cell(numel(files), 1);
GroupDesignMatrixL = [];
cnt = 1;
for subj = 1:numel(files)

    load(files{subj}, 'lpfcDesign','tfr','LPFC','frex','newt')
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    dlpfcChans = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));

    if isempty(dlpfcChans)
        continue
    else
        power = squeeze(mean(tfr.powspctrm(:,ismember(tfr.label,dlpfcChans),:,:),2));

        lpfcDesign([lpfcDesign.Electrode] > numel(dlpfcChans),:) = [];
    
        shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);

        shortDesign.Subject = repmat(subj,height(shortDesign),1);

        GroupDesignMatrixL = cat(1,GroupDesignMatrixL,shortDesign);

        GroupPowerMatrixL{cnt} = power;
       
        cnt = cnt + 1;
    end
end

GroupPowerMatrixL = GroupPowerMatrixL(1:cnt-1);

RTs = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);
meanRT = mean(RTs.RT);

for i = 1:length(GroupPowerMatrixL)
    comp_power(:,:,i) = squeeze(mean(GroupPowerMatrixL{i}));
end

contourf(newt,frex,mean(comp_power,3),50,'EdgeColor','none')
set(gca,'YScale','log')
yticks([2 4 8 12 30 75 150])
xticks([0 0.6 1.2]) % S-locked
% xticks([-0.5 0 0.5]) % R-locked
clim([-0.6 0.6])
% xline(0,'--k','LineWidth',1.25) % R-locked
xline(meanRT,'--k','LineWidth',1.25) % S-locked

title('dlPFC: Cue Aligned Grand Average')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
%%
% Plot the color bar
cb = colorbar('Location','southoutside');

caxis([-0.6 0.6])         
set(cb,'XTick',[-0.6 0 0.6])
cb.Label.String = '\bfPower (z)\rm';
axis off

%% Median Split by RT
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

GroupPowerMatrixL_low = cell(numel(files), 1);
GroupPowerMatrixL_high = cell(numel(files), 1);
GroupDesignMatrixL_low = [];
GroupDesignMatrixL_high = [];

cnt = 1;
for subj = 1:numel(files)

    load(files{subj}, 'lpfcDesign','tfr','LPFC','frex','newt')
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])

    % Conflict encoding
    dlpfcChans = intersect(...
        elec.label(ismember(elec.ROI,'dlPFC')), ...
        w2{1,1}.label(w2{1,1}.sig_chans) ...
    );

    if isempty(dlpfcChans)
        continue
    else
        % Extract theta power (or whichever band you prefer)
        theta_power = squeeze(mean(tfr.powspctrm(:, ...
                                 ismember(tfr.label,dlpfcChans), ...
                                 frex >= 4 & frex <= 8, :),3));

        lpfcDesign([lpfcDesign.Electrode] > numel(dlpfcChans), :) = [];
    
        rowsToremove = strcmpi(lpfcDesign.PreviousType,'None');
        shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

        lpfcDesign(rowsToremove,:) = [];
        theta_power(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];


        rowsToremove = ~(strcmpi(lpfcDesign.PreviousType,'neu') | ...
                         strcmpi(lpfcDesign.PreviousType,'con')) | ...
                          strcmpi(lpfcDesign.CurrentType,'inc');
        lpfcDesign(rowsToremove,:) = [];

        rowsToremove_power = ~(strcmpi(shortDesign.PreviousType,'neu') | ...
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
        GroupDesignMatrixL_low = cat(1, GroupDesignMatrixL_low, shortDesign_low);
        GroupDesignMatrixL_high = cat(1, GroupDesignMatrixL_high, shortDesign_high);

        % Depending on how your power data are shaped, do the same logic
        if numel(size(theta_power)) == 2
            GroupPowerMatrixL_low{cnt} = theta_power_low;
            GroupPowerMatrixL_high{cnt} = theta_power_high;
        else
            GroupPowerMatrixL_low{cnt}  = squeeze(mean(theta_power_low, 2));
            GroupPowerMatrixL_high{cnt} = squeeze(mean(theta_power_high,2));
        end
   
    cnt = cnt + 1;
    end
end

GroupPowerMatrixL_high = GroupPowerMatrixL_high(1:cnt-1);
GroupPowerMatrixL_low = GroupPowerMatrixL_low(1:cnt-1);