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

total_n = 0;
total_active = 0;
cnt = 1;
for subj = 1:numel(files)

    load(files{subj})
    
    tempstr = strsplit(files{subj},'.');
    SBJ = tempstr{1};    

    atlas_id = 'Dx';
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
    load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
    load([SBJ_vars.dirs.recon,SBJ,'_elec_',proc_id,'_pat_',atlas_id,'_final_paper.mat']);
    load([SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat'])
        
    % Task encoding
    dlpfcChans = intersect(elec.label(ismember(elec.ROI,'dlPFC')),w2{1,1}.label(w2{1,1}.sig_chans));
    
    total_n = total_n + length(elec.label(ismember(elec.ROI,'dlPFC')));

    if isempty(dlpfcChans)
        continue
    else
        nChans = numel(dlpfcChans);
        total_active = total_active + nChans;

        dlpfc_idx = ismember(hfa.label,dlpfcChans);
        
    
        [lpfcDesign,~] = fn_create_LMM_design(trial_info,nChans,1,'power');
    
        rowsToremove = strcmpi(lpfcDesign.PreviousType,'None');
        shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');
    
    
        lpfcDesign(rowsToremove,:) = [];
        hfa.ztrial(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];

        % rowsToremove = ~strcmpi(lpfcDesign.PreviousType,'inc') | strcmpi(lpfcDesign.CurrentType,'inc');
        % rowsToremove_power = ~strcmpi(shortDesign.PreviousType,'inc') | strcmpi(shortDesign.CurrentType,'inc');
        % 
        % lpfcDesign(rowsToremove,:) = [];
        % hfa.ztrial(rowsToremove_power,:,:) = [];
        % shortDesign(rowsToremove_power,:) = [];

        rowsToremove = strcmpi(lpfcDesign.PreviousType,'inc') | strcmpi(lpfcDesign.PreviousType,'con');
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'inc') | strcmpi(shortDesign.PreviousType,'con');

        lpfcDesign(rowsToremove,:) = [];
        hfa.ztrial(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];

        hfa.ztrial(:,~dlpfc_idx,:) = [];
    
        shortDesign.Subject = repmat(subj,size(shortDesign,1),1);

        GroupDesignMatrixL = cat(1,GroupDesignMatrixL,shortDesign);
        GroupPowerMatrixL{cnt} = squeeze(mean(hfa.ztrial,2));
        cnt = cnt + 1;
    end
end

GroupPowerMatrixL = GroupPowerMatrixL(1:cnt-1);

GroupDesignMatrixL.RT = round(GroupDesignMatrixL.RT*1000);
GroupDesignMatrixL.PC = GroupDesignMatrixL.PC - 0.33;

% Create CurrentConflict column
GroupDesignMatrixL.CurrentConflict = repmat({'NoConflict'}, height(GroupDesignMatrixL), 1); % Default to 'NoConflict'
GroupDesignMatrixL.CurrentConflict(strcmp(GroupDesignMatrixL.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc'

% Create PreviousConflict column
GroupDesignMatrixL.PreviousConflict = repmat({'NoConflict'}, height(GroupDesignMatrixL), 1); % Default to 'NoConflict'
GroupDesignMatrixL.PreviousConflict(strcmp(GroupDesignMatrixL.PreviousType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where PreviousType is 'inc'


GroupDesignMatrixL.CurrentConflict = categorical(GroupDesignMatrixL.CurrentConflict);
GroupDesignMatrixL.CurrentConflict = reordercats(GroupDesignMatrixL.CurrentConflict, {'NoConflict', 'Conflict'});

GroupDesignMatrixL.PreviousConflict = categorical(GroupDesignMatrixL.PreviousConflict);
GroupDesignMatrixL.PreviousConflict = reordercats(GroupDesignMatrixL.PreviousConflict, {'NoConflict', 'Conflict'});

GroupDesignMatrixL.BlockType = categorical(GroupDesignMatrixL.BlockType);
GroupDesignMatrixL.BlockType = reordercats(GroupDesignMatrixL.BlockType, {'minc', 'same', 'mcon'});


RTs = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);
meanRT = mean(RTs.RT)/1000;
