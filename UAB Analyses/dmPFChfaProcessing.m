%% HFA Stats

% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];
proc_id = 'main_ft';

% Initialize a cell arrays to store power data for each subject
GroupPowerMatrixM = cell(numel(files), 1);
GroupDesignMatrixM = [];

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
    daccChans = intersect(elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC')),w2{1,1}.label(w2{1,1}.sig_chans));

    total_n = total_n + length(elec.label(ismember(elec.ROI,'aMCC') | ismember(elec.ROI,'SMC')));

    if isempty(daccChans)
        continue
    else
        nChans = numel(daccChans);

        total_active = total_active + nChans;

        dacc_idx = ismember(hfa.label,daccChans);
        
    
        [mpfcDesign,~] = fn_create_LMM_design(trial_info,nChans,1,'power');
    
        rowsToremove = strcmpi(mpfcDesign.PreviousType,'None');
        shortDesign = mpfcDesign(mpfcDesign.Electrode == 1,:);
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');
    
    
        mpfcDesign(rowsToremove,:) = [];
        hfa.ztrial(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];

        % rowsToremove = ~strcmpi(mpfcDesign.PreviousType,'inc') | strcmpi(mpfcDesign.CurrentType,'inc');
        % rowsToremove_power = ~strcmpi(shortDesign.PreviousType,'inc') | strcmpi(shortDesign.CurrentType,'inc');
        % 
        % 
        % mpfcDesign(rowsToremove,:) = [];
        % hfa.ztrial(rowsToremove_power,:,:) = [];
        % shortDesign(rowsToremove_power,:) = [];

        rowsToremove = strcmpi(mpfcDesign.PreviousType,'inc') | strcmpi(mpfcDesign.PreviousType,'con');
        rowsToremove_power = strcmpi(shortDesign.PreviousType,'inc') | strcmpi(shortDesign.PreviousType,'con');


        mpfcDesign(rowsToremove,:) = [];
        hfa.ztrial(rowsToremove_power,:,:) = [];
        shortDesign(rowsToremove_power,:) = [];

        hfa.ztrial(:,~dacc_idx,:) = [];

        shortDesign.Subject = repmat(subj,size(shortDesign,1),1);

        GroupDesignMatrixM = cat(1,GroupDesignMatrixM,shortDesign);
        GroupPowerMatrixM{cnt} = squeeze(mean(hfa.ztrial,2));
        cnt = cnt + 1;
    end
end

GroupPowerMatrixM = GroupPowerMatrixM(1:cnt-1);
%
GroupDesignMatrixM.RT = round(GroupDesignMatrixM.RT*1000);
GroupDesignMatrixM.PC = GroupDesignMatrixM.PC - 0.33;
%
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