restoredefaultpath;
root_dir = '/Users/anaskhan/Desktop/';
ft_dir = [root_dir 'Apps/fieldtrip/'];
proc_id = 'main_ft';
atlas_id = 'Dx';

addpath(ft_dir);
ft_defaults;

addpath([root_dir 'PRJ_Stroop/scripts']);
addpath([root_dir 'PRJ_Stroop/scripts/utils']);
%%
% Define subjects and directories (adjust as needed)
subjectList = {'IR31','IR32','IR35','IR41','IR48','IR52','IR54',...
               'IR57','IR61','IR65','IR67','IR68','IR72','IR74','IR79','IR82','IR84'};

dataDir = fullfile(root_dir, 'PRJ_Stroop', 'results', 'Power', 'HFA', 'Stim'); 

% Get subject files
files = dir(fullfile(dataDir, '*.mat'));
files = {files.name};

% Initialize a structure to hold the prepared data
preparedData.subjectDesignMatrices = struct();
preparedData.dlPFC = {};  
preparedData.dmPFC = {};  

for f = 1:numel(files)

    currFile = fullfile(dataDir, files{f});
    load(currFile);  
    
    % Extract subject ID from filename (e.g., 'IR31.mat' → 'IR31')
    tempStr = strsplit(files{f}, '.');
    SBJ = tempStr{1};
    
    if ~ismember(SBJ, subjectList)
        continue;
    end
    
    % Run subject-specific variable file
    run(fullfile(root_dir, 'PRJ_Stroop', 'scripts', 'SBJ_vars', [SBJ '_vars.m']));
    
    % Load additional subject data (adjust paths as needed)
    load(fullfile(SBJ_vars.dirs.events, [SBJ, '_trial_info_final.mat']));
    load(fullfile(SBJ_vars.dirs.recon, [SBJ, '_elec_', proc_id, '_pat_', atlas_id, '_final_paper.mat']));
    load(fullfile(SBJ_vars.dirs.proc, [SBJ, '_smANOVA_ROI_paper_pCNI_PC_HGh_S2t2_zbt_WL05_WS002.mat']));
    
    % Identify channels for each ROI (using your criteria)
    dlPFC_chans = intersect(elec.label(ismember(elec.ROI, 'dlPFC')), ...
                            w2{1,1}.label(w2{1,1}.sig_chans));
    dmPFC_chans = intersect(elec.label(ismember(elec.ROI, 'aMCC') | ismember(elec.ROI, 'SMC')), ...
                            w2{1,1}.label(w2{1,1}.sig_chans));
    
    if isempty(dlPFC_chans) || isempty(dmPFC_chans)
        continue;
    end
    
    [designMatrix, ~] = fn_create_LMM_design(trial_info, 1, 1, 'power');

    rowsToRemove = strcmpi(designMatrix.PreviousType, 'None');

    designMatrix(rowsToRemove, :) = [];
    designMatrix.Subject = repmat({SBJ}, height(designMatrix), 1);
    
    preparedData.subjectDesignMatrices.(SBJ) = designMatrix;    
    %% Extract and store dlPFC electrode power data
    dlpfc_idx = ismember(hfa.label, dlPFC_chans);
    dlPFC_power = hfa.ztrial(~rowsToRemove, dlpfc_idx, :);  % [nTrials x nChannels x nTimePoints]

    nElectrodes = sum(dlpfc_idx);
    for elecIdx = 1:nElectrodes
        entry = struct();
        entry.Subject = SBJ;
        entry.Electrode = elecIdx;  % or store the channel label if preferred
        entry.PowerData = squeeze(dlPFC_power(:, elecIdx, :));  % [nTrials x nTimePoints]
        preparedData.dlPFC{end+1} = entry;
    end
    
    %% Extract and store dmPFC electrode power data
    dmPFC_idx = ismember(hfa.label, dmPFC_chans);
    dmPFC_power = hfa.ztrial(~rowsToRemove, dmPFC_idx, :);
    nElectrodes = sum(dmPFC_idx);
    for elecIdx = 1:nElectrodes
        entry = struct();
        entry.Subject = SBJ;
        entry.Electrode = elecIdx;
        entry.PowerData = squeeze(dmPFC_power(:, elecIdx, :));
        preparedData.dmPFC{end+1} = entry;
    end
end

% Save the prepared data
save('PreparedData.mat', 'preparedData', '-v7.3');
