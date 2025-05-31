files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

% Initialize an empty array to store results
results = [];
conditionList = {};

% Initialize a counter for results indexing
resultIdx = 0;

% Loop through each subject file
for idx = 1:length(files)
    filename = files{idx};
    subjectID = erase(filename, '.mat'); % Extract subject ID from filename
    
    % Load the .mat file
    data = load(filename);
    
    % Flags to check existence
    hasPACMtL = isfield(data, 'PACMtL');
    hasPACLtM = isfield(data, 'PACLtM');
    
    % Proceed only if PACMtL or PACLtM exists
    if hasPACMtL || hasPACLtM
        % Increment result index
        resultIdx = resultIdx + 1;
        
        % Initialize storage for averages for this subject
        avg_zMI_MtL = [];
        avg_zMI_LtM = [];
        
        % Process PACMtL if it exists
        if hasPACMtL
            PACMtL = data.PACMtL;
            nConditions = length(PACMtL); % Assuming 3 conditions
            
            % Retrieve condition names if not already set
            if isempty(conditionList)
                for cond = 1:nConditions
                    conditionList{cond} = PACMtL{cond}{1}.Condition;
                end
            end
            
            % Initialize cell array to hold average zMI per condition
            avg_zMI_MtL_cond = cell(nConditions, 1);
            
            for cond = 1:nConditions
                % Get the cell array of electrode pairs for this condition
                electrodePairs = PACMtL{cond};
                nPairs = length(electrodePairs); % Number of electrode pairs
                
                % Check if there are electrode pairs
                if nPairs == 0
                    continue; % Skip if no electrode pairs
                end
                
                % Get the size of zMI to preallocate the 3D matrix
                zMI_size = size(electrodePairs{1}.zMI); % Assuming all zMI are the same size
                zMI_matrix = zeros([zMI_size, nPairs]);
                
                for p = 1:nPairs
                    % Extract the zMI matrix for each electrode pair
                    zMI = electrodePairs{p}.zMI; % Size [low_freq_phase x high_freq_amp]
                    
                    % Concatenate zMI matrices along the third dimension (electrode pairs)
                    zMI_matrix(:, :, p) = zMI;
                end
                
                % Average over electrode pairs for this condition
                avg_zMI = mean(zMI_matrix, 3); % Mean over the third dimension (electrode pairs)
                avg_zMI_MtL_cond{cond} = avg_zMI;
            end
        else
            avg_zMI_MtL_cond = [];
        end
        
        % Process PACLtM if it exists
        if hasPACLtM
            PACLtM = data.PACLtM;
            nConditions = length(PACLtM); % Assuming 3 conditions
            
            % Initialize cell array to hold average zMI per condition
            avg_zMI_LtM_cond = cell(nConditions, 1);
            
            for cond = 1:nConditions
                % Get the cell array of electrode pairs for this condition
                electrodePairs = PACLtM{cond};
                nPairs = length(electrodePairs); % Number of electrode pairs
                
                % Check if there are electrode pairs
                if nPairs == 0
                    continue; % Skip if no electrode pairs
                end
                
                % Get the size of zMI to preallocate the 3D matrix
                zMI_size = size(electrodePairs{1}.zMI); % Assuming all zMI are the same size
                zMI_matrix = zeros([zMI_size, nPairs]);
                
                for p = 1:nPairs
                    % Extract the zMI matrix for each electrode pair
                    zMI = electrodePairs{p}.zMI; % Size [low_freq_phase x high_freq_amp]
                    
                    % Concatenate zMI matrices along the third dimension (electrode pairs)
                    zMI_matrix(:, :, p) = zMI;
                end
                
                % Average over electrode pairs for this condition
                avg_zMI = mean(zMI_matrix, 3); % Mean over the third dimension (electrode pairs)
                avg_zMI_LtM_cond{cond} = avg_zMI;
            end
        else
            avg_zMI_LtM_cond = [];
        end
        
        % Store results for this subject
        results(resultIdx).subjectID = subjectID;
        results(resultIdx).average_zMI_MtL = avg_zMI_MtL_cond;
        results(resultIdx).average_zMI_LtM = avg_zMI_LtM_cond;
    end
end


conditionListLower = lower(conditionList);

% Find indices for 'mcon' and 'minc'
mconIdx = find(strcmp(conditionListLower, 'mcon'));
mincIdx = find(strcmp(conditionListLower, 'minc'));

% Check that the conditions were found
if isempty(mconIdx) || isempty(mincIdx)
    error('Condition "mcon" or "minc" not found in conditionList.');
end

% Initialize arrays to hold aggregated data
numSubjects = length(results);

[nLowFreq, nHighFreq] = size(results(1).average_zMI_MtL{mconIdx});

aggregate_mcon_MtL = zeros(nLowFreq, nHighFreq, numSubjects);
aggregate_minc_MtL = zeros(nLowFreq, nHighFreq, numSubjects);

aggregate_mcon_LtM = zeros(nLowFreq, nHighFreq, numSubjects);
aggregate_minc_LtM = zeros(nLowFreq, nHighFreq, numSubjects);

% Loop through each subject in 'results'
for idx = 1:numSubjects
    % Extract data for MtL
    avg_zMI_MtL = results(idx).average_zMI_MtL;
    
    % Extract data for LtM
    avg_zMI_LtM = results(idx).average_zMI_LtM;
    
    aggregate_mcon_MtL(:, :, idx) = avg_zMI_MtL{mconIdx};
    aggregate_minc_MtL(:, :, idx) = avg_zMI_MtL{mincIdx};
    
    % For LtM
    aggregate_mcon_LtM(:, :, idx) = avg_zMI_LtM{mconIdx};
    aggregate_minc_LtM(:, :, idx) = avg_zMI_LtM{mincIdx};
end


% Compute the grand average across subjects for each condition and direction
% For MtL
grand_avg_mcon_MtL = mean(aggregate_mcon_MtL, 3);
grand_avg_minc_MtL = mean(aggregate_minc_MtL, 3);

% For LtM
grand_avg_mcon_LtM = mean(aggregate_mcon_LtM, 3);
grand_avg_minc_LtM = mean(aggregate_minc_LtM, 3);

% Save the aggregated data
save('aggregated_zMI_conditions.mat', ...
    'aggregate_mcon_MtL', 'aggregate_minc_MtL', ...
    'aggregate_mcon_LtM', 'aggregate_minc_LtM', ...
    'grand_avg_mcon_MtL', 'grand_avg_minc_MtL', ...
    'grand_avg_mcon_LtM', 'grand_avg_minc_LtM', ...
    'conditionList');

n_bins = 18;
LF_steps = 2:2:12;
HF_steps = 30:5:150;
contourf(LF_steps,HF_steps,grand_avg_mcon_MtL',50,'EdgeColor','none')
cb = colorbar;
clim([0 0.3])
xlim([4 12])
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')
set(cb,'YTick',[0 1 2])
cb.Label.String = ['\bf' 'PAC (z)' '\rm']; 

figure(2)
contourf(LF_steps,HF_steps,grand_avg_mcon_LtM',50,'EdgeColor','none')
cb = colorbar;
clim([0 0.3])
xlim([4 12])
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')
set(cb,'YTick',[0 1 2])
cb.Label.String = ['\bf' 'PAC (z)' '\rm']; 

for iLow = 1:length(LF_steps)
    for iHigh = 1:length(HF_steps)
        
         [~,p(iLow,iHigh),~,stats] = ttest(squeeze(aggregate_mcon_MtL(iLow,iHigh,:)),squeeze(aggregate_minc_MtL(iLow,iHigh,:)));
         t_stat(iLow,iHigh) = stats.tstat;
    end
end

figure(3)
contourf(LF_steps,HF_steps,t_stat',50,'EdgeColor','none')
cb = colorbar;
clim([-2 2])
xlim([4 12])
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')
set(cb,'YTick',[0 1 2])
cb.Label.String = ['\bf' 'PAC (z)' '\rm']; 
%%
% List of subject files
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

% Initialize an empty array to store results
results = [];
conditionList = {};

% Initialize a counter for results indexing
resultIdx = 0;

% Loop through each subject file
for idx = 1:length(files)
    filename = files{idx};
    subjectID = erase(filename, '.mat'); % Extract subject ID from filename
    
    % Load the .mat file
    data = load(filename);
    
    % Flags to check existence
    hasPACMtL = isfield(data, 'PACMtL');
    hasPACLtM = isfield(data, 'PACLtM');
    
    % Proceed only if PACMtL or PACLtM exists
    if hasPACMtL || hasPACLtM
        % Increment result index
        resultIdx = resultIdx + 1;
        
        % Initialize storage for averages for this subject
        avg_MI_MtL = [];
        avg_MI_LtM = [];
        
        % Process PACMtL if it exists
        if hasPACMtL
            PACMtL = data.PACMtL;
            nConditions = length(PACMtL); % Assuming conditions are consistent
            
            % Retrieve condition names if not already set
            if isempty(conditionList)
                for cond = 1:nConditions
                    conditionList{cond} = PACMtL{cond}{1}.Condition;
                end
            end
            
            % Initialize cell array to hold average MI per condition
            avg_MI_MtL_cond = cell(nConditions, 1);
            
            for cond = 1:nConditions
                % Get the cell array of electrode pairs for this condition
                electrodePairs = PACMtL{cond};
                nPairs = length(electrodePairs); % Number of electrode pairs
                
                % Check if there are electrode pairs
                if nPairs == 0
                    continue; % Skip if no electrode pairs
                end
                
                % Get the size of MI to preallocate the 3D matrix
                MI_size = size(electrodePairs{1}.MI); % Assuming all MI are the same size
                MI_matrix = zeros([MI_size, nPairs]);
                
                for p = 1:nPairs
                    % Extract the MI matrix for each electrode pair
                    MI = electrodePairs{p}.MI; % Size [low_freq_phase x high_freq_amp]
                    
                    % Concatenate MI matrices along the third dimension (electrode pairs)
                    MI_matrix(:, :, p) = MI;
                end
                
                % Average over electrode pairs for this condition
                avg_MI = mean(MI_matrix, 3); % Mean over the third dimension (electrode pairs)
                avg_MI_MtL_cond{cond} = avg_MI;
            end
        else
            avg_MI_MtL_cond = [];
        end
        
        % Process PACLtM if it exists
        if hasPACLtM
            PACLtM = data.PACLtM;
            nConditions = length(PACLtM); % Assuming conditions are consistent
            
            % Initialize cell array to hold average MI per condition
            avg_MI_LtM_cond = cell(nConditions, 1);
            
            for cond = 1:nConditions
                % Get the cell array of electrode pairs for this condition
                electrodePairs = PACLtM{cond};
                nPairs = length(electrodePairs); % Number of electrode pairs
                
                % Check if there are electrode pairs
                if nPairs == 0
                    continue; % Skip if no electrode pairs
                end
                
                % Get the size of MI to preallocate the 3D matrix
                MI_size = size(electrodePairs{1}.MI); % Assuming all MI are the same size
                MI_matrix = zeros([MI_size, nPairs]);
                
                for p = 1:nPairs
                    % Extract the MI matrix for each electrode pair
                    MI = electrodePairs{p}.MI; % Size [low_freq_phase x high_freq_amp]
                    
                    % Concatenate MI matrices along the third dimension (electrode pairs)
                    MI_matrix(:, :, p) = MI;
                end
                
                % Average over electrode pairs for this condition
                avg_MI = mean(MI_matrix, 3); % Mean over the third dimension (electrode pairs)
                avg_MI_LtM_cond{cond} = avg_MI;
            end
        else
            avg_MI_LtM_cond = [];
        end
        
        % Store results for this subject
        results(resultIdx).subjectID = subjectID;
        results(resultIdx).average_MI_MtL = avg_MI_MtL_cond;
        results(resultIdx).average_MI_LtM = avg_MI_LtM_cond;
    end
end

% Now, process the aggregated data for the 'mcon' and 'minc' conditions

% Convert conditionList to lower case for case-insensitive comparison
conditionListLower = lower(conditionList);

% Find indices for 'mcon' and 'minc'
mconIdx = find(strcmp(conditionListLower, 'mcon'));
mincIdx = find(strcmp(conditionListLower, 'minc'));

% Check that the conditions were found
if isempty(mconIdx) || isempty(mincIdx)
    error('Condition "mcon" or "minc" not found in conditionList.');
end

% Initialize arrays to hold aggregated data
numSubjects = length(results);

% Assuming all MI matrices have the same dimensions, get size from first subject
[nLowFreq, nHighFreq] = size(results(1).average_MI_MtL{mconIdx});

% Initialize aggregation arrays
aggregate_mcon_MtL = zeros(nLowFreq, nHighFreq, numSubjects);
aggregate_minc_MtL = zeros(nLowFreq, nHighFreq, numSubjects);

aggregate_mcon_LtM = zeros(nLowFreq, nHighFreq, numSubjects);
aggregate_minc_LtM = zeros(nLowFreq, nHighFreq, numSubjects);

% Loop through each subject in 'results'
for idx = 1:numSubjects
    % Extract data for MtL
    avg_MI_MtL = results(idx).average_MI_MtL;
    
    % Extract data for LtM
    avg_MI_LtM = results(idx).average_MI_LtM;
    
    % For MtL
    aggregate_mcon_MtL(:, :, idx) = avg_MI_MtL{mconIdx};
    aggregate_minc_MtL(:, :, idx) = avg_MI_MtL{mincIdx};
    
    % For LtM
    aggregate_mcon_LtM(:, :, idx) = avg_MI_LtM{mconIdx};
    aggregate_minc_LtM(:, :, idx) = avg_MI_LtM{mincIdx};
end

% Compute the grand average across subjects for each condition and direction
% For MtL
grand_avg_mcon_MtL = mean(aggregate_mcon_MtL, 3);
grand_avg_minc_MtL = mean(aggregate_minc_MtL, 3);

% For LtM
grand_avg_mcon_LtM = mean(aggregate_mcon_LtM, 3);
grand_avg_minc_LtM = mean(aggregate_minc_LtM, 3);


n_bins = 18;
LF_steps = 2:2:12;
HF_steps = 30:5:150;

contourf(LF_steps,HF_steps,grand_avg_mcon_MtL',50,'EdgeColor','none')
cb = colorbar;
clim([0 0.0006])
xlim([4 12])
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')
% set(cb,'YTick',[0 0.0005 0.001])
cb.Label.String = ['\bf' 'Modulation Index' '\rm']; 
title('Low Conflict Load: dACC -> dlPFC PAC')

figure(2)
contourf(LF_steps,HF_steps,grand_avg_mcon_LtM',50,'EdgeColor','none')
cb = colorbar;
clim([0 0.0006])
xlim([4 12])
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')
% set(cb,'YTick',[0 0.0005 0.001])
cb.Label.String = ['\bf' 'Modulation Index' '\rm']; 
title('Low Conflict Load: dlPFC -> dACC PAC')

% Stats

for iLow = 1:length(LF_steps)
    for iHigh = 1:length(HF_steps)
        dCFC(iLow,iHigh,:) = squeeze((aggregate_minc_LtM(iLow,iHigh,:) - aggregate_minc_MtL(iLow,iHigh,:))...
            ./(aggregate_minc_LtM(iLow,iHigh,:) + aggregate_minc_MtL(iLow,iHigh,:)));
        % isNorm(iLow,iHigh) = swtest(dCFC(iLow,iHigh,:));
        [~,p(iLow,iHigh)] = ttest(squeeze(dCFC(iLow,iHigh,:)));
    end
end