%% dlPFC Theta/Beta Conflict
unique_sbjs = unique(GroupDesignMatrixL.Subject);
for subj = 1:numel(GroupPowerMatrixL)
    temptrials = GroupDesignMatrixL([GroupDesignMatrixL.Subject] == unique_sbjs(subj),:);
    
    tempcon = squeeze(mean(GroupPowerMatrixL{subj}([temptrials.CurrentConflict] == 'NoConflict',:,:)));
    tempinc = squeeze(mean(GroupPowerMatrixL{subj}([temptrials.CurrentConflict] == 'Conflict',:,:)));

   
    con(:,subj) = tempcon;
    inc(:,subj) = tempinc;

end

[hl, ~] = boundedline(newt,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',newt,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Beta Power (z)') % Beta
% ylabel('Theta Power (z)') % Theta
xlabel('Time (s)')

% xline(meanRT,'--k','LineWidth',1.25)
% xlim([-0.103 1.253])
% xticks([0 0.6 1.2])

xline(0,'--k','LineWidth',1.25)
xlim([-0.753 0.753])
xticks([-0.5 0 0.5])

% yticks([-0.4 0 0.4 0.8])
% yticks([-0.8 -0.4 0 0.4]) % Beta S
yticks([-0.3 0 0.3]) % Beta R


title('dlPFC: Response Aligned Beta Power')
set(gca,'FontSize',18)

% Include stat line

cluster_thresh = prctile(permConflict, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empConflictTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

% p-values
max_p = max((sum(abs(empConflictTFCE(sigIdx)) < permConflict,2)) ./ 1000)

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1);         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos); sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -0.7; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 5);
end
set(gca,'LineWidth',1)
set(gca,'Box','off')
set(gca,'FontSize',18)
ylim([-0.72 0.6])
%% dlPFC HFA Conflict Block
unique_sbjs = unique(GroupDesignMatrixL.Subject);
blocktypes = {'minc', 'same', 'mcon'};

for subj = 1:numel(GroupPowerMatrixL)
    % Select trials for the current subject
    temptrials = GroupDesignMatrixL(GroupDesignMatrixL.Subject == unique_sbjs(subj), :);

    for b = 1:length(blocktypes)
        blocktype = blocktypes{b};

        % Get indices of trials for the current block type
        block_trials = temptrials.BlockType == blocktype;

        % Indices for 'NoConflict' and 'Conflict' within the current block type
        con_trials = temptrials.CurrentConflict == 'NoConflict' & block_trials;
        inc_trials = temptrials.CurrentConflict == 'Conflict' & block_trials;

        % Compute the mean over trials for 'NoConflict' and 'Conflict'
        tempcon = squeeze(mean(GroupPowerMatrixL{subj}(con_trials, :), 1));
        tempinc = squeeze(mean(GroupPowerMatrixL{subj}(inc_trials, :), 1));


        % Store the results in a 3D array (time x subject x block type)
        con(:, subj, b) = tempcon;
        inc(:, subj, b) = tempinc;

        clear smoothcon smoothinc tempcon tempinc
    end
end
con(:,end,:) = [];
inc(:,end,:) = [];

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
newt = hfa.newtime;
for b = 1:length(blocktypes)
    figure(b)
    [hl, ~] = boundedline(newt,mean(con(:,:,b),2),std(con(:,:,b),[],2)./sqrt(size(con(:,:,b),2)),'b',newt,mean(inc(:,:,b),2),std(inc(:,:,b),[],2)./sqrt(size(inc(:,:,b),2)),'r');
    hl(1).LineWidth = 1.25;
    hl(2).LineWidth = 1.25;
    ylabel('HFA (z)')
    xlabel('Time (s)')

    % xline(meanRT,'--k','LineWidth',1.25)
    % xlim([-0.103 1.253])
    % xticks([0 0.6 1.2])
    % ylim([-0.5 4.5])
    % yticks([0 2.25 4.5])


    xline(0,'--k','LineWidth',1.25)
    xlim([-0.753 0.753])
    xticks([-0.5 0 0.5])
    ylim([-1 4.5])
    yticks([0 2.25 4.5])
   
    if strcmpi(blocktypes{b},'mcon')
        title('dlPFC: Response Aligned HFA 50%')
    elseif strcmpi(blocktypes{b},'same')
        title('dlPFC: Response Aligned HFA 33%')
    else
        title('dlPFC: Response Aligned HFA 16%')
    end
    set(gca,'FontSize',24)
    set(gca,'LineWidth',1)
    set(gca,'Box','off')
end
% Include stat line

cluster_thresh = prctile(permConflictBlock, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empCBTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

% p-values
max_p = max((sum(abs(empCBTFCE(sigIdx)) < permConflictBlock,2)) ./ 1000)

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1),         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos), sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -0.5; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 5);
end
ylim([-0.7 4.5])

%% dlPFC Theta/Beta Conflict Block
unique_sbjs = unique(GroupDesignMatrixL.Subject);
blocktypes = {'minc', 'same', 'mcon'};

for subj = 1:numel(GroupPowerMatrixL)
    % Select trials for the current subject
    temptrials = GroupDesignMatrixL(GroupDesignMatrixL.Subject == unique_sbjs(subj), :);

    for b = 1:length(blocktypes)
        blocktype = blocktypes{b};

        % Get indices of trials for the current block type
        block_trials = temptrials.BlockType == blocktype;

        % Indices for 'NoConflict' and 'Conflict' within the current block type
        con_trials = temptrials.CurrentConflict == 'NoConflict' & block_trials;
        % con_trials = strcmpi(temptrials.CurrentType,'con') & block_trials;

        inc_trials = temptrials.CurrentConflict == 'Conflict' & block_trials;

        % Compute the mean over trials for 'NoConflict' and 'Conflict'
        tempcon = squeeze(mean(GroupPowerMatrixL{subj}(con_trials, :), 1));
        tempinc = squeeze(mean(GroupPowerMatrixL{subj}(inc_trials, :), 1));


        % Store the results in a 3D array (time x subject x block type)
        con(:, subj, b) = tempcon;
        inc(:, subj, b) = tempinc;

        clear smoothcon smoothinc tempcon tempinc
    end
end
con(:,end,:) = [];
inc(:,end,:) = [];

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

for b = 1:length(blocktypes)
    figure(b)
    [hl, ~] = boundedline(newt,mean(con(:,:,b),2),std(con(:,:,b),[],2)./sqrt(size(con(:,:,b),2)),'b',newt,mean(inc(:,:,b),2),std(inc(:,:,b),[],2)./sqrt(size(inc(:,:,b),2)),'r');
    hl(1).LineWidth = 1.25;
    hl(2).LineWidth = 1.25;
    ylabel('Beta Power (z)')
    % ylabel('Theta Power (z)')

    xlabel('Time (s)')
    
    xlim([-0.103 1.253])
    xticks([0 0.6 1.2]);
    xline(meanRT,'--k','LineWidth',1.25)

    % xlim([-0.7503 0.75])
    % xticks([-0.5 0 0.5]);
    % xline(0,'--k','LineWidth',1.25)

    % Beta S lims
    ylim([-1.25 0.75])
    yticks([-1 -0.5 0 0.5])

    % Theta S Lims
    % ylim([-1.4 1])
    % yticks([-1 0 1])

    % Theta R lims
    % ylim([-1.2 0.75])
    % yticks([-0.75 0 0.75])

    % Beta R lims
    % ylim([-0.9 0.4])
    % yticks([-0.5 -0.25 0 0.25])

    if strcmpi(blocktypes{b},'mcon')
        title('dlPFC: Stimulus Aligned Beta 50%')
    elseif strcmpi(blocktypes{b},'same')
        title('dlPFC: Stimulus Aligned Beta 33%')
    else
        title('dlPFC: Stimulus Aligned Beta 16%')
    end
    set(gca,'FontSize',24)
    set(gca,'LineWidth',1)
    set(gca,'Box','off')
end

% Include stat line

cluster_thresh = prctile(permConflictBlock, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empCBTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

% p-values
% max_p = max((sum(abs(empCBTFCE(sigIdx)) < permConflictBlock,2)) ./ 1000)

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1);         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos); sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
% yLevel = -0.75; %Theta S
% yLevel = -1; %Theta R

% yLevel = -0.8; %Beta R
yLevel = -1.25; %Beta S
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 5);
end

%% dmPFC HFA Conflict Block
unique_sbjs = unique(GroupDesignMatrixM.Subject);
blocktypes = {'minc', 'same', 'mcon'};

for subj = 1:numel(GroupPowerMatrixM)
    % Select trials for the current subject
    temptrials = GroupDesignMatrixM(GroupDesignMatrixM.Subject == unique_sbjs(subj), :);

    for b = 1:length(blocktypes)
        blocktype = blocktypes{b};

        % Get indices of trials for the current block type
        block_trials = temptrials.BlockType == blocktype;

        % Indices for 'NoConflict' and 'Conflict' within the current block type
        con_trials = temptrials.CurrentConflict == 'NoConflict' & block_trials;
        inc_trials = temptrials.CurrentConflict == 'Conflict' & block_trials;

        % Compute the mean over trials for 'NoConflict' and 'Conflict'
        tempcon = squeeze(mean(GroupPowerMatrixM{subj}(con_trials, :), 1));
        tempinc = squeeze(mean(GroupPowerMatrixM{subj}(inc_trials, :), 1));


        % Store the results in a 3D array (time x subject x block type)
        con(:, subj, b) = tempcon;
        inc(:, subj, b) = tempinc;

        clear smoothcon smoothinc tempcon tempinc
    end
end
con(:,end,:) = [];
inc(:,end,:) = [];

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
newt = hfa.newtime;
for b = 1:length(blocktypes)
    figure(b)
    [hl, ~] = boundedline(newt,mean(con(:,:,b),2),std(con(:,:,b),[],2)./sqrt(size(con(:,:,b),2)),'b',newt,mean(inc(:,:,b),2),std(inc(:,:,b),[],2)./sqrt(size(inc(:,:,b),2)),'r');
    hl(1).LineWidth = 1.25;
    hl(2).LineWidth = 1.25;
    ylabel('HFA (z)')
    xlabel('Time (s)')
    
    
    % xlim([-0.103 1.253])
    % xticks([0 0.6 1.2]);
    % ylim([-0.7 3.5])
    % yticks([0 1.75 3.5])
    % xline(meanRT,'--k','LineWidth',1.25)

    xlim([-0.753 0.753])
    xticks([-0.5 0 0.5])
    xline(0,'--k','LineWidth',1.25)

    ylim([-1 3.5])
    yticks([0 1.75 3.5])
    
    if strcmpi(blocktypes{b},'mcon')
        title('dmPFC: Response Aligned HFA 50%')
    elseif strcmpi(blocktypes{b},'same')
        title('dmPFC: Response Aligned HFA 33%')
    else
        title('dmPFC: Response Aligned HFA 16%')
    end
    set(gca,'FontSize',24)
    set(gca,'LineWidth',1)
    set(gca,'Box','off')
end

% Include stat line

cluster_thresh = prctile(permConflictBlock, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empCBTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1),         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos), sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -0.4; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 5);
end

ylim([-1.2 3.5])
%% dmPFC Theta/Beta Conflict Block
unique_sbjs = unique(GroupDesignMatrixM.Subject);
blocktypes = {'minc', 'same', 'mcon'};

for subj = 1:numel(GroupPowerMatrixM)
    % Select trials for the current subject
    temptrials = GroupDesignMatrixM(GroupDesignMatrixM.Subject == unique_sbjs(subj), :);

    for b = 1:length(blocktypes)
        blocktype = blocktypes{b};

        % Get indices of trials for the current block type
        block_trials = temptrials.BlockType == blocktype;

        % Indices for 'NoConflict' and 'Conflict' within the current block type
        con_trials = temptrials.CurrentConflict == 'NoConflict' & block_trials;
        inc_trials = temptrials.CurrentConflict == 'Conflict' & block_trials;

        % Compute the mean over trials for 'NoConflict' and 'Conflict'
        tempcon = squeeze(mean(GroupPowerMatrixM{subj}(con_trials, :), 1));
        tempinc = squeeze(mean(GroupPowerMatrixM{subj}(inc_trials, :), 1));


        % Store the results in a 3D array (time x subject x block type)
        con(:, subj, b) = tempcon;
        inc(:, subj, b) = tempinc;

        clear smoothcon smoothinc tempcon tempinc
    end
end
con(:,end,:) = [];
inc(:,end,:) = [];

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

for b = 1:length(blocktypes)
    figure(b)
    [hl, ~] = boundedline(newt,mean(con(:,:,b),2),std(con(:,:,b),[],2)./sqrt(size(con(:,:,b),2)),'b',newt,mean(inc(:,:,b),2),std(inc(:,:,b),[],2)./sqrt(size(inc(:,:,b),2)),'r');
    hl(1).LineWidth = 1.25;
    hl(2).LineWidth = 1.25;
    ylabel('Beta Power (z)')
    % ylabel('Theta Power (z)')
    xlabel('Time (s)')

    % xlim([-0.103 1.253])
    % xticks([0 0.6 1.2]);
    % xline(meanRT,'--k','LineWidth',1.25)

    xlim([-0.753 0.753])
    xticks([-0.5 0 0.5]);
    xline(0,'--k','LineWidth',1.25)

    % Theta R
    % ylim([-0.5 1.25])
    % yticks([0 0.625 1.25])

    % Theta S
    % ylim([-1.25 1.5])
    % yticks([-1 0 1])

    % Beta R
    ylim([-0.6 0.5])
    yticks([-0.4 0 0.4])

    % Beta S
    % ylim([-1.2 1])
    % yticks([-0.8 0 0.8])


    if strcmpi(blocktypes{b},'mcon')
        title('dmPFC: Cue Aligned Beta 50%')
    elseif strcmpi(blocktypes{b},'same')
        title('dmPFC: Cue Aligned Beta 33%')
    else
        title('dmPFC: Cue Aligned Beta 16%')
    end
    set(gca,'FontSize',24)
    set(gca,'LineWidth',1)
    set(gca,'Box','off')
end

% Include stat line

cluster_thresh = prctile(permConflictBlock, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empCBTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1),         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos), sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -1; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 5);
end
%% dlPFC Block type

unique_sbjs = unique(GroupDesignMatrixL.Subject);
for subj = 1:numel(GroupPowerMatrixL)
    temptrials = GroupDesignMatrixL([GroupDesignMatrixL.Subject] == unique_sbjs(subj),:);
    
    tempmc = squeeze(mean(GroupPowerMatrixL{subj}([temptrials.BlockType] == 'mcon',:,:)));
    tempmi = squeeze(mean(GroupPowerMatrixL{subj}([temptrials.BlockType] == 'minc',:,:)));
    
    mc(:,subj) = tempmc;
    mi(:,subj) = tempmi;

end
clear tempinc tempcon smoothinc smoothcon

mi(:,end) = [];
mc(:,end) = [];

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

[hl, ~] = boundedline(newt,mean(mc,2),std(mc,[],2)./sqrt(size(mc,2)),'b',newt,mean(mi,2),std(mi,[],2)./sqrt(size(mi,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Theta Power (z)')
xlabel('Time (s)')
xline(meanRT,'--k','LineWidth',1.25)
xlim([-0.203 1.25])
xlim([statTime(1)-0.003 statTime(end)+0.003])
ylim([-0.55 0.65])
yticks([-0.5 0 0.5])
xticks([0 0.6 1.2])
title('dlPFC: Cue Aligned Theta Power')
set(gca,'FontSize',16)

%% dmPFC Theta/Beta Conflict

unique_sbjs = unique(GroupDesignMatrixM.Subject);
for subj = 1:numel(GroupPowerMatrixM)
    temptrials = GroupDesignMatrixM([GroupDesignMatrixM.Subject] == unique_sbjs(subj),:);

    tempcon = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.CurrentConflict] == 'NoConflict',:,:)));
    tempinc = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.CurrentConflict] == 'Conflict',:,:)));

    con(:,subj) = tempcon;
    inc(:,subj) = tempinc;

end

[hl, ~] = boundedline(newt,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',newt,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Beta Power (z)') % Beta
% ylabel('Theta Power (z)') % Theta
xlabel('Time (s)')

% xlim([-0.103 1.253])
% xticks([0 0.6 1.2])
% xline(meanRT,'--k','LineWidth',1.25)

xlim([-0.753 0.753])
xticks([-0.5 0 0.5])
xline(0,'--k','LineWidth',1.25)

% ylim([-0.7 0.4]); yticks([-0.6 -0.3 0 0.3]) % Beta S
ylim([-0.4 0.2]); yticks([-0.4 -0.2 0 0.2]) % Beta R

% ylim([-0.4 1]); yticks([0 0.6 1.2]) %Theta
title('dmPFC: Response Aligned Beta Power')
% title('dmPFC: Response Aligned Theta Power')
set(gca,'FontSize',18)

% Include stat line

cluster_thresh = prctile(permConflict, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empConflictTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1),         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos), sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -0.35; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 2.5);
end
ylim([-0.6 1.2])

set(gca,'LineWidth',1)
set(gca,'Box','off')
set(gca,'FontSize',18)
%% dmPFC Block type

unique_sbjs = unique(GroupDesignMatrixM.Subject);
for subj = 1:numel(GroupPowerMatrixM)
    temptrials = GroupDesignMatrixM([GroupDesignMatrixM.Subject] == unique_sbjs(subj),:);
    
    tempmc = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.BlockType] == 'mcon' & [temptrials.PreviousType] == 'neu',:,:)));
    tempmi = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.BlockType] == 'minc' & [temptrials.PreviousType] == 'neu',:,:)));
    
    mc(:,subj) = tempmc;
    mi(:,subj) = tempmi;

end
clear tempinc tempcon smoothinc smoothcon

mi(:,end) = [];
mc(:,end) = [];

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

[hl, ~] = boundedline(newt,mean(mc,2),std(mc,[],2)./sqrt(size(mc,2)),newt,mean(mi,2),std(mi,[],2)./sqrt(size(mi,2)),'cmap',[236 0 140; 20 170 75]./255,'alpha');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Theta Power (z)')
xlabel('Time (s)')
xline(meanRT,'--k','LineWidth',1.25)
xlim([-0.203 1.25])
xlim([statTime(1)-0.003 statTime(end)+0.003])
ylim([-0.55 0.65])
yticks([-0.5 0 0.5])
xticks([0 0.6 1.2])
title('dlPFC: Cue Aligned Theta Power')
% Include stat line

cluster_thresh = prctile(permBlock, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empBlockTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

% p-values
max_p = max((sum(abs(empBlockTFCE(sigIdx)) < permBlock,2)) ./ 1000)


sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1);         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos); sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -2.25; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 5);
end
ylim([-0.6 1.2])


%% dmPFC Previous type

unique_sbjs = unique(GroupDesignMatrixM.Subject);
for subj = 1:numel(GroupPowerMatrixM)
    temptrials = GroupDesignMatrixM([GroupDesignMatrixM.Subject] == unique_sbjs(subj),:);
    
    tempmc = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.PreviousConflict] == 'NoConflict',:,:)));
    tempmi = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.PreviousConflict] == 'Conflict',:,:)));
    
    mc(:,subj) = tempmc;
    mi(:,subj) = tempmi;

end
clear tempinc tempcon smoothinc smoothcon

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

[hl, ~] = boundedline(newt,mean(mc,2),std(mc,[],2)./sqrt(size(mc,2)),'b',newt,mean(mi,2),std(mi,[],2)./sqrt(size(mi,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Theta Power (z)')
xlabel('Time (s)')
xline(meanRT,'--k','LineWidth',1.25)
xlim([-0.203 1.25])
xlim([statTime(1)-0.003 statTime(end)+0.003])
ylim([-0.55 0.65])
yticks([-0.5 0 0.5])
xticks([0 0.6 1.2])
title('dlPFC: Cue Aligned Theta Power')
set(gca,'FontSize',16)
%% dmPFC HFA Conflict
unique_sbjs = unique(GroupDesignMatrixM.Subject);
for subj = 1:numel(GroupPowerMatrixM)
    temptrials = GroupDesignMatrixM([GroupDesignMatrixM.Subject] == unique_sbjs(subj),:);
    
    tempcon = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.CurrentConflict] == 'NoConflict',:)));
    tempinc = squeeze(mean(GroupPowerMatrixM{subj}([temptrials.CurrentConflict] == 'Conflict',:)));

   
    con(:,subj) = tempcon;
    inc(:,subj) = tempinc;

end
clear tempinc tempcon 
newt = hfa.newtime;
addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
[hl, ~] = boundedline(newt,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',newt,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('HFA (z)')
xlabel('Time (s)')

xlim([-0.103 1.253])
xticks([0 0.6 1.2])
xline(meanRT,'--k','LineWidth',1.25)
yticks([0 1.5 3])

xlim([-0.753 0.753])
xticks([-0.5 0 0.5])
xline(0,'--k','LineWidth',1.25)
yticks([0 1.5 3])

title('dmPFC: Response Aligned HFA')

set(gca,'LineWidth',1)
set(gca,'Box','off')
set(gca,'FontSize',18)
% Include stat line
cluster_thresh = prctile(permConflict, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empConflictTFCE) > cluster_thresh;

sigIdx = find(sigTimes);

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1),         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos), sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -0.5; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 5);
end
ylim([-0.75 3])
%% dlPFC HFA Conflict

unique_sbjs = unique(GroupDesignMatrixL.Subject);
for subj = 1:numel(GroupPowerMatrixL)
    temptrials = GroupDesignMatrixL([GroupDesignMatrixL.Subject] == unique_sbjs(subj),:);
    
    tempcon = squeeze(mean(GroupPowerMatrixL{subj}([temptrials.CurrentConflict] == 'NoConflict',:)));
    tempinc = squeeze(mean(GroupPowerMatrixL{subj}([temptrials.CurrentConflict] == 'Conflict',:)));

   
    con(:,subj) = tempcon;
    inc(:,subj) = tempinc;

end
clear tempinc tempcon 

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
newt = hfa.newtime;
[hl, ~] = boundedline(newt,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',newt,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('HFA (z)')
xlabel('Time (s)')


xlim([-0.103 1.253])
xticks([0 0.6 1.2])
xline(meanRT,'--k','LineWidth',1.25)
yticks([0 1.5 3])

xlim([-0.753 0.753])
xticks([-0.5 0 0.5])
xline(0,'--k','LineWidth',1.25)
yticks([0 1.75 3.5])


title('dlPFC: Response Aligned HFA')
set(gca,'FontSize',18)

% Include stat line

cluster_thresh = prctile(permConflict, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empConflictTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

sigIdx = sort(sigIdx);  
d      = diff(sigIdx);

boundaryPos = find(d > 1);  % Indices where a gap > 1 occurs

if isempty(boundaryPos)
    % === Only one cluster ===
    clusterStarts = sigIdx(1);
    clusterEnds   = sigIdx(end);
else
    % === Multiple clusters ===
    clusterStarts = [sigIdx(1),         sigIdx(boundaryPos+1)];
    clusterEnds   = [sigIdx(boundaryPos), sigIdx(end)];
end

% Plot horizontal lines for each cluster
hold on;
yLevel = -0.5; 
for i = 1:numel(clusterStarts)
    x1 = newt(clusterStarts(i));
    x2 = newt(clusterEnds(i));
    plot([x1 x2], [yLevel yLevel], 'k', 'LineWidth', 2.5);
end

set(gca,'LineWidth',1)
set(gca,'Box','off')
set(gca,'FontSize',18)

ylim([-0.7 3.5])
%% dmPFC Median Split RT Gratton

for subj = 1:numel(GroupPowerMatrixM_high)
    
    temphigh = squeeze(mean(GroupPowerMatrixM_high{subj}));
    templow = squeeze(mean(GroupPowerMatrixM_low{subj}));

    high(:,subj) = temphigh;
    low(:,subj) = templow;

end

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
[hl, ~] = boundedline(newt,mean(low,2),std(low,[],2)./sqrt(size(low,2)),newt,mean(high,2),std(high,[],2)./sqrt(size(high,2)),'cmap',[0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560]);
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Theta Power (z)')
xlabel('Time (s)')
title({'dmPFC: Cue Aligned', 'Median-Split Theta C-NC'})

xlim([-0.103 1.253])
xticks([0 0.6 1.2])

yticks([-0.4 0 0.4 0.8])

set(gca,'LineWidth',1)
set(gca,'Box','off')
set(gca,'FontSize',18)

%% RT Stats

cluster_thresh = prctile(permRT, 100 - (100*0.025));

% Identify significant times
sigTimes = abs(empRTTFCE) > cluster_thresh;
sigIdx = find(sigTimes);

% p-values
max_p = max((sum(abs(empRTTFCE(sigIdx)) < permRT,2)) ./ 1000)

figure;
plot(newt, empRTT, 'r'); 
hold on;

% Find differences between successive indices
diffSig = diff(sigIdx);

% Find where a break (i.e., non-contiguous jump) occurs. A break is any gap > 1.
breakIndices = find(diffSig > 1);

% Determine the starting and ending indices of each contiguous segment.
% The first segment always starts at the first element of sigIdx.
segStart = [sigIdx(1); sigIdx(breakIndices + 1)];
% The last segment always ends at the last element of sigIdx.
segEnd = [sigIdx(breakIndices); sigIdx(end)];

% Now, plot each segment with a thicker line.
for i = 1:length(segStart)

    curPoints = segStart(i):segEnd(i);
    
    % Plot this segment in red with a thicker line (e.g., LineWidth of 2 or 3)
    plot(newt(curPoints), empRTT(curPoints), 'r', 'LineWidth', 3);
end

hold off;
yline(0,'--k','LineWidth',1.25)
title('dmPFC: Stimulus Aligned Theta ~ RT')
ylabel('Wald t')
xlabel('Time (s)')
yticks([-6 -3 0 3 6])
xlim([-0.103 1.253])
xticks([0 0.6 1.2])
set(gca,'FontSize',22)
%% dlPFC Median Split RT Gratton

for subj = 1:numel(GroupPowerMatrixL_high)
    
    temphigh = squeeze(mean(GroupPowerMatrixL_high{subj}));
    templow = squeeze(mean(GroupPowerMatrixL_low{subj}));

    high(:,subj) = temphigh;
    low(:,subj) = templow;

end

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
[hl, ~] = boundedline(newt,mean(low,2),std(low,[],2)./sqrt(size(low,2)),newt,mean(high,2),std(high,[],2)./sqrt(size(high,2)),'cmap',[0.3010 0.7450 0.9330; 0.4940 0.1840 0.5560]);
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Theta Power (z)')
xlabel('Time (s)')
title({'dlPFC: Cue Aligned', 'Median-Split Theta NC-NC'})

xlim([-0.103 1.253])
xticks([0 0.6 1.2])

yticks([-0.4 0 0.4 0.8])

set(gca,'LineWidth',1)
set(gca,'Box','off')
set(gca,'FontSize',18)

%% ITI ISPC

trialsidx = GroupDesign([GroupDesign.ElectrodePair] == 1,:);

mc = [];
mi = [];

unique_sbjs = unique(trialsidx.Subject);
for i = 1:numel(GroupISPC)
    temptrials = trialsidx([trialsidx.Subject] == unique_sbjs(i),:);
    mc(:,i) = squeeze(mean(GroupISPC{i}([temptrials.BlockType] == 'mcon' & strcmpi(temptrials.PreviousType,'neu'),:,:)));
    mi(:,i) = squeeze(mean(GroupISPC{i}([temptrials.BlockType] == 'minc' & strcmpi(temptrials.PreviousType,'neu'),:,:)));
end

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
[hl, ~] = boundedline(newt,mean(mc,2),std(mc,[],2)./sqrt(size(mc,2)),'b',newt,mean(mi,2),std(mi,[],2)./sqrt(size(mi,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Theta ISPC (z)')
xlabel('Time (s)')
title('dACC-dlPFC: Cue Aligned Theta ISPC')
xline(meanRT,'--k')
yticks([0 1.5 3.5])
xticks([0 0.625 1.25])
xlim([-0.203 1.25])
set(gca,'FontSize',16)