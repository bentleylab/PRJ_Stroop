% Get a list of all files and folders in the current directory
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

% Initialize a cell arrays to store power data for each subject
GroupISPC = cell(numel(files), 1);
GroupDesign = [];

for subj = 1:numel(files)

    load(files{subj})
    ztISPC = extract_theta_ISPC(zISPC, frex);
    clear zISPC
    
    ztISPC = permute(cat(3,ztISPC{:}),[1 3 2]);
    ztISPC = squeeze(mean(ztISPC,2));

    rowsToremove = strcmpi(coherenceDesign.PreviousType,'None');
    shortDesign = coherenceDesign(coherenceDesign.ElectrodePair == 1,:);
    rowsToremove_ispc = strcmpi(shortDesign.PreviousType,'None');

    coherenceDesign(rowsToremove,:) = [];
    ztISPC(rowsToremove_ispc,:,:) = [];
    shortDesign(rowsToremove_ispc,:) = [];

    % rowsToremove = strcmpi(coherenceDesign.CurrentType,'neu') | strcmpi(coherenceDesign.PreviousType,'neu');% | strcmpi(coherenceDesign.CurrentType,'inc');
    % coherenceDesign(rowsToremove,:) = [];
    % 
    % rowsToremove_ispc = strcmpi(shortDesign.CurrentType,'neu') | strcmpi(shortDesign.PreviousType,'neu');% | strcmpi(shortDesign.CurrentType,'inc');
    % ztISPC(rowsToremove_ispc,:,:) = [];
    % shortDesign(rowsToremove_ispc,:) = [];

    % rowsToremove = ~strcmpi(coherenceDesign.BlockType,'same');
    % coherenceDesign(rowsToremove,:) = [];
    % 
    % rowsToremove_ispc = ~strcmpi(shortDesign.BlockType,'same');
    % ztISPC(rowsToremove_ispc,:,:) = [];
    % shortDesign(rowsToremove_ispc,:) = [];

    % Focus on only pI trials
    % rowsToremove = strcmpi(coherenceDesign.PreviousType,'con');
    % coherenceDesign(rowsToremove,:) = [];
    % 
    % rowsToremove_ispc = strcmpi(shortDesign.PreviousType,'con');
    % ztISPC(rowsToremove_ispc,:,:) = [];
    % shortDesign(rowsToremove_ispc,:) = [];

    shortDesign.Subject = repmat(subj,height(shortDesign),1);

    GroupDesign = cat(1,GroupDesign,shortDesign);
    GroupISPC{subj} = ztISPC;
end

%%
GroupDesign.RT = round(GroupDesign.RT*1000);
GroupDesign.PC = GroupDesign.PC - 0.33;

% Create CurrentConflict column
GroupDesign.CurrentConflict = repmat({'NoConflict'}, height(GroupDesign), 1); % Default to 'NoConflict'
GroupDesign.CurrentConflict(strcmp(GroupDesign.CurrentType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where CurrentType is 'inc'

% Create PreviousConflict column
GroupDesign.PreviousConflict = repmat({'NoConflict'}, height(GroupDesign), 1); % Default to 'NoConflict'
GroupDesign.PreviousConflict(strcmp(GroupDesign.PreviousType, 'inc')) = {'Conflict'}; % Change to 'Conflict' where PreviousType is 'inc'

%%
GroupDesign.CurrentConflict = categorical(GroupDesign.CurrentConflict);
GroupDesign.CurrentConflict = reordercats(GroupDesign.CurrentConflict, {'NoConflict', 'Conflict'});

GroupDesign.PreviousConflict = categorical(GroupDesign.PreviousConflict);
GroupDesign.PreviousConflict = reordercats(GroupDesign.PreviousConflict, {'NoConflict', 'Conflict'});

GroupDesign.BlockType = categorical(GroupDesign.BlockType);
GroupDesign.BlockType = reordercats(GroupDesign.BlockType, {'minc', 'same', 'mcon'});

%% ISPC by current trial type

trialsidx = GroupDesign([GroupDesign.ElectrodePair] == 1,:);

con = [];
inc = [];

unique_sbjs = unique(trialsidx.Subject);
for i = 1:numel(GroupISPC)
    temptrials = trialsidx([trialsidx.Subject] == unique_sbjs(i),:);
    con(:,i) = squeeze(mean(GroupISPC{i}([temptrials.CurrentConflict] == 'NoConflict',:,:)));
    inc(:,i) = squeeze(mean(GroupISPC{i}([temptrials.CurrentConflict] == 'Conflict',:,:)));
end
% con = con(1:59,:);
% inc = inc(1:59,:);


addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))
[hl, ~] = boundedline(newt,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',newt,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Theta ISPC (z)')
xlabel('Time (s)')
title('dmPFC-dlPFC: Cue Aligned Theta ISPC')
xline(meanRT,'--k','LineWidth',1.25)
yticks([0 1.5 3.5])
xticks([0 0.625 1.25])
xlim([-0.203 1.25])
set(gca,'FontSize',16)
%% Intercept only incongruent

% Shrink length of trial based on RT (S)
IT(60:end) = [];
IB(60:end) = [];
IP(60:end) = [];
IF(60:end) = [];
newt(60:end) = [];

% Shrink length of trial based on RT (R)
% IT = IT(11:71);
% IB = IB(11:71);
% IP = IP(11:71);
% IF = IF(11:71);
% newt = newt(11:71);

threshP = IP < 0.05;
islands = bwconncomp(threshP);

sig_clusters = find(cellfun(@(x) length(x),islands.PixelIdxList) > 5);

for ii = 1:length(sig_clusters)
    sig_idxs{ii} = islands.PixelIdxList{sig_clusters(ii)};
end

% Different way of plotting: plot the t-stats

plot(newt,IF,'r')
hold on
for ii = 1:length(sig_idxs)
    plot(newt(sig_idxs{ii}),IF(sig_idxs{ii}),'r','LineWidth',3)
end
hold off

% title('dlPFC Theta Power Conflict Main Effect')
title('dlPFC Beta Power Conflict Main Effect')

ylabel('F-Statistic')
xlabel('Time (s)')
xlim([-0.2 1.25])
xticks([0 0.625 1.25])
% xlim([-0.7503 0.75]) % response locked
% xticks([-0.5 0 0.5]) % response locked
yticks([0 17 35])
% yticks([0 40 80]) % response locked

RTs = GroupDesignMatrixL([GroupDesignMatrixL.Electrode] == 1,:);
meanRT = mean(RTs.RT)/1000;
% xline(meanRT,'--k')
xline(0,'--k')
set(gca,'FontSize',16)
%% ISPC grand average

for i = 1:numel(files)
    
    % Load the sigChans and the zISPC file
    load(files{i}, 'sigChans')
    load(['/Users/anaskhan/Desktop/PRJ_Stroop/results/newROIs/ISPC/Stim/' files{i}])
    
    % Initialize variable to store extracted ISPC values for significant channels
    sigISPC = [];
    
    % Loop through each significant channel pair in sigChans
    for j = 1:size(sigChans, 1)
        mpfcChan = sigChans{j, 1};
        lpfcChan = sigChans{j, 2};
        
        % Find the indices of the MPFC and LPFC channels
        mpfcIdx = find(strcmp(MPFC, mpfcChan));
        lpfcIdx = find(strcmp(LPFC, lpfcChan));
        
        % If both channels are found in the lists, extract the corresponding zISPC data
        if ~isempty(mpfcIdx) && ~isempty(lpfcIdx)
            % Extract zISPC values for the current significant channel pair
            sigISPC(:,:,j) = squeeze(mean(zISPC(:, mpfcIdx, lpfcIdx, :, :)));
        end
    end

    Alltrials(:,:,i) = mean(sigISPC, 3); 
end

[~,fidx] = arrayfun(@(x) min(abs(x-frex)), [2 4 8 12 30]);
frexticks = frex(fidx);
clear fidx

contourf(newt,frex,mean(Alltrials,3),50,'EdgeColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 30])
ylim([2 30])
xlim([-0.203 1.25])
xticks([0 0.625 1.25])
% xlim([-0.7503 0.75]) % response aligned lims
% xticks([-0.5 0 0.5]) % response aligned lims
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('dACC-dlPFC: Cue Aligned Grand Average')
% title('dACC-dlPFC: Response Aligned ISPC')
cb = colorbar;
clim([-1 1])
lims = clim;
set(cb,'YTick',[lims(1) 0 lims(2)])
cb.Label.String = ['\bf' 'ISPC (z)' '\rm']; 
fontsize(gcf,16,'points')
