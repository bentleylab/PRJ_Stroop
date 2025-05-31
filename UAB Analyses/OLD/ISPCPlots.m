files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

%%
allISPC = [];

for i = 1:numel(files)
    load(files{i})
    allISPC(:,:,i) = squeeze(mean(mean(mean(zISPC),3),2));
end

[~,fidx] = arrayfun(@(x) min(abs(x-frex)), [2 8 12 30]);
frexticks = frex(fidx);
clear fidx

contourf(newt,frex,mean(allISPC,3),50,'EdgeColor','none')
set(gca,'YScale','log')
cb = colorbar;
clim([-1 1])
set(cb,'yTick',[-1 0 1])
cb.Label.String = ['\bf' 'ISPC (z)' '\rm'];
yticks(frexticks)
yticklabels([2 8 12 30])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
title('dACC-dlPFC Phase Coherence Grand Average')
%%

cC = [];
iC = [];

% for i = 1:numel(files)
%     load(files{i})
%     zISPC = squeeze(mean(zISPC(:,:,:,thetafrex,:),4));
%     shortDesign = coherenceDesign([coherenceDesign.ElectrodePair] == 1,:);
%     zISPC(strcmpi(shortDesign.PreviousType,'None'),:,:,:) = [];
%     shortDesign(strcmpi(shortDesign.PreviousType,'None'),:) = [];
%     zISPC(strcmpi(shortDesign.PreviousType,'neu'),:,:,:) = [];
%     shortDesign(strcmpi(shortDesign.PreviousType,'neu'),:) = [];
%     zISPC(~strcmpi(shortDesign.CurrentType,'con'),:,:,:) = [];
%     shortDesign(~strcmpi(shortDesign.CurrentType,'con'),:) = [];
% end

GroupISPC = [];
GroupDesign = [];

for i = 1:numel(files)
    load(files{i})
    thetafrex = frex >= 2 & frex <= 6;
    zISPC = mean(zISPC(:,:,:,thetafrex,:),4);
    zISPC = squeeze(mean(mean(zISPC,3),2));
    shortDesign = coherenceDesign([coherenceDesign.ElectrodePair] == 1,:);
    zISPC(strcmpi(shortDesign.PreviousType,'None'),:) = [];
    shortDesign(strcmpi(shortDesign.PreviousType,'None'),:) = [];
    zISPC(strcmpi(shortDesign.PreviousType,'neu'),:) = [];
    shortDesign(strcmpi(shortDesign.PreviousType,'neu'),:) = [];
    zISPC(~strcmpi(shortDesign.CurrentType,'con'),:) = [];
    shortDesign(~strcmpi(shortDesign.CurrentType,'con'),:) = [];
    GroupDesign = [GroupDesign;shortDesign];
    GroupISPC = [GroupISPC;zISPC];
    cC(:,i) = mean(zISPC(strcmpi(shortDesign.PreviousType,'con'),:));
    iC(:,i) = mean(zISPC(strcmpi(shortDesign.PreviousType,'inc'),:));
end

%%
GroupDesign.PreviousType = categorical(GroupDesign.PreviousType);
GroupDesign.PreviousType = reordercats(GroupDesign.PreviousType, {'con', 'inc'});
GroupDesign.RT = GroupDesign.RT.*1000;
%% Prev Type Stats
[~,numTimePoints] = size(GroupISPC);

previousB = NaN(numTimePoints,1);
previousP = NaN(numTimePoints,1);

for timei = 1:numTimePoints
        
        % Add these power values to the design matrix
        tempDesignMatrix = GroupDesign;
        tempDesignMatrix.ispc = GroupISPC(:,timei);
        
        lme = fitlme(tempDesignMatrix, 'ispc ~ PreviousType + (1 | Subject)', 'FitMethod', 'REML');
        aov = anova(lme,'DFMethod','Satterthwaite');

        previousP(timei) = aov.pValue(strcmpi(aov.Term,'PreviousType'));
        previousB(timei) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'PreviousType_inc'));
end

[h,~,~,~] = fdr_bh(previousP);
plot(newt,previousB)

%% RT Stats
[~,numTimePoints] = size(GroupISPC);

RTB = NaN(numTimePoints,1);
RTP = NaN(numTimePoints,1);

for timei = 1:numTimePoints
        
        % Add these power values to the design matrix
        tempDesignMatrix = GroupDesign;
        tempDesignMatrix.ispc = GroupISPC(:,timei);
        
        lme = fitlme(tempDesignMatrix, 'ispc ~ RT + (1 | Subject)', 'FitMethod', 'REML');
        aov = anova(lme,'DFMethod','Satterthwaite');

        RTP(timei) = aov.pValue(strcmpi(aov.Term,'RT'));
        RTB(timei) = lme.Coefficients.Estimate(strcmpi(lme.Coefficients.Name,'RT'));
end

[h,~,~,~] = fdr_bh(RTP);
plot(newt,RTB)
hold on
y = NaN(size(newt)); 
y(h) = -0.0025; 
plot(newt, y, '-k','LineWidth',5); 
ylabel('RT Coefficient (beta)')
xlabel('Time (s)')
title('ispc ~ RT + (1 | Subject)')
%% Plot cC vs iC phase coherence

ccSEM = std(cC,[],2)./sqrt(size(cC,2));
icSEM = std(iC,[],2)./sqrt(size(iC,2));

[hl,hp] = boundedline(newt,mean(iC,2),icSEM,'r',newt,mean(cC,2),ccSEM,'b');
xlim([-0.2 2])
hl(1).LineWidth = 1.5;
hl(2).LineWidth = 1.5;
legend({'iC','cC'})
xlabel('Time (s)')
ylabel('2-6 Hz ISPC (z)')
%% Granger
cCMtL = [];
cCLtM = [];

iCMtL = [];
iCLtM = [];

% for i = 1:numel(files)
%     load(files{i})
%     zISPC = squeeze(mean(zISPC(:,:,:,thetafrex,:),4));
%     shortDesign = coherenceDesign([coherenceDesign.ElectrodePair] == 1,:);
%     zISPC(strcmpi(shortDesign.PreviousType,'None'),:,:,:) = [];
%     shortDesign(strcmpi(shortDesign.PreviousType,'None'),:) = [];
%     zISPC(strcmpi(shortDesign.PreviousType,'neu'),:,:,:) = [];
%     shortDesign(strcmpi(shortDesign.PreviousType,'neu'),:) = [];
%     zISPC(~strcmpi(shortDesign.CurrentType,'con'),:,:,:) = [];
%     shortDesign(~strcmpi(shortDesign.CurrentType,'con'),:) = [];
% end

% GroupGrangerMtL = [];
% GroupGrangerLtM = [];
% GroupDesign = [];

for i = 1:numel(files)
    load(files{i})

    m_to_l_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'None')) = [];
    l_to_m_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'None')) = [];
    coherenceDesign(strcmpi(coherenceDesign.PreviousType,'None'),:) = [];

    m_to_l_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'neu')) = [];
    l_to_m_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'neu')) = [];
    coherenceDesign(strcmpi(coherenceDesign.PreviousType,'neu'),:) = [];

    m_to_l_GC(:,:,:,~strcmpi(coherenceDesign.CurrentType,'con')) = [];
    l_to_m_GC(:,:,:,~strcmpi(coherenceDesign.CurrentType,'con')) = [];
    coherenceDesign(~strcmpi(coherenceDesign.CurrentType,'con'),:) = [];

    % GroupDesign = [GroupDesign;coherenceDesign];
    % GroupISPC = [GroupISPC;zISPC];
    cCMtL(:,i) = squeeze(mean(mean(mean(m_to_l_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'con')),4),2),1));
    iCMtL(:,i) = squeeze(mean(mean(mean(m_to_l_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'inc')),4),2),1));
    cCLtM(:,i) = squeeze(mean(mean(mean(l_to_m_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'con')),4),2),1));
    iCLtM(:,i) = squeeze(mean(mean(mean(l_to_m_GC(:,:,:,strcmpi(coherenceDesign.PreviousType,'inc')),4),2),1));
end

cCMtLSEM = std(cCMtL,[],2)./sqrt(size(cCMtL,2));
iCMtLSEM = std(iCMtL,[],2)./sqrt(size(iCMtL,2));

cCLtMSEM = std(cCLtM,[],2)./sqrt(size(cCLtM,2));
iCLtMSEM = std(iCLtM,[],2)./sqrt(size(iCLtM,2));

[hl,hp] = boundedline(frex,mean(cCMtL,2),cCMtLSEM,'r',frex,mean(cCLtM,2),cCLtMSEM,'k');
hl(1).LineWidth = 1.5;
hl(2).LineWidth = 1.5;
xlabel('Frequency (Hz)')
ylabel('GC')
legend({'cCMtL','cCLtM'})

figure(2)
[hl,hp] = boundedline(frex,mean(iCMtL,2),iCMtLSEM,'r',frex,mean(iCLtM,2),iCLtMSEM,'k');
hl(1).LineWidth = 1.5;
hl(2).LineWidth = 1.5;
xlabel('Frequency (Hz)')
ylabel('GC')
legend({'iCMtL','iCLtM'})

figure(3)
[hl,hp] = boundedline(frex,mean(cCMtL,2),cCMtLSEM,'r',frex,mean(iCMtL,2),iCMtLSEM,'k');
hl(1).LineWidth = 1.5;
hl(2).LineWidth = 1.5;
xlabel('Frequency (Hz)')
ylabel('GC')
legend({'cCMtL','iCMtL'})

figure(4)
[hl,hp] = boundedline(frex,mean(cCLtM,2),cCLtMSEM,'r',frex,mean(iCLtM,2),iCLtMSEM,'k');
hl(1).LineWidth = 1.5;
hl(2).LineWidth = 1.5;
xlabel('Frequency (Hz)')
ylabel('GC')
legend({'cCLtM','iCLtM'})