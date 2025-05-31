%% Import data lateral -> medial - medial -> lateral (using the theta only pipeline)
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_cI = [];
comp_iI = [];
for filei = 1:length(files)
    load(files{filei})
    comp_cI = cat(3,comp_cI,squeeze(l_to_m_GC_cI(3:33,:)-m_to_l_GC_cI(3:33,:)));
    comp_iI = cat(3,comp_iI,squeeze(l_to_m_GC_iI(3:33,:)-m_to_l_GC_iI(3:33,:)));
end

%% Import data lateral -> medial - medial -> lateral
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_cI = [];
comp_iI = [];
for filei = 1:length(files)
    load(files{filei})
    comp_cI = cat(3,comp_cI,squeeze(mean(mean(l_to_m_GC_cI(:,:,3:33,:)-m_to_l_GC_cI(:,:,3:33,:),2),1)));
    comp_iI = cat(3,comp_iI,squeeze(mean(mean(l_to_m_GC_iI(:,:,3:33,:)-m_to_l_GC_iI(:,:,3:33,:),2),1)));
end
%% Import data lateral -> medial
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_cI = [];
comp_iI = [];
for filei = 1:length(files)
    load(files{filei})
    comp_cI = cat(3,comp_cI,squeeze(mean(mean(l_to_m_GC_cI(:,:,3:33,:),2),1)));
    comp_iI = cat(3,comp_iI,squeeze(mean(mean(l_to_m_GC_iI(:,:,3:33,:),2),1)));
end

%% Medial -> lateral data
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_cI = [];
comp_iI = [];
for filei = 1:length(files)
    load(files{filei})
    comp_cI = cat(3,comp_cI,squeeze(mean(mean(m_to_l_GC_cI(:,:,3:33,:),2),1)));
    comp_iI = cat(3,comp_iI,squeeze(mean(mean(m_to_l_GC_iI(:,:,3:33,:),2),1)));
end
%% Import and normalize data
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_cI = [];
comp_iI = [];
for filei = 1:length(files)
    load(files{filei})
    BL = mean(cat(4,l_to_m_GC_cI(:,:,:,1:13),l_to_m_GC_iI(:,:,:,1:13)),4);

    normcI = 100*((l_to_m_GC_cI - BL)./BL);
    normiI = 100*((l_to_m_GC_iI - BL)./BL);

    comp_cI = cat(3,comp_cI,squeeze(mean(mean(normcI(:,:,3:33,:),2),1)));
    comp_iI = cat(3,comp_iI,squeeze(mean(mean(normiI(:,:,3:33,:),2),1)));
end
%% Collapsed across theta
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_cI = [];
comp_iI = [];
for filei = 1:length(files)
    load(files{filei})
    comp_cI = cat(2,comp_cI,squeeze(mean(mean(sum(l_to_m_GC_cI(:,:,5:9,:),3),2),1)));
    comp_iI = cat(2,comp_iI,squeeze(mean(mean(sum(l_to_m_GC_iI(:,:,5:9,:),3),2),1)));
end


%% Collapsed across theta and normalized
files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_cI = [];
comp_iI = [];
for filei = 1:length(files)
    load(files{filei})
    
    thetacI = sum(l_to_m_GC_cI(:,:,5:9,:),3);
    thetaiI = sum(l_to_m_GC_iI(:,:,5:9,:),3);

    BL = mean(cat(4,thetacI(:,:,:,1:13),thetaiI(:,:,:,1:13)),4);

    normcI = 100*((thetacI - BL)./BL);
    normiI = 100*((thetaiI - BL)./BL);

    comp_cI = cat(4,comp_cI,squeeze(mean(mean(normcI,2),1)));
    comp_iI = cat(4,comp_iI,squeeze(mean(mean(normiI,2),1)));

end

comp_cI = squeeze(comp_cI);
comp_iI = squeeze(comp_iI);
%% Plot theta granger

meaniI = mean(comp_iI,2)';
meancI = mean(comp_cI,2)';

semcI = (std(comp_cI,[],2)./sqrt(size(comp_cI,2)))';
semiI = (std(comp_iI,[],2)./sqrt(size(comp_iI,2)))';

% Must be column vectors
plot(grangerTime,meaniI,'r','LineWidth',2)
patch([grangerTime fliplr(grangerTime)],[meaniI-semiI fliplr(meaniI+semiI)],'r','FaceAlpha',0.35,'EdgeColor','none')
hold on
plot(grangerTime,meancI,'b','LineWidth',2)
patch([grangerTime fliplr(grangerTime)],[meancI-semcI fliplr(meancI+semcI)],'b','FaceAlpha',0.35,'EdgeColor','none')

xlabel('Time (s)')
ylabel('Theta Granger Causality')
legend('iI','','cI','','')
%% Garbage stats attempt

for fi = 1:size(comp_cI,1)
    for ti = 1:size(comp_cI,2)
        [~, p(fi,ti)] = ttest2(squeeze(comp_cI(fi,ti,:)),squeeze(comp_iI(fi,ti,:)));
    end
end

%% Plot raw GC lateral to medial

figure(1)
contourf(grangerTime,frex(3:33),mean(comp_cI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-0.04 0.04])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-0.04 0.04])
cb.Label.String = ['\bf' 'LPFC -> MPFC GC' '\rm'];
title('cI: Granger Causality')

figure(2)
contourf(grangerTime,frex(3:33),mean(comp_iI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-0.04 0.04])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-0.04 0.04])
cb.Label.String = ['\bf' 'LPFC -> MPFC GC' '\rm'];
title('iI:  Granger Causality')

figure(3)
contourf(grangerTime,frex(3:33),mean(comp_iI - comp_cI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-0.02 0.02])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-0.02 0.02])
cb.Label.String = ['\bf' 'LPFC -> MPFC GC' '\rm'];
title('iI - cI:  Granger Causality')


%% Plot raw GC medial to lateral

figure(1)
contourf(grangerTime,frex(3:33),mean(comp_cI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-0.04 0.04])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-0.04 0.04])
cb.Label.String = ['\bf' 'MPFC -> LPFC GC' '\rm'];
title('cI: Granger Causality')

figure(2)
contourf(grangerTime,frex(3:33),mean(comp_iI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-0.04 0.04])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-0.04 0.04])
cb.Label.String = ['\bf' 'MPFC -> LPFC GC' '\rm'];
title('iI:  Granger Causality')

figure(3)
contourf(grangerTime,frex(3:33),mean(comp_iI - comp_cI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-0.005 0.005])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-0.005 0.005])
cb.Label.String = ['\bf' 'MPFC -> LPFC GC' '\rm'];
title('iI - cI:  Granger Causality')
%% Plot normalized GC

figure(1)
contourf(grangerTime,frex(3:33),mean(comp_cI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-80 80])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-80 80])
cb.Label.String = ['\bf' 'LPFC -> MPFC GC' '\rm'];
title('cI: Granger Causality')

figure(2)
contourf(grangerTime,frex(3:33),mean(comp_iI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-80 80])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-80 80])
cb.Label.String = ['\bf' 'LPFC -> MPFC GC' '\rm'];
title('iI:  Granger Causality')

figure(3)
contourf(grangerTime,frex(3:33),mean(comp_iI - comp_cI,3),50,'EdgeColor','none')
cb = colorbar;
clim([-40 40])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(cb,'YTick',[-40 40])
cb.Label.String = ['\bf' 'LPFC -> MPFC GC' '\rm'];
title('iI - cI:  Granger Causality')