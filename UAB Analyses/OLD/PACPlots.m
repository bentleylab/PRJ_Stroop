files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

GAPACML = [];
for subj = 1:length(files)
    load(files{subj})
    GAPACML(:,:,subj) = mean(cat(3,PACMtLALL.zMI{:}),3);
end


GAPACLM = [];
for subj = 1:length(files)
    load(files{subj})
    GAPACLM(:,:,subj) = mean(cat(3,PACLtMALL.zMI{:}),3);
end

LF_steps = 2:2:12;

HF_steps = 30:5:150;

% Grand Average zPAC M->L
contourf(LF_steps,HF_steps,mean(GAPACML,3)',50,'EdgeColor','none')
colorbar
clim([0 2])

% Grand Average zPAC L->M
contourf(LF_steps,HF_steps,mean(GAPACLM,3)',50,'EdgeColor','none')
colorbar
clim([0 2])
%%
groupPAC = [];
for subj = 1:length(files)
    load(files{subj})
    zPAC = cat(3,PACMtLALL.zMI{:});
    zPAC = squeeze(mean(mean(zPAC(find(LF_steps==2):find(LF_steps==8),find(HF_steps==70):find(HF_steps==140),:),1),2));
    TG_chans = zPAC > 1.64;
    PACMtLALL.zMI(~TG_chans) = [];
    PACMtLALL.MI(~TG_chans) = [];
    PACMtLALL.PPHist(~TG_chans) = [];
    groupPAC{subj}.data = PACMtLALL;
end
%%
groupPAC = [];
for subj = 1:length(files)
    load(files{subj})
    zPAC = cat(3,PACLtMALL.zMI{:});
    zPAC = squeeze(mean(mean(zPAC(find(LF_steps==2):find(LF_steps==8),find(HF_steps==70):find(HF_steps==140),:),1),2));
    TG_chans = zPAC > 1.64;
    PACLtMALL.zMI(~TG_chans) = [];
    PACLtMALL.MI(~TG_chans) = [];
    PACLtMALL.PPHist(~TG_chans) = [];
    groupPAC{subj}.data = PACLtMALL;
end
%%
%MPFC -> LPFC
subplot(1,2,1)
allPac = mean(cat(3,PACMtL{1,1}.zMI{:}),3);
contourf(LF_steps,HF_steps,allPac',50,'EdgeColor','none')
colorbar
clim([0 2])
title('MC Block: MPFC Phase Modulating LPFC Amp')
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')

subplot(1,2,2)
allPac = mean(cat(3,PACMtL{1,2}.zMI{:}),3);
contourf(LF_steps,HF_steps,allPac',50,'EdgeColor','none')
colorbar
clim([0 2])
title('MI Block: MPFC Phase Modulating LPFC Amp')
xlabel('Amplitude Frequency (Hz)')

figure(2)
%LPFC -> MPFC
subplot(1,2,1)
allPac = mean(cat(3,PACLtM{1,1}.zMI{:}),3);
contourf(LF_steps,HF_steps,allPac',50,'EdgeColor','none')
colorbar
clim([0 2])
title('MC Block: LPFC Phase Modulating MPFC Amp')
xlabel('Phase Frequency (Hz)')
ylabel('Amplitude Frequency (Hz)')

subplot(1,2,2)
allPac = mean(cat(3,PACLtM{1,2}.zMI{:}),3);
contourf(LF_steps,HF_steps,allPac',50,'EdgeColor','none')
colorbar
clim([0 2])
title('MI Block: LPFC Phase Modulating MPFC Amp')
xlabel('Amplitude Frequency (Hz)')