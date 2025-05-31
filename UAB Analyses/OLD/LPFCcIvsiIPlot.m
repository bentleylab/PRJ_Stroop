files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

compcI = [];
compiI = [];

for subj = 1:numel(files)

    load(files{subj},"lpfcDesign","LPFC_pows","newt","frex")

    % Remove rows from design matrix and power matrix that have previous
    % type "None"
    rowsToremove = strcmpi(lpfcDesign.PreviousType,'None');
    shortDesign = lpfcDesign(lpfcDesign.Electrode == 1,:);
    rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

    lpfcDesign(rowsToremove,:) = [];
    shortDesign(rowsToremove_power,:) = [];
    LPFC_pows(rowsToremove_power,:,:,:) = [];

    % Remove rows from design matrix and power matrix that have previous or
    % current type "neu"
    rowsToremove = strcmpi(lpfcDesign.PreviousType,'neu') | strcmpi(lpfcDesign.CurrentType,'neu');
    rowsToremove_power = strcmpi(shortDesign.PreviousType,'neu') | strcmpi(shortDesign.CurrentType,'neu');

    lpfcDesign(rowsToremove,:) = [];
    shortDesign(rowsToremove_power,:) = [];
    LPFC_pows(rowsToremove_power,:,:,:) = [];

    % Finally remove all rows from design matrix and power matrix that are
    % current type "con"
    rowsToremove = strcmpi(lpfcDesign.CurrentType,'con');
    rowsToremove_power = strcmpi(shortDesign.CurrentType,'con');

    lpfcDesign(rowsToremove,:) = [];
    shortDesign(rowsToremove_power,:) = [];
    LPFC_pows(rowsToremove_power,:,:,:) = [];

    compcI(:,:,subj) = squeeze(mean(mean(LPFC_pows(strcmpi(shortDesign.PreviousType,'con'),:,:,:)),2));
    compiI(:,:,subj) = squeeze(mean(mean(LPFC_pows(strcmpi(shortDesign.PreviousType,'inc'),:,:,:)),2));
    diffs(:,:,subj) = squeeze(mean(mean(LPFC_pows(strcmpi(shortDesign.PreviousType,'con'),:,:,:)) ...
        - mean(LPFC_pows(strcmpi(shortDesign.PreviousType,'inc'),:,:,:)),2));
end
%%
[~,fidx] = arrayfun(@(x) min(abs(x-frex)), [2 4 8 12 20 32]);
frexticks = frex(fidx);
clear fidx

% cI
figure(1)
contourf(newt,frex,mean(compcI,3),50,'LineColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 20 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-1 1]);
lims = clim;
set(cb,'YTick',[lims(1) 0 lims(2)])
cb.Label.String = ['\bf' 'Power (z)' '\rm']; 
title('LPFC: Cue Aligned cI')
fontsize(gcf,16,'points')

% iI
figure(2)
contourf(newt,frex,mean(compiI,3),50,'LineColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 20 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([lims(1) lims(2)]);
set(cb,'YTick',[lims(1) 0 lims(2)])
cb.Label.String = ['\bf' 'Power (z)' '\rm']; 
title('LPFC: Cue Aligned iI')
fontsize(gcf,16,'points')

% Diff
figure(3)
contourf(newt,frex,mean(diffs,3),50,'LineColor','none')
set(gca,'YScale','log')
yticks(frexticks)
yticklabels([2 4 8 12 20 32])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
cb = colorbar;
clim([-0.5 0.5]);
set(cb,'YTick',[-0.5 0 0.5])
cb.Label.String = ['\bf' 'Power (z)' '\rm']; 
title('LPFC: Cue Aligned cI - iI')
fontsize(gcf,16,'points')
