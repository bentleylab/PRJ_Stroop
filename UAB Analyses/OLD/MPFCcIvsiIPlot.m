files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};


compcI = [];
compiI = [];

for subj = 1:numel(files)

    load(files{subj},"mpfcDesign","MPFC_pows","newt","frex")

    % Remove rows from design matrix and power matrix that have previous
    % type "None"
    rowsToremove = strcmpi(mpfcDesign.PreviousType,'None');
    shortDesign = mpfcDesign(mpfcDesign.Electrode == 1,:);
    rowsToremove_power = strcmpi(shortDesign.PreviousType,'None');

    mpfcDesign(rowsToremove,:) = [];
    shortDesign(rowsToremove_power,:) = [];
    MPFC_pows(rowsToremove_power,:,:,:) = [];

    % Remove rows from design matrix and power matrix that have previous or
    % current type "neu"
    rowsToremove = strcmpi(mpfcDesign.PreviousType,'neu') | strcmpi(mpfcDesign.CurrentType,'neu');
    rowsToremove_power = strcmpi(shortDesign.PreviousType,'neu') | strcmpi(shortDesign.CurrentType,'neu');

    mpfcDesign(rowsToremove,:) = [];
    shortDesign(rowsToremove_power,:) = [];
    MPFC_pows(rowsToremove_power,:,:,:) = [];

    % Finally remove all rows from design matrix and power matrix that are
    % current type "con"
    rowsToremove = strcmpi(mpfcDesign.CurrentType,'con');
    rowsToremove_power = strcmpi(shortDesign.CurrentType,'con');

    mpfcDesign(rowsToremove,:) = [];
    shortDesign(rowsToremove_power,:) = [];
    MPFC_pows(rowsToremove_power,:,:,:) = [];

    compcI(:,:,subj) = squeeze(mean(mean(MPFC_pows(strcmpi(shortDesign.PreviousType,'con'),:,:,:)),2));
    compiI(:,:,subj) = squeeze(mean(mean(MPFC_pows(strcmpi(shortDesign.PreviousType,'inc'),:,:,:)),2));
    diffs(:,:,subj) = squeeze(mean(mean(MPFC_pows(strcmpi(shortDesign.PreviousType,'con'),:,:,:)) ...
        - mean(MPFC_pows(strcmpi(shortDesign.PreviousType,'inc'),:,:,:)),2));
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
title('MPFC: Cue Aligned cI')
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
title('MPFC: Cue Aligned iI')
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
title('MPFC: Cue Aligned cI - iI')
fontsize(gcf,16,'points')