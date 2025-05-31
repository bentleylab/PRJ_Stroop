SBJs = {'CP24','IR21','IR26','IR31','IR32','IR35','IR39','IR41',...
        'IR52','IR54','IR57','IR61','IR65','IR67','IR68','IR72','IR74'};

cnt = 1;
for sbj = 1:numel(SBJs)
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJs{sbj} '_vars.m']);
    load(strcat(SBJ_vars.dirs.preproc,SBJs{sbj},'_preproc_',proc_id,'.mat'));

    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJs{sbj},'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);
    load([root_dir 'PRJ_Stroop/data/' SBJs{sbj} '/04_proc/' SBJs{sbj} '_smANOVA_ROI_pCNI_PC_HGh_R5t1_zbt_WL05_WS002.mat'])


    MPFC_chans = fn_select_elec_lab_match(elec, 'b', 'Dx', 'MPFC');
    LPFC_chans = fn_select_elec_lab_match(elec, 'b', 'Dx', 'LPFC');

    if isempty(data.label(ismember(data.label,MPFC_chans) & w2{1}.sig_chans)) | isempty(data.label(ismember(data.label,LPFC_chans) & w2{1}.sig_chans))
        continue
    else
        groupw2MA(:,cnt) = mean(w2{1}.trial(ismember(data.label,MPFC_chans) & w2{1}.sig_chans,:).*100);
        groupw2MI(:,cnt) = mean(w2{1}.trial(ismember(data.label,MPFC_chans) & ~w2{1}.sig_chans,:).*100);
    
        groupw2LA(:,cnt) = mean(w2{1}.trial(ismember(data.label,LPFC_chans) & w2{1}.sig_chans,:).*100);
        groupw2LI(:,cnt) = mean(w2{1}.trial(ismember(data.label,LPFC_chans) & ~w2{1}.sig_chans,:).*100);
        cnt = cnt + 1;
    end

end

meanWMA = mean(groupw2MA,2);
meanWLA = mean(groupw2LA,2);
meanWMI = mean(groupw2MI,2);
meanWLI = mean(groupw2LI,2);

SEMWMA = std(groupw2MA,[],2)./sqrt(size(groupw2MA,2));
SEMWLA = std(groupw2LA,[],2)./sqrt(size(groupw2LA,2));
SEMWMI = std(groupw2MI,[],2)./sqrt(size(groupw2MI,2));
SEMWLI = std(groupw2LI,[],2)./sqrt(size(groupw2LI,2));

t = w2{1}.time';

figure(1)
boundedline(t, meanWMA, SEMWMA, '-b', t, meanWMI, SEMWMI, '-k')
ylabel('[%EV]')
xlabel('Time (s)')
title('MPFC')

figure(2)
boundedline(t, meanWLA, SEMWLA, '-r', t, meanWLI, SEMWLI, '-k')
ylabel('[%EV]')
xlabel('Time (s)')
title('LPFC')

cnt = 1;
for sbj = 1:numel(SBJs)
    eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJs{sbj} '_vars.m']);
    load(strcat(SBJ_vars.dirs.preproc,SBJs{sbj},'_preproc_',proc_id,'.mat'));

    atlas_id = 'Dx';
    load([SBJ_vars.dirs.recon,SBJs{sbj},'_elec_',proc_id,'_pat_',atlas_id,'_final.mat']);
    load([root_dir 'PRJ_Stroop/data/' SBJs{sbj} '/04_proc/' SBJs{sbj} '_smANOVA_ROI_pCNI_PC_HGh_R5t1_zbt_WL05_WS002.mat'])


    MPFC_chans = fn_select_elec_lab_match(elec, 'b', 'Dx', 'MPFC');
    LPFC_chans = fn_select_elec_lab_match(elec, 'b', 'Dx', 'LPFC');

    if isempty(data.label(ismember(data.label,MPFC_chans) & w2{2}.sig_chans)) | isempty(data.label(ismember(data.label,LPFC_chans) & w2{2}.sig_chans))
        continue
    else
        groupw2MA(:,cnt) = mean(w2{2}.trial(ismember(data.label,MPFC_chans) & w2{2}.sig_chans,:).*100);
        groupw2MI(:,cnt) = mean(w2{2}.trial(ismember(data.label,MPFC_chans) & ~w2{2}.sig_chans,:).*100);
    
        groupw2LA(:,cnt) = mean(w2{2}.trial(ismember(data.label,LPFC_chans) & w2{2}.sig_chans,:).*100);
        groupw2LI(:,cnt) = mean(w2{2}.trial(ismember(data.label,LPFC_chans) & ~w2{2}.sig_chans,:).*100);
        cnt = cnt + 1;
    end

end

meanWMA = mean(groupw2MA,2);
meanWLA = mean(groupw2LA,2);
meanWMI = mean(groupw2MI,2);
meanWLI = mean(groupw2LI,2);

SEMWMA = std(groupw2MA,[],2)./sqrt(size(groupw2MA,2));
SEMWLA = std(groupw2LA,[],2)./sqrt(size(groupw2LA,2));
SEMWMI = std(groupw2MI,[],2)./sqrt(size(groupw2MI,2));
SEMWLI = std(groupw2LI,[],2)./sqrt(size(groupw2LI,2));

t = w2{2}.time';

figure(3)
boundedline(t, meanWMA, SEMWMA, '-b', t, meanWMI, SEMWMI, '-k')
ylabel('[%EV]')
xlabel('Time (s)')
title('MPFC')

figure(4)
boundedline(t, meanWLA, SEMWLA, '-r', t, meanWLI, SEMWLI, '-k')
ylabel('[%EV]')
xlabel('Time (s)')
title('LPFC')