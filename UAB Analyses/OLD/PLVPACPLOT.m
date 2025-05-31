files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};

comp_con = [];
comp_ncon = [];
for i = 1:length(files)
    load(files{i})
    comp_con(:,:,i) = conflict;
    comp_ncon(:,:,i) = noconflict;
end

contourf(time,lowFreqs,mean(comp_con-comp_ncon,3),50,'EdgeColor','none')

contourf(time,lowFreqs,mean(comp_con,3),50,'EdgeColor','none')

contourf(time,lowFreqs,mean(comp_ncon,3),50,'EdgeColor','none')

for ti = 1:size(comp_con,2)
    for fi = 1:size(comp_con,1)
        [~, ~, ~, stats] = ttest(comp_con(fi,ti,:),comp_ncon(fi,ti,:));
        t(fi,ti) = stats.tstat;
    end
end

comp_tcon = squeeze(mean(comp_con(2:6,:,:)));
comp_tncon = squeeze(mean(comp_ncon(2:6,:,:)));

[hl, ~] = boundedline(time,mean(comp_tcon,2),std(comp_tcon,[],2)./sqrt(size(comp_tcon,2)),'r',time,mean(comp_tncon,2),std(comp_tncon,[],2)./sqrt(size(comp_tncon,2)),'b')
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('PAC (% Baseline)')
xlabel('Time (s)');
xlim([-0.1 1.25])
title('dmPFC\theta - dlPFC\gamma PAC')
legend({'','','Conflict','NoConflict'})