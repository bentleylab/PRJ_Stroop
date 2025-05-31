files = dir;
files = {files(~(strcmpi('.',{files.name}) | strcmpi('..',{files.name}) | strcmpi('.DS_Store',{files.name}))).name};
%% Stim
for subj = 1:length(files)
    load(files{subj})
    net_incLtM = granger_inc(:,:,2,1)-granger_inc(:,:,1,2);
    net_incLtM = net_incLtM(:,nearest(time,-0.2):nearest(time,1.25));
    
    net_conLtM = granger_con(:,:,2,1)-granger_con(:,:,1,2);
    net_conLtM = net_conLtM(:,nearest(time,-0.2):nearest(time,1.25));

    con(:,:,subj) = net_conLtM;
    inc(:,:,subj) = net_incLtM;
end
time = time(nearest(time,-0.2):nearest(time,1.25));

contourf(time,frex,mean(inc-con,3),50,'EdgeColor','none')
colorbar
clim([-0.02 0.02])
ylim([4 30])
set(gca,'YScale','log')
xticks([0 0.625 1.25])
yticks([4 8 12 30])
cb = colorbar;
clim([-0.02 0.02])
lims = clim;
set(cb,'YTick',[lims(1) 0 lims(2)])
cb.Label.String = ['\bf' 'Net dlPFC -> dACC GC' '\rm']; 
title('Incongruent - Congruent: Net dlPFC -> dACC GC')
fontsize(gcf,24,'points')

%% Stim new labels
for subj = 1:length(files)
    load(files{subj})
    net_incLtM = granger_conflict(:,:,2,1)-granger_conflict(:,:,1,2);
    net_incLtM = net_incLtM(:,nearest(time,-0.2):nearest(time,1.25));
    
    net_conLtM = granger_NoConflict(:,:,2,1)-granger_NoConflict(:,:,1,2);
    net_conLtM = net_conLtM(:,nearest(time,-0.2):nearest(time,1.25));

    con(:,:,subj) = net_conLtM;
    inc(:,:,subj) = net_incLtM;
end
time = time(nearest(time,-0.2):nearest(time,1.25));

contourf(time,frex,mean(inc-con,3),50,'EdgeColor','none')
colorbar
ylim([4 30])
set(gca,'YScale','log')
xticks([0 0.625 1.25])
yticks([4 8 12 30])
cb = colorbar;
clim([-0.0125 0.0125])
lims = clim;
set(cb,'YTick',[-0.01 0 0.01])
cb.Label.String = ['\bf' 'Net dlPFC -> dACC GC' '\rm']; 
title('Conflict - NoConflict: Net dlPFC -> dACC GC')
set(gca,'FontSize',16)
%% Conflict x Block Plots
% Initialize arrays to store data across subjects

mcon = [];
minc = [];
same = [];

for subj = 1:length(files)
    load(files{subj}); % Load subject-specific data
    
    % Extract time indices
    time_idx = nearest(time,-0.2):nearest(time,1.25);
    
    % Compute net GC for Conflict and NoConflict conditions for 'mcon' block type
    net_conflict_mcon = granger_conflict_mcon(:,:,2,1) - granger_conflict_mcon(:,:,1,2);
    net_noconflict_mcon = granger_NoConflict_mcon(:,:,2,1) - granger_NoConflict_mcon(:,:,1,2);
    net_diff_mcon = net_conflict_mcon(:, time_idx) - net_noconflict_mcon(:, time_idx);
    mcon(:,:,subj) = net_diff_mcon;
    
    % Compute net GC for 'minc' block type
    net_conflict_minc = granger_conflict_minc(:,:,2,1) - granger_conflict_minc(:,:,1,2);
    net_noconflict_minc = granger_NoConflict_minc(:,:,2,1) - granger_NoConflict_minc(:,:,1,2);
    net_diff_minc = net_conflict_minc(:, time_idx) - net_noconflict_minc(:, time_idx);
    minc(:,:,subj) = net_diff_minc;
    
    % Compute net GC for 'same' block type
    net_conflict_same = granger_conflict_same(:,:,2,1) - granger_conflict_same(:,:,1,2);
    net_noconflict_same = granger_NoConflict_same(:,:,2,1) - granger_NoConflict_same(:,:,1,2);
    net_diff_same = net_conflict_same(:, time_idx) - net_noconflict_same(:, time_idx);
    same(:,:,subj) = net_diff_same;
end

% Update 'time' vector to match the selected time range
time = time(nearest(time,-0.2):nearest(time,1.25));

% Compute the mean across subjects for each block type
mean_mcon = mean(mcon, 3);
mean_minc = mean(minc, 3);
mean_same = mean(same, 3);

% Plotting the net GC difference (Conflict - NoConflict) for each block type
block_types = {'mcon', 'minc', 'same'};
mean_data = {mean_mcon, mean_minc, mean_same};

for i = 1:length(block_types)
    figure;
    contourf(time, frex, mean_data{i}, 50, 'EdgeColor', 'none');
    colorbar;
    ylim([4 30]);
    set(gca, 'YScale', 'log');
    xticks([0 0.625 1.25]);
    yticks([4 8 12 30]);
    cb = colorbar;
    clim([-0.0125 0.0125]);
    set(cb, 'YTick', [-0.01 0 0.01]);
    cb.Label.String = '\bfNet dlPFC -> dACC GC\rm';
    title(['Conflict - NoConflict: ', block_types{i}, ' Net dlPFC -> dACC GC']);
    set(gca, 'FontSize', 16);
end


%% Resp
for subj = 1:length(files)
    load(files{subj})
    net_incLtM = granger_inc(:,:,2,1)-granger_inc(:,:,1,2);
    net_incLtM = net_incLtM(:,nearest(time,-0.75):nearest(time,0.75));
    
    net_conLtM = granger_con(:,:,2,1)-granger_con(:,:,1,2);
    net_conLtM = net_conLtM(:,nearest(time,-0.75):nearest(time,0.75));

    con(:,:,subj) = net_conLtM;
    inc(:,:,subj) = net_incLtM;
end
time = time(nearest(time,-0.75):nearest(time,0.75));

contourf(time,frex,mean(inc-con,3),50,'EdgeColor','none')
colorbar
clim([-0.025 0.025])
ylim([4 30])
set(gca,'YScale','log')
cb = colorbar;
clim([-0.025 0.025])
ylim([4 30])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
set(cb,'yTick',[-0.025 0 0.025])
cb.Label.String = ['\bf' 'Net dlPFC -> dACC GC' '\rm'];
%% Theta Resp
for subj = 1:length(files)
    load(files{subj})
    net_incLtM = granger_inc(:,:,2,1)-granger_inc(:,:,1,2);
    net_incLtM = net_incLtM(:,nearest(time,-0.75):nearest(time,0.75));
    
    net_conLtM = granger_con(:,:,2,1)-granger_con(:,:,1,2);
    net_conLtM = net_conLtM(:,nearest(time,-0.75):nearest(time,0.75));

    con(:,subj) = mean(net_conLtM(frex >= 4 & frex <= 8,:));
    inc(:,subj) = mean(net_incLtM(frex >= 4 & frex <= 8,:));
end
time = time(nearest(time,-0.75):nearest(time,0.75));

addpath(genpath('/Users/anaskhan/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/boundedline.m'))

[hl, ~] = boundedline(time,mean(con,2),std(con,[],2)./sqrt(size(con,2)),'b',time,mean(inc,2),std(inc,[],2)./sqrt(size(inc,2)),'r');
hl(1).LineWidth = 1.25;
hl(2).LineWidth = 1.25;
ylabel('Net dlPFC -> dACC Theta GC')
xlabel('Time (s)')


%%
for subj = 1:length(files)
    load(files{subj})
    net_iCLtM = granger_iC(:,:,2,1)-granger_iC(:,:,1,2);
    net_iCLtM = net_iCLtM(:,dsearchn(time',-0.2):dsearchn(time',1.25));
    
    net_cCLtM = granger_cC(:,:,2,1)-granger_cC(:,:,1,2);
    net_cCLtM = net_cCLtM(:,dsearchn(time',-0.2):dsearchn(time',1.25));

    cC(:,:,subj) = net_cCLtM;
    iC(:,:,subj) = net_iCLtM;
end
time = time(dsearchn(time',-0.2):dsearchn(time',1.25));

contourf(time,frex,mean(iC-cC,3),50,'EdgeColor','none')
set(gca,'YScale','log')
cb = colorbar;
clim([-0.05 0.05])
ylim([4 30])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
set(cb,'yTick',[-0.05 0 0.05])
cb.Label.String = ['\bf' 'Net dlPFC -> dACC GC' '\rm'];