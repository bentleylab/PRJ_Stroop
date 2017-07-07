%% RT Comparison
% Compare automatically and manually marked RTs; save png with differences
% Also corrects trial_info variable and saves out corrected and old versions
% New spreadsheet from Rana should be converted into csv with format:
%   block_n, trial_n, error (0/1), resp_onset
close all
clear all

%% Parameters
SBJ = 'IR32';
preproc_id = strcat(SBJ,'_IH_CA');
outlier_thresh = 0.1;                 % in ms
% outlier_sd_thresh = 10;             % in standard deviations

event_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/03_events/');
csv_filename = strcat(event_dir,SBJ,'_Rana_RTs_only.csv');%IR31: '_noB1',IR21: '_notail'
trial_info_filename = strcat(event_dir,preproc_id,'_trial_info.mat');
trial_info_auto_filename = strcat(event_dir,preproc_id,'_trial_info_auto.mat');
trigger_plot_filename = strcat(event_dir,preproc_id,'_event_trigger_plot.fig');

%% Load Data
if exist(trial_info_auto_filename)      % if has been run before, start from the orig/auto data
    load(trial_info_auto_filename);
    trial_info = trial_info_auto;
else
    load(trial_info_filename);
end
man_data = csvread(csv_filename,1,0);

if trial_info.sample_rate~=1000
    error('ERROR!!! sample_rate ~= 1000, so new trial_info sample idxs will be wrong!')
end

%% Calculate RTs, Differences
if ~isempty(trial_info.ignore_trials)
    man_data(trial_info.ignore_trials,:) = [];
    fprintf('Ignoring %i trials\n',length(trial_info.ignore_trials));
end
skip_acc = find(man_data(:,3)==-1);
man_data(skip_acc,3) = NaN;
skip_rt = find(man_data(:,4)==-1);
man_data(skip_rt,4) = NaN;
fprintf('Skipping %i trials for accuracy\n',length(skip_acc));
fprintf('Skipping %i trials for RTs\n',length(skip_rt));
%onset_auto = trial_info.resp_onset;
rt_auto    = trial_info.response_time;
onset_man  = man_data(:,4);
rt_man     = onset_man-trial_info.word_onset./1000;
dif       = rt_auto-rt_man;
dif_mean  = nanmean(dif);
dif_sd    = nanstd(dif);

%% Plot RTs, Differences
figure
subplot(2,1,1); hold on;
plot(rt_auto,'r');
plot(rt_man,'k');
legend('auto','manual');
title(strcat(SBJ,'- RT Comparison'));

subplot(2,1,2);
plot(dif);
legend('auto-man');
title(strcat('Mean=',num2str(dif_mean),' +/- ',num2str(dif_sd)));

%% Find worst discrepancies
% thresh = outlier_sd_thresh*std(dif);
thresh = outlier_thresh;
worst_dif_ix = find(abs(dif)>thresh);
worst_trials = NaN([length(worst_dif_ix),3]);       % block_n, trial_n, dif
worst_trials(:,1) = trial_info.block_n(worst_dif_ix);
block_starts = [1 find(diff(trial_info.block_n)~=0)'+1];
% block_len = unique(diff(block_starts));
if strcmp(SBJ,'IR21')
    block_len = 24;
else
    block_len = 36;
end
worst_trials(:,2) = trial_info.trial_n(worst_dif_ix)-(trial_info.block_n(worst_dif_ix)-1).*block_len;
worst_trials(:,3) = dif(worst_dif_ix);
%openfig(trigger_plot_filename);
% disp(strcat('Threshold =',num2str(outlier_sd_thresh),' SDs = ',num2str(thresh),' ms'));
disp(strcat('Threshold = ',num2str(thresh),' ms'));
if isempty(worst_trials)
    disp('All trials under threshold!');
else
    fprintf('Bad trials detected:\n');
    worst_trials
end

%% Update trial_info, save results
trial_info_auto = trial_info;
trial_info.resp_onset = man_data(:,4).*1000;
trial_info.response_time = rt_man;
trial_info.error = man_data(:,3);

save(trial_info_filename, 'trial_info');
save(trial_info_auto_filename, 'trial_info_auto');