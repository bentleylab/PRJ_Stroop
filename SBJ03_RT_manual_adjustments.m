function SBJ03_RT_manual_adjustments(SBJ,block,outlier_thresh,save_plot,save_it)
%% RT Comparison
% Compare automatically and manually marked RTs; save png with differences
% Also corrects trial_info variable and saves out corrected and old versions
% New spreadsheet from Rana should be converted into csv with format:
%   block_n, trial_n, error (0/1), resp_onset
%% Parameters
if isempty(outlier_thresh)
    fprintf('WARNING: No outlier_threshold given; using default = 0.15s\n');
    outlier_thresh = 0.15;                 % in ms
%     outlier_sd_thresh = 10;             % in standard deviations
end

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
if numel(SBJ_vars.raw_file)>1
    block_suffix = strcat('_',SBJ_vars.block_name{block});
else
    block_suffix = SBJ_vars.block_name{block};   % should just be ''
end

csv_filename = strcat(SBJ_vars.dirs.events,SBJ,'_RT_manual',block_suffix,'.csv');
trial_info_man_filename = strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_manual',block_suffix,'.mat');
trial_info_auto_filename = strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_auto',block_suffix,'.mat');

%% Load Data
load(trial_info_auto_filename);
man_data = csvread(csv_filename,1,0);

if trial_info.sample_rate~=1000
    error('ERROR!!! sample_rate ~= 1000, so new trial_info sample idxs will be wrong!')
end

%% Calculate RTs, Differences
% Remove ignore_trials to match length of lists
if ~isempty(trial_info.ignore_trials)
    man_data(trial_info.ignore_trials,:) = [];
    fprintf('Ignoring %i trials\n',length(trial_info.ignore_trials));
end
% Remove trials with artifacts, bad RTs, or errors (keeping list length equal)
skip_bad = find(man_data(:,3)==-1);
skip_err = find(man_data(:,3)==1);
skip_rt  = find(man_data(:,4)==-1);

skip_trial_ix = unique(vertcat(skip_bad,skip_rt,skip_err));
resp_onset_man = man_data(:,4);
resp_onset_man(skip_trial_ix) = NaN;
fprintf('Skipping %i trials for artifacts\n',length(skip_bad));
fprintf('Skipping %i trials for errors\n',length(skip_err));
fprintf('Skipping %i trials for no RTs\n',length(skip_rt));
fprintf('Trials remaining = %i (tossed %i)\n\n',sum(~isnan(resp_onset_man)),length(skip_trial_ix));

% Compute the difference in seconds
rt_auto     = trial_info.response_time;
rt_man      = resp_onset_man-trial_info.word_onset./1000;
dif         = rt_auto-rt_man;
dif_mean    = nanmean(dif);
dif_sd      = nanstd(dif);

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

if save_plot
    saveas(gcf, strcat(SBJ_vars.dirs.events,SBJ,'_RT_comparison',block_suffix,'.png'));
end

%% Identify worst discrepancies
% Find worst trials
% thresh = outlier_sd_thresh*std(dif);
thresh = outlier_thresh;
worst_dif_ix = find(abs(dif)>thresh);

% Build info on worst trials
worst_trials = NaN([length(worst_dif_ix),4]);       % block_n, trial_n, dif
worst_trials(:,1) = trial_info.block_n(worst_dif_ix);
worst_trials(:,3) = trial_info.trial_n(worst_dif_ix);
block_starts = [1 find(diff(trial_info.block_n)~=0)'+1];
% block_len = unique(diff(block_starts));
if strcmp(SBJ,'IR21')   % adjust for original paradigm design in IR21 only (IR26,27?)
    block_len = 24;
else
    block_len = 36;
end
worst_trials(:,2) = trial_info.trial_n(worst_dif_ix)-(trial_info.block_n(worst_dif_ix)-1).*block_len;
worst_trials(:,4) = dif(worst_dif_ix);
%openfig(trigger_plot_filename);
% disp(strcat('Threshold =',num2str(outlier_sd_thresh),' SDs = ',num2str(thresh),' ms'));
disp(strcat('Threshold = ',num2str(thresh),' ms'));
if isempty(worst_trials)
    disp('All trials under threshold!');
else
    fprintf('Bad trials detected: (B#, T# per B, total T#, RT difference)\n');
    worst_trials
end

%% Update trial_info, save results
% Add all manual data to trial_info (without NaNs to keep all info for later)
trial_info.resp_onset    = man_data(:,4).*1000;
trial_info.response_time = man_data(:,4)-trial_info.word_onset./1000;
trial_info.error         = man_data(:,3);

% Keep -1 RT info
trial_info.resp_onset(skip_rt)    = -1;
trial_info.response_time(skip_rt) = -1;

if save_it
    save(trial_info_man_filename, 'trial_info');
end

end
