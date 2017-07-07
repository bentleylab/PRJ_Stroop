%% Save Data with only clean trials
clear all
SBJ       = 'IR35';
data_id   = strcat(SBJ,'_LAC_WM');
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
epoch_id  = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];

SBJ_dir    = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
epoch_dir  = strcat(SBJ_dir,'06_epochs/');

%% Load data
load(strcat(SBJ_dir,'04_proc/',data_id,'.mat'));
% need only ok_epochs index, the rest is all in 04_proc file (variables seem same as of 10/10/16...)
load(strcat(epoch_dir,data_id,'_',epoch_id,'.mat'),'ok_epochs');

%% Modify data
% log all dropped trials
trial_info.orig_n_trials = length(trial_info.word_onset) + length(trial_info.ignore_trials);
trial_info.reject_trials = trial_info.trial_n(setdiff(1:length(trial_info.word_onset),ok_epochs));
trial_info.bad_trials    = vertcat(trial_info.ignore_trials', trial_info.reject_trials);
trial_info.ok_epochs     = ok_epochs;
trial_info.epoch_lim     = epoch_lim;
% cut out bad trials from all info fields
trial_info.condition_n   = trial_info.condition_n(ok_epochs);
trial_info.block_n       = trial_info.block_n(ok_epochs);
trial_info.trial_n       = trial_info.trial_n(ok_epochs);
trial_info.response_time = trial_info.response_time(ok_epochs);
trial_info.marker_time   = trial_info.marker_time(ok_epochs);
trial_info.onset_time    = trial_info.onset_time(ok_epochs);
trial_info.word_onset    = trial_info.word_onset(ok_epochs);
trial_info.resp_onset    = trial_info.resp_onset(ok_epochs);
trial_info.word          = trial_info.word(ok_epochs);
trial_info.color         = trial_info.color(ok_epochs);
trial_info.trialtype     = trial_info.trialtype(ok_epochs);
trial_info.blocktype     = trial_info.blocktype(ok_epochs);

%% Save Data
out_filename = strcat(epoch_dir,data_id,'_clean.mat');
save(out_filename,'data_ecog','header_ecog','trial_info');
