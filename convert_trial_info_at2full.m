%% Convert epochs to non-analysis_time format

clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));

SBJ      = 'IR21';%'IR21','IR31','IR32','IR35','IR39'};%'IR39'};%};%
data_id   = 'RC_WM';%{'RC_WM'},{'LAC_WM'},{'IH_CA'},{'LAC_WM'},{'RAC_WM'}};%,'RAC_BP_data'}};%{'RAC_WM'}};

epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
SBJ_dir   = fullfile('/home/knight/hoycw/PRJ_Stroop/data/',SBJ);
event_dir = fullfile(SBJ_dir,'03_events/');

%% Load Data
load(strcat(event_dir,SBJ,'_KLA_analysis_time.mat'));
load(strcat(event_dir,SBJ,'_',data_id,'_trial_info.mat'));
trial_info_KLA = trial_info;

%% Convert times
trial_info.word_onset = fn_convert_evnt_times_at2full(trial_info.word_onset,...
                                                            analysis_time,trial_info.sample_rate);
trial_info.resp_onset = fn_convert_evnt_times_at2full(trial_info.resp_onset,...
                                                            analysis_time,trial_info.sample_rate);

%% Save Output
out_filename = strcat(event_dir,SBJ,'_trial_info_full.mat');
save(out_filename,'trial_info');
