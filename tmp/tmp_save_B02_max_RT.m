%% Stroop analysis - High Gamma Power on Single Trials
clear all; %close all
addpath(genpath('/home/knight/hoycw/PRJ_Stroop/scripts/'));
addpath(genpath('/home/knight/hoycw/Apps/export_fig-master/'));
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults
% addpath(genpath('/home/knight/hoycw/Apps/EEGlab/eeglab12_0_2_6b/functions/sigprocfunc/'));

% Analysis Parameters
SBJs       = {'IR21','IR31','IR35','IR39'};%'IR32'};%};%
ch_ids     = {{'RC1-RC2'},{'LAC1-LAC2'},{'LAC1-LAC2'},{'RAC1-RAC2'}};%};%,{'RC_WM'},{'LAC_WM'},{'IH_CA'},{'LAC_WM'},{'RAC_WM'}};%,'RAC_BP_data'}};%{'RAC_WM'}};
cond_name  = 'CI';              % conditions to average
analyses   = {'ERP','theta','HG'};
%   CI   = only con and inc trial types for all blocks
%   pcon = proportion congruency (con and inc for both mcon and minc)
%   !!!conseq = congruency sequence effects (cC, cI, iC, iI)
epoch_lim = [-200 2000];           % ID of data used for trial rejection (these epochs are NOT used)
buff_lim  = [1000, 1000];            % buffers are for cutting time series and then plotting
sig_lim   = [1000, 1000];

% Analysis parameters
HG_type      = 'wideband';                     % only for HG: 'wideband', 'multiband'
event        = 'resp';                          % 'stim'/'resp': event to lock trials
bsln_type    = {'demean','zscore','zscore'};    % 'zscore', 'demean', 'none'
bsln_event   = 's';                    % 's'/'r': event to lock baselining
bsln_lim     = [250, -50];             % ms before and after event
smooth_it    = 1;                      % smooth the result before averaging (0/1)
smooth_freq  = 10;
trial_type   = 'datapad';               % 'nanpad', 'datapad'
sig_win_len  = 50;                      % Length of significance testing window in ms
sig_win_step = 25;                      % Step size to slide teh sig window forward
sig_clust_len= 100;                     % length in ms over which consecutive windows must be significant
sig_nboots   = 1000;                    % number of iterations for permutation testing

% Plotting parameters
save_fig      = 1;
vis_fig       = 'off';                % 'on' or 'off' to determine if plot is shown
sem_alpha     = 0.5;                  % transparency of sem shading (0:1)
clim_perc     = [5 95];              % colormap percentile limits
plot_lim      = [500 500];            % ms to plot before the event and **after RT**
x_step        = 250;                  % step of x tick marks
event_ln_width= 2;
fig_type      = 'png';

%% Process parameters
% Condition Parameters
[cond_lab, cond_colors, cond_style] = fn_condition_label_styles(cond_name);

% Analysis Parameters
all_analysis_id = '';
all_bsln_type = '';
for an_ix = 1:length(analyses)
    [analysis_id{an_ix}, env_it{an_ix}, filt_lim{an_ix}] =...
        fn_B02_analysis_params(analyses{an_ix},HG_type);
    all_analysis_id = strcat(all_analysis_id, analysis_id{an_ix}(1));
    all_bsln_type = strcat(all_bsln_type, bsln_type{an_ix}(1));
end

% Filtering
if smooth_it == 0
    smooth_id = '';
elseif smooth_it == 1
    smooth_id = strcat('_sm',num2str(smooth_freq));
else
    error('smooth_it not in [0,1]');
end

% Epoching/Averaging
epoch_id = [num2str(floor(epoch_lim(1)),'%+d') '.' num2str(floor(epoch_lim(2)),'%+d')];
bsln_id  = ['_' all_bsln_type '.' bsln_event num2str(bsln_lim(1)) '.' num2str(bsln_lim(2))];
if strcmp(trial_type,'nanpad')
    trl_id = 'nan';
elseif strcmp(trial_type,'datapad')
    trl_id = 'data';
else
    error(strcat('Unknown trial averaging type: ',trial_type));
end
if length(bsln_type)~=length(analyses)
    error('Mismatched number of bsln_type and analyses');
end

% Plotting (Y axis)
if (plot_lim(1)>buff_lim(1)) || (plot_lim(2)>buff_lim(2))
    error('plot_lim exceed buff_lim')
end
plot_lim_id = strcat('_xlim',num2str(plot_lim(1)),'.',num2str(plot_lim(2)));
y_scale_id = strcat('_c',num2str(clim_perc(1)),'.',num2str(clim_perc(2)));

%% Analysis Loop
for SBJ_ix = 1:length(SBJs)
    %     for data_id_ix = 1:length(data_ids{SBJ_ix})
    SBJ = SBJs{SBJ_ix};
    data_id = strcat(SBJ,'_ft_preproc');%strcat(SBJ,'_',data_ids{SBJ_ix}{data_id_ix});
    SBJ_dir = strcat('/home/knight/hoycw/PRJ_Stroop/data/',SBJ,'/');
    fig_dir  = strcat('/home/knight/hoycw/PRJ_Stroop/results/multifreq/',SBJ,'/',cond_name,'/',event,'/');
    if ~exist(fig_dir,'dir')
        mkdir(fig_dir);
    end
    
    %% Load data
    load(strcat(SBJ_dir,'02_preproc/',data_id,'.mat'));
    cfg = [];
    cfg.channel = ch_ids{SBJ_ix};
    data = ft_selectdata(cfg,data);
    [data_ecog, header_ecog] = fn_format_data_ft2KLA(data);
    %         load(strcat(SBJ_dir,'04_proc/',data_id,'.mat'),'data_ecog','header_ecog');
    load(strcat(SBJ_dir,'03_events/',SBJ,'_trial_info_full.mat'));
    
    % need only ok_epochs index, the rest is all in 04_proc file (variables seem same as of 10/10/16...)
    bob_bad_epochs = load(strcat(SBJ_dir,'03_events/',SBJ,'_bob_bad_epochs.mat'),'bad_epochs');
    bob_bad_epochs = bob_bad_epochs.bad_epochs;
    %         if strcmp(data_id,'IR39_RAC_BP_data') || strcmp(data_id,'IR39_RAC_ft_KLA')  %!!! SHITTY WORK AROUND!!!
    %             KLA_ok_epochs = load(strcat(SBJ_dir,'06_epochs/IR39_RAC_WM_',epoch_id,'.mat'),'ok_epochs');
    %         else
    %             KLA_ok_epochs = load(strcat(SBJ_dir,'06_epochs/',data_id,'_',epoch_id,'.mat'),'ok_epochs');
    %         end
    %         KLA_ok_epochs = KLA_ok_epochs.ok_epochs;
    
    % Correct all info to only good trials
    if header_ecog.sample_rate~=1000
        error('\n!!!WARNING: Sampling rate is not 1 kHz, timing will be off!!\n')
    end
    % Toss epochs that overlap with bad_epochs from Bob
    bad_samples = [];
    time_epoch_idx = fn_epoch_cuts_datapad(1:size(data_ecog,2),trial_info.word_onset,...
        trial_info.resp_onset,buff_lim);
    for bad_epoch_ix = 1:size(bob_bad_epochs,1)
        bad_samples = [bad_samples bob_bad_epochs(bad_epoch_ix,1):bob_bad_epochs(bad_epoch_ix,2)];
    end
    bob_ok_epochs = [];
    for epoch_ix = 1:size(time_epoch_idx,1)
        if isempty(intersect(time_epoch_idx(epoch_ix,:),bad_samples))
            bob_ok_epochs = [bob_ok_epochs epoch_ix];
        end
    end
    ok_epochs = bob_ok_epochs;%intersect(bob_ok_epochs,KLA_ok_epochs);
    fprintf('Num Bob ok epochs: %i\n',length(bob_ok_epochs));
    fprintf('Num KLA ok epochs: %i\n',length(ok_epochs));
    fprintf('Overlap          : %i\n',length(ok_epochs));
    
    % Toss trials with no RT
    no_RT_epochs = find(isnan(trial_info.resp_onset)==1);
    good_epochs = setdiff(ok_epochs,no_RT_epochs);
    fprintf('Num trials w/o RT: %i\n',sum(no_RT_epochs));
    fprintf('FINAL NUM GOOD TRIALS: %i\n',length(good_epochs));
    % con = 1-3, neu = 4-6, inc = 7-9
    % within those: same, mic, mcon
    % trial_type = NaN(size(trial_info.condition_n));
    word_onset = trial_info.word_onset(good_epochs);
    resp_onset = trial_info.resp_onset(good_epochs);
    RTs        = round(1000*trial_info.response_time(good_epochs)); % converts sec to ms
    [RTs_sorted,RTs_sort_idx] = sort(RTs);
    max_RT     = max(RTs);
    save(strcat(SBJ_dir,'03_events/tmp_',SBJ,'_B02_max_RT.mat'),'max_RT');
    %     end
end
