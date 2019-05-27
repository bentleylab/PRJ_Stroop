function SBJ07a_ERP_actv(SBJ,proc_id,an_id)
% Calculates ERPs and inter-trial phase coherence to identify which ERPs have
%   consistent phase locking and are likely true evoked responses, and should
%   therefore be used in future analyses
% clear all; %close all;
rng('shuffle'); % seed randi with time

%% Data Preparation
%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Preprocess the data
cfgpp = [];
cfgpp.demean    = 'yes';
cfgpp.hpfilter  = erp_vars.hp_yn;
cfgpp.hpfreq    = erp_vars.hp_freq;
if strcmp(erp_vars.hp_yn,'yes')
    cfgpp.hpfiltord = 4;            % Leaving blank causes instability error, 1 or 2 works 
end
cfgpp.lpfilter  = erp_vars.lp_yn;
cfgpp.lpfreq    = erp_vars.lp_freq;
roi_filt = ft_preprocessing(cfgpp,roi);

%% Cut into Trials
if strcmp(erp_vars.bsln_evnt,'stim')
    bsln_events = trial_info.word_onset;
    if strcmp(event_type,'stim')
        roi_filt_trl = fn_ft_cut_trials_equal_len(roi_filt,bsln_events,trial_info.condition_n',trial_lim_s*roi_filt.fsample);
    elseif strcmp(event_type,'resp')
        max_RT  = max(trial_info.response_time);
        % Check that baseline will be included in trial_lim_s
        if trial_lim_s(1)>erp_vars.bsln_lim(1)
            trial_cut = [erp_vars.bsln_lim(1) max_RT+trial_lim_s(2)];
        else
            trial_cut = trial_lim_s;
        end
        % Cut out to max_RT+trial_lim_s(2)+max(cfg_hfa.t_ftimwin)
        roi_filt_trl = fn_ft_cut_trials_equal_len(roi_filt,bsln_events,trial_info.condition_n',...
            round(trial_cut*roi.fsample));
    end
elseif strcmp(erp_vars.bsln_evnt,'resp')
    if strcmp(event_type,'resp')
        bsln_events = trial_info.resp_onset;
    else
        error('Why are you using resp as baseline when locking to stim???');
    end
    roi_filt_trl = fn_ft_cut_trials_equal_len(roi_filt,bsln_events,trial_info.condition_n',trial_lim_s*roi_filt.fsample);
end

%% Compute ERPs
% Baseline Correction
cfg_bsln = [];
cfg_bsln.demean         = erp_vars.demean_yn;
cfg_bsln.baselinewindow = erp_vars.bsln_lim;
roi_filt_trl_bsln = ft_preprocessing(cfg_bsln,roi_filt_trl);

% Average ERP
cfg_avg = [];
cfg_avg.keeptrials = 'yes';
roi_erp = ft_timelockanalysis(cfg_avg,roi_filt_trl_bsln);

% Trim to stat window
cfg_trim = [];
cfg_trim.latency = stat_lim;
erp_stat = ft_selectdata(cfg_trim,roi_erp);

%% Permute Trials for Null ERP
null_erp = NaN(n_boots, numel(roi_filt_trl.label), size(erp_stat.time,2));
for boot_ix = 1:n_boots
    % Circular shift
    shift_data = roi_filt_trl;
    shift_amt = randi(size(roi_filt_trl.time{1},2),numel(roi_filt_trl.trial),1);
    for t_ix = 1:numel(roi_filt_trl.trial)
        shift_data.trial{t_ix} = circshift(roi_filt_trl.trial{t_ix},shift_amt(t_ix),2);
    end
    
    % Baseline correct
    shift_data_bsln = ft_preprocessing(cfg_bsln,shift_data);
    
    % Average ERP
    cfg_avg.keeptrials = 'no';
    shift_erp = ft_timelockanalysis(cfg_avg,shift_data_bsln);
    
    % Trim and stash
    shift_erp = ft_selectdata(cfg_trim,shift_erp);
    null_erp(boot_ix,:,:) = shift_erp.avg;
end

%% Compute Statistics
actv_ch = {};
actv_ch_epochs = {};
actv_win_dp = actv_win*roi_filt_trl.fsample;
for ch_ix = 1:numel(erp_stat.label)
    % Compute t-test per time point
    n_tests   = size(erp_stat.time,2);
    pvals     = NaN([1 n_tests]);
    for time_ix = 1:n_tests
        null_amp_dist = sort(abs(null_erp(:,ch_ix,time_ix)));                       % Sort absolute values
        amp_logical = logical(abs(mean(erp_stat.trial(:,ch_ix,time_ix))) > null_amp_dist);    % Find n_boots less than real value
        pvals(time_ix) = (numel(null_amp_dist)-numel(find(amp_logical)))/n_boots;  % convert to proportion of n_boots
    end

    % Find epochs with significant task activations
%     [~, qvals] = mafdr(pvals); % Errors on some random ch during SBJ08a_HFA_actv (e.g., ROF8-9 in IR32), so try this:
%     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [~, ~, ~, qvals] = fdr_bh(pvals);
    actv_mask = qvals<=0.05;
    actv_chunks = fn_find_chunks(actv_mask);
    actv_chunks(actv_mask(actv_chunks(:,1))==0,:) = [];
    actv_chunk_sz = diff(actv_chunks,1,2)+1;
    actv_epochs = actv_chunks(actv_chunk_sz > actv_win_dp,:);
    if ~isempty(actv_epochs)
        actv_ch = {actv_ch{:} erp_stat.label{ch_ix}};
        actv_ch_epochs = {actv_ch_epochs{:}, erp_stat.time(actv_epochs)};
    end
    
%     % Plot ERP vs. null distribution
%     plot(squeeze(null_erp(:,ch_ix,:))','Color',[0.7 0.7 0.7]); hold on;
%     shadedErrorBar([],squeeze(mean(erp_stat.trial(:,ch_ix,:),1)),...
%         squeeze(std(erp_stat.trial(:,ch_ix,:),1))./sqrt(size(erp_stat.trial,1)),...
%         {'LineWidth',3,'Color','k'},1);
% %     plot(squeeze(mean(erp_stat.trial(:,ch_ix,:),1)),'LineWidth',3,'Color','k');
%     plot(pvals,'LineWidth',2,'Color','r');
%     pause;
%     hold off;
end


%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_ROI_',an_id,'.mat');
fprintf('Saving %s\n',data_out_filename);
save(data_out_filename,'-v7.3','roi_erp','actv_ch','actv_ch_epochs');

end
