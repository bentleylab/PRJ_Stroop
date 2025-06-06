function SBJ13a_FOOOF_PSD_save(SBJ,an_id)
% Calculates PSDs for given dataset, saves frequencies and PSDs for FOOOF
% in python

% clear all; %close all;
%% Check which root directory
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set Up Directories
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
% Load Data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);

if strcmp(proc_id,'main_ft')
    load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',proc_id,'.mat'));
    % Load elec
    if strcmp(atlas_id,'Yeo7') || strcmp(atlas_id,'Yeo17')
        elec_space = 'mni_v';
    else
        elec_space = 'pat';
    end
    elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' proc_id '_' elec_space '_' atlas_id '.mat'];
    load(elec_fname);
    % Select ROI-specific elecs
    [roi_list, ~] = fn_roi_label_styles(roi_id);
    % Get ROI info per elec
    elec.roi = fn_atlas2roi_labels(elec.atlas_label,atlas_id,roi_id);
    atlas_in_elecs = {};
    for roi_ix = 1:numel(roi_list)
        atlas_in_elecs = [atlas_in_elecs; elec.label(strcmp(elec.roi,roi_list{roi_ix}))];
    end
    cfgs = []; cfgs.channel = atlas_in_elecs;
    elec = fn_select_elec(cfgs, elec);
    data = ft_selectdata(cfgs,data);
elseif strcmp(proc_id,'srate')
    if any(SBJ_vars.low_srate)
        load(strcat(SBJ_vars.dirs.import,SBJ,'_',num2str(SBJ_vars.low_srate(1)),'hz.mat'));
    else
        load(strcat(SBJ_vars.dirs.import,SBJ,'_1000hz.mat'));
    end
    elec = fn_load_elec_orig(SBJ,'pat','');
end
if ~all(strcmp(data.label,elec.label))
    error('Mismatch between data and elec labels!');
end
% load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
% 
% % Select Conditions of Interest
% [cond_lab, ~, ~] = fn_condition_label_styles(conditions);
% cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
% for cond_ix = 1:length(cond_lab)
%     % Get binary condition index
%     cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
%         trial_info.condition_n));
% end

%% Compute PSDs
[tmp_psd,freqs] = pwelch(data.trial{1}(1,:),wn_len,[],[],data.fsample);   % as in fooof_mat docs
psds = zeros(numel(elec.label),size(tmp_psd,1));
psds(1,:) = tmp_psd;
for ch_ix = 2:numel(data.label)
    % pwelch(signal, window, noverlap, f, sample_freq)
    %   computed per column, returned in corresponding columns
    %   window (int) - length of data segment to comptue PSDs
    %   noverlap (int) - # data points overlapping between windows, default is 50% of window
%    [fft_data,freqs] = pwelch(data.trial{1}(ch_ix,:),2048,0,2048,sample_freq);
    [psds(ch_ix,:),~] = pwelch(data.trial{1}(ch_ix,:),wn_len,[],[],data.fsample);
end

%% Save out in python readable format
labels = data.label;
rois = elec.roi;
out_fname = [SBJ_vars.dirs.preproc,SBJ,'_',an_id,'_PSDs.mat'];
save(out_fname,'freqs','psds','labels','rois','roi_id','atlas_id','foi','pk_bw_lim','mn_pk_amp','max_n_pks');

end
