function SBJ07a_ERP_ITPC(SBJ,conditions,pipeline_id,an_id)
% Calculates ERPs and inter-trial phase coherence to identify which ERPs have
%   consistent phase locking and are likely true evoked responses, and should
%   therefore be used in future analyses
% clear all; %close all;

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
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

% Select Conditions of Interest
[cond_lab, ~, ~] = fn_condition_label_styles(conditions);
cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n));
end

% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);

%% Compute ERPs
% Preprocess the data
cfgpp = [];
cfgpp.hpfilter  = erp_vars.hp_yn;
cfgpp.hpfreq    = erp_vars.hp_freq;
if strcmp(erp_vars.hp_yn,'yes')
    cfgpp.hpfiltord = 4;            % Leaving blank causes instability error, 1 or 2 works 
end
cfgpp.lpfilter  = erp_vars.lp_yn;
cfgpp.lpfreq    = erp_vars.lp_freq;
roi_filt = ft_preprocessing(cfgpp,roi);

% Cut into Trials
if ~strcmp(erp_vars.bsln_evnt,event_type)
    error(['an_vars ' an_id ' is trying to baseline ERPs to a different event than it is locked to!']);
end
if strcmp(event_type,'stim')
    events = trial_info.word_onset;
elseif strcmp(event_type,'resp')
    events = trial_info.resp_onset;
else
    error(stract('ERROR: unknown event_type ',event_type));
end
roi_filt_trl = fn_ft_cut_trials_equal_len(roi_filt,events,trial_info.condition_n',trial_lim_s*roi_filt.fsample);

% Baseline Correction
cfg = [];
cfg.demean         = erp_vars.demean_yn;
cfg.baselinewindow = erp_vars.bsln_lim;
roi_filt_trl = ft_preprocessing(cfg,roi_filt_trl);

% Average ERP
cfgavg = [];
cfgavg.keeptrials = 'yes';
roi_erp = ft_timelockanalysis(cfgavg,roi_filt_trl);

%% Compute ITPC
% cut into trials
roi_trl = fn_ft_cut_trials_equal_len(roi,events,trial_info.condition_n',trial_lim_s*roi.fsample);
% frequency decomposition
roi_freq = ft_freqanalysis(cfg_tfr,roi_trl);

% ITPC computation
itc = [];
itc.label = roi_freq.label;
itc.freq  = roi_freq.freq;
itc.time  = roi_freq.time;
itc.dimord = 'chan_freq_time';

F = roi_freq.fourierspctrm;   % copy the fourier spectrum
N = size(F,1);                  % number of trials

itc.itpc = F./abs(F);           % divide by amplitude
itc.itpc = sum(itc.itpc,1);     % sum angles
itc.itpc = abs(itc.itpc)/N;     % take absolute value and normalize
itc.itpc = squeeze(itc.itpc);   % remove first singleton dimension

% %% Plot
% figure
% imagesc(itc.time, itc.freq, squeeze(itc.itpc(1,:,:)));
% axis xy
% title('ITPC');

%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.SBJ,'04_proc/',SBJ,'_',conditions,'_ROI_',an_id,'.mat');
fprintf('Saving %s\n',data_out_filename);
save(data_out_filename,'-v7.3','roi_erp','stat');

end
