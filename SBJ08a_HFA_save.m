function SBJ08a_HFA_save(SBJ,pipeline_id,an_id)
% Calculates high frequency activity, computes cluster-based statistics, and plots the results
% clear all; %close all;

% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
% eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);

% Load Data
load(strcat(SBJ_vars.dirs.preproc,SBJ,'_preproc_',pipeline_id,'.mat'));
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

%% Select Channel(s)
cfgs = [];
cfgs.channel = SBJ_vars.ch_lab.ROI;
roi = ft_selectdata(cfgs,data);
roi_fsample = roi.fsample;
clear data;

%% Cut into Trials
% Padding:
%   At a minimum, trial_lim_s must extend 1/2*max(filter_window) prior to 
%   the first data point to be estimated to avoid edges in the filtering window.
%   (ft_freqanalysis will return NaN for partially empty windows, e.g. an edge pre-trial,
%   but ft_preprocessing would return a filtered time series with an edge artifact.)
%   Also note this padding buffer should be at least 3x the slowest cycle
%   of interest.
if strcmp(HFA_type,'multiband')
    pad_len = 0.5*max(cfg_hfa.t_ftimwin)*3;
elseif any(strcmp(HFA_type,{'broadband','hilbert'}))
    % add 250 ms as a rule of thumb, or onger if necessary
    pad_len = 0.5*max([1/min(fois)*3 0.25]);
end
% Add extra 10 ms just because trimming back down to trial_lim_s exactly leave
%   one NaN on the end (smoothing that will NaN out everything)
trial_lim_s_pad = [trial_lim_s(1)-pad_len trial_lim_s(2)+pad_len+0.01];
%!!! BEWARE: This is suboptimal for R-locked because trial_lim_s(1)=-0.5,
%which adds 250ms to the time series that isn't necessary; the realign_tfr
%function should still cut to the desired data, but it'll take longer.

% Always normalize to pre-stimulus baseline for HFA
bsln_events = trial_info.word_onset;
if strcmp(event_type,'stim')
    % Cut to desired trial_lim_s
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,trial_info.condition_n',...
        round(trial_lim_s_pad*roi_fsample));
elseif strcmp(event_type,'resp')
    % Check that baseline will be included in trial_lim_s
    if trial_lim_s(1)>bsln_lim(1)
        error(['ERROR: trial_lim_s does not include bsln_lim for an_id = ' an_id]);
    end
    % Cut out to max_RT+trial_lim_s(2)+max(cfg_hfa.t_ftimwin)
    max_RT  = max(trial_info.response_time);
    roi_trl = fn_ft_cut_trials_equal_len(roi,bsln_events,trial_info.condition_n',...
        round([trial_lim_s_pad(1) max_RT+trial_lim_s_pad(2)]*roi_fsample));
else
    error(['Unknown event_type: ' event_type]);
end
clear roi;

%% Compute HFA
fprintf('===================================================\n');
fprintf('------------------ HFA Calculations ---------------\n');
fprintf('===================================================\n');
if strcmp(HFA_type,'multiband')
    cfg_hfa.trials = 'all';
    hfa = ft_freqanalysis(cfg_hfa, roi_trl);
elseif strcmp(HFA_type,'hilbert')
    % Hilbert method to extract power
    hfas = cell(size(fois));
    orig_lab = roi_trl.label;
    for f_ix = 1:numel(fois)
        cfg_hfa.bpfreq = bp_lim(f_ix,:);
        fprintf('\n------> %s filtering: %.03f - %.03f\n', HFA_type, bp_lim(f_ix,1), bp_lim(f_ix,2));
        hfas{f_ix} = ft_preprocessing(cfg_hfa,roi_trl);
        hfas{f_ix}.label = strcat(hfas{f_ix}.label,[':' num2str(fois(f_ix),4)]);
    end
    % Treat different freqs as channels
    hfa_tmp = hfas{1};      % Save to plug in averaged data
    hfa = ft_appenddata([], hfas{:}); clear hfas;
elseif strcmp(HFA_type,'broadband')
    error('Stop using broadband and use filter-hilbert or multitapers you dummy!');
    %         % Filter to single HFA band
    %         cfgpp=[];.hpfilter='yes';.hpfreq=70;.lpfilter='yes';lpfreq=150;
    %         roi = ft_preprocessing(cfgpp,roi);
else
    error('Unknown HFA_type provided');
end

% Trim back down to original trial_lim_s to exclude NaNs
cfg_trim = [];
cfg_trim.latency = trial_lim_s;
if strcmp(event_type,'resp')
    cfg_trim.latency(2) = max_RT+trial_lim_s(2);
end
hfa = ft_selectdata(cfg_trim,hfa);
hfa_tmp = ft_selectdata(cfg_trim,hfa_tmp);

%% Baseline Correction
fprintf('===================================================\n');
fprintf('---------------- Baseline Correction --------------\n');
fprintf('===================================================\n');
switch bsln_type
    case {'zboot', 'zscore'}
        if strcmp(HFA_type,'multiband')
            hfa = fn_bsln_ft_tfr(hfa,bsln_lim,bsln_type,n_boots);
        elseif any(strcmp(HFA_type,{'broadband','hilbert'}))
            hfa = fn_bsln_ft_filtered(hfa,bsln_lim,bsln_type,n_boots);
        else
            error('Unknown HFA_type provided');
        end
    case {'relchange', 'demean', 'my_relchange'}
        error(['bsln_type ' bsln_type ' is not compatible with one-sample t test bsln activation stats']);
%         cfgbsln = [];
%         cfgbsln.baseline     = bsln_lim;
%         cfgbsln.baselinetype = bsln_type;
%         cfgbsln.parameter    = 'powspctrm';
%         hfa = ft_freqbaseline(cfgbsln,hfa);
    otherwise
        error(['No baseline implemented for bsln_type: ' bsln_type]);
end

%% Smooth Power Time Series
if smooth_pow_ts
    % error catches
    if ~strcmp(lp_yn,'yes')
        if strcmp(hp_yn,'yes')
            error('Why are you only high passing?');
        else
            error('Why is smooth_pow_ts yes but no lp or hp?');
        end
    end
    fprintf('===================================================\n');
    fprintf('----------------- Filtering Power -----------------\n');
    fprintf('===================================================\n');
    if isfield(hfa,'powspctrm')
        for ch_ix = 1:numel(hfa.label)
            for f_ix = 1:numel(hfa.freq)
                if strcmp(lp_yn,'yes') && strcmp(hp_yn,'yes')
                    hfa.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_bandpass(...
                        hfa.powspctrm(:,ch_ix,f_ix,:), roi_fsample, hp_freq, lp_freq);
                elseif strcmp(lp_yn,'yes')
                    hfa.powspctrm(:,ch_ix,f_ix,:) = fn_EEGlab_lowpass(...
                        hfa.powspctrm(:,ch_ix,f_ix,:), roi_fsample, lp_freq);
                else
                    error('weird non-Y/N filtering options!');
                end
            end
        end
    else
        if strcmp(lp_yn,'yes') && strcmp(hp_yn,'yes')
            hfa.trial{1} = fn_EEGlab_bandpass(hfa.trial{1}, roi_fsample, hp_freq, lp_freq);
        elseif strcmp(lp_yn,'yes')
            hfa.trial{1} = fn_EEGlab_lowpass(hfa.trial{1}, roi_fsample, lp_freq);
        else
            error('weird non-Y/N filtering options!');
        end
    end
end

%% Merge multiple bands
if strcmp(HFA_type,'multiband')
    cfg_avg = [];
    cfg_avg.freq = 'all';
    cfg_avg.avgoverfreq = 'yes';
    hfa = ft_selectdata(cfg_avg,hfa);
elseif strcmp(HFA_type,'hilbert')
    hfa_tmp.label = orig_lab;
    lab_ix = 1;
    ch_check = zeros(size(hfa.label));
    for ch_ix = 1:numel(hfa_tmp.label)
        ch_lab_ix = find(~cellfun(@isempty,strfind(hfa.label,hfa_tmp.label{ch_ix})));
        ch_check(ch_lab_ix) = ch_check(ch_lab_ix)+ch_ix;
        for t_ix = 1:numel(hfa_tmp.trial)
            hfa_tmp.trial{t_ix}(ch_ix,:) = mean(hfa.trial{t_ix}(ch_lab_ix,:),1);
        end
        lab_ix = ch_lab_ix(end)+1;
    end
    hfa = hfa_tmp; clear hfa_tmp;
end

%% Re-align to event of interest if necessary (e.g., response)
if strcmp(event_type,'resp')
    if strcmp(HFA_type,'multiband')
        hfa = fn_realign_tfr_s2r(hfa,trial_info.response_time,trial_lim_s);
    elseif any(strcmp(HFA_type,{'broadband','hilbert'}))
        hfa = fn_realign_filt_s2r(hfa,trial_info.response_time,trial_lim_s);
    end
elseif ~strcmp(event_type,'stim')
    error(['ERROR: unknown event_type ' event_type]);
end

%% Downsample
if resample_ts && hfa.fsample~=resample_freq
    cfgrs = [];
    cfgrs.resamplefs = resample_freq;
    cfgrs.detrend = 'no';
    hfa = ft_resampledata(cfgrs, hfa);
end

%% Save Results
data_out_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',data_out_filename);
fprintf('===================================================\n');
save(data_out_filename,'-v7.3','hfa');

end
