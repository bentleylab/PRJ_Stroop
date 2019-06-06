function SBJ10a_crRT_smANOVA(SBJ,an_id,stat_id)
%% function SBJ10a_crRT_smANOVA(SBJ,an_id,stat_id)
%   Run ANOVA with smart sliding windows for given time-frequency analysis
%   Sliding windows are smart and will not average pre- and post-RT data
%   Correlation with RT and regression of RT as confound supported
% INPUTS:
%   SBJ [str] - subject ID
%   an_id [str] - HFA analysis to run stats
%   stat_id [str] - ID of the statistical parameters and design
% OUTPUTS:
%   w2 [struct] - pseudo-FT structure with main ANOVA output
%   stat [FT struct] - output of correlation with RT if st.rt_corr==1
%   st [struct] - stat params loaded via stat_id

%% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

rng('shuffle'); % seed randi with time

%% Load Data
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);

load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
load(strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat'));

% Check if more than one frequency, error for now
if numel(hfa.freq)>1
    error('HFA has more than one frequency, can''t run on that for now!');
end

%% Check for RTs out of range
if st.min_rt
    % Find bad trials
    bad_rt_idx = trial_info.response_time<st.min_rt;
    st.bad_trial_n = trial_info.trial_n(bad_rt_idx);
    % Exclude those trials
    ti_fields = fieldnames(trial_info);
    orig_n_trials = numel(trial_info.trial_n);
    for f_ix = 1:numel(ti_fields)
        if numel(trial_info.(ti_fields{f_ix}))==orig_n_trials
            trial_info.(ti_fields{f_ix}) = trial_info.(ti_fields{f_ix})(~bad_rt_idx);
        end
    end
end

%% Build Design Matrix
design = cell([1 numel(st.groups)]);
levels = cell([1 numel(st.groups)]);
for grp_ix = 1:numel(st.groups)
    [levels{grp_ix}, ~, ~] = fn_condition_label_styles(st.groups{grp_ix});
    design{grp_ix} = nan([numel(trial_info.trial_n) 1]);
    for level_ix = 1:numel(levels{grp_ix})
        trl_idx = fn_condition_index(levels{grp_ix}{level_ix}, trial_info.condition_n, 'trial_info', trial_info);
        design{grp_ix}(trl_idx) = level_ix;
    end
end

%% Select data in stat window
cfg_trim = [];
% Remove excluded trials
if st.min_rt
    cfg_trim.trials = ~bad_rt_idx;
else
    cfg_trim.trials = 'all';
end

% Adjust stat_lim if needed
cfg_trim.latency = st.stat_lim;
for lim_ix = find(~cellfun(@isempty,st.lim_adj))
    if strcmp(st.lim_adj{lim_ix},'min(RT)')
        % S-locked to min(RT) or custom window
        cfg_trim.latency(lim_ix) = cfg_trim.latency(lim_ix)+min(trial_info.response_time);
    elseif lim_ix==2 && strcmp(st.lim_adj{lim_ix},'RT')
        % Align end to max(RT) to alter trim for custom windows
        cfg_trim.latency(lim_ix) = cfg_trim.latency(lim_ix)+max(trial_info.response_time);
    else
        % Nothing else makes sense...
        error(['What are you trying to do with st.lim_adj{' num2str(lim_ix) '} = ' st.lim_adj{lim_ix}]);
    end
end
hfa = ft_selectdata(cfg_trim,hfa);

%% Run correlations with RT
if st.rt_corr
    fprintf('================== Running Correlation with RT =======================\n');
    cfg_rt.design           = zscore(trial_info.response_time);
    cfg_rt.ivar             = 1;
    stat = ft_freqstatistics(cfg_rt, hfa);
end
% OUTPUT:
%   .rho   = correlation coefficients
%   .stat  = t values
%   .prob  = p values

%% Regress off Reaction Time
if st.regress_rt
    fprintf('================== Regressing RT =======================\n');
    cfg_conf = [];
    cfg_conf.output = 'residual';
    cfg_conf.confound = zscore(trial_info.response_time);
    hfa = ft_regressconfound(cfg_conf, hfa);
end
% OUTPUT:
%   As of August 2017, can no longer get residuals, model, and beta; now choose via .output argument
%   'beta'  = weights
%       .beta = 4D double
%   'model' = confounds * weights = X * X\Y
%       output is a .model struct that contains .powspctrm, .time, .trialinfo, etc.
%   'residual' = Yclean = Y - X * X\Y (i.e., the residuals after subtracting the model off
%       output is .powspctrm as normal
% NOTE: .stat and .prob output if .statistics = 'yes', but DO NOT TRUST THEM!!! (e.g., p values of 2...)

%% Average HFA in Sliding Windows
fprintf('================== Averaging HFA within Windows =======================\n');
% Sliding window parameters
if st.cust_win
    if ~strcmp(st.evnt_lab,'S') || ~strcmp(st.lim_adj{2},'RT')
        error('This is only meant for D analyses, check your options!');
    end
    % Define custom window per trial based on RT
    win_lim = zeros([numel(trial_info.trial_n) 2]);
    hfa_win = zeros([size(hfa.powspctrm,1) size(hfa.powspctrm,2)]);%size(hfa.powspctrm,3) size(win_lim,1)
    for trl_ix = 1:numel(trial_info.trial_n)
        pre_rt_ix = find(hfa.time<trial_info.response_time(trl_ix));
        win_lim(trl_ix,:) = [1 pre_rt_ix(end)];
        hfa_win(trl_ix,:) = squeeze(nanmean(hfa.powspctrm(trl_ix,:,1,win_lim(trl_ix,1):win_lim(trl_ix,2)),4));
    end
    win_center = round(mean(mean(win_lim,2)));
else
    win_lim    = fn_sliding_window_lim(squeeze(hfa.powspctrm(1,1,1,:)),...
        round(st.win_len*trial_info.sample_rate),...
        round(st.win_step*trial_info.sample_rate));
    
    % Average in windows
    hfa_win = zeros([size(hfa.powspctrm,1) size(hfa.powspctrm,2) size(win_lim,1)]);%size(hfa.powspctrm,3)
    for w_ix = 1:size(win_lim,1)
        hfa_win(:,:,w_ix) = squeeze(nanmean(hfa.powspctrm(:,:,1,win_lim(w_ix,1):win_lim(w_ix,2)),4));
    end
    win_center = round(mean(win_lim,2));
end

%% Run ANOVA
fprintf('================== Running ANOVA =======================\n');
% Create structure for w2 in fieldtrip style
w2.design    = design;
w2.cond      = st.groups;
w2.time      = hfa.time(win_center);
w2.win_lim   = win_lim;
w2.win_lim_s = hfa.time(win_lim);
w2.label     = hfa.label;
w2.dimord    = 'rpt_chan_time';
% w2.trial     = zeros([numel(w2.cond) length(hfa.label) length(w2.time)]);
w2.boot      = zeros([numel(w2.cond) length(hfa.label) length(w2.time) st.n_boots]);
% w2.pval      = w2.trial;


% Compute ANOVA and Explained Variance for real model
w2.trial = fn_mass_ANOVA(hfa_win,design);

% Compute ANOVA for permuted data
rand_design = design;
% b = '';
fprintf('boot #: ');
for boot_ix = 1:st.n_boots
%     m = sprintf(' permutation %d/%d', boot_ix, n_boots);
%     fprintf([b m]); b = repmat('\b',[1 length(m)]);
    fprintf('%i..',boot_ix);
    for grp_ix = 1:numel(st.groups)
        rand_design{grp_ix} = design{grp_ix}(randperm(size(design{grp_ix},1)));
    end
    w2.boot(:,:,:,boot_ix) = fn_mass_ANOVA(hfa_win,rand_design);
    if mod(boot_ix,20)==0
        fprintf('\n');
    end
end

% Compute statistics
w2.pval = sum(bsxfun(@ge,w2.boot,w2.trial),4)./st.n_boots; % sum(boots>real)/n_boots
w2.zscore   = norminv(1-w2.pval,0,1);
w2.bootmean = mean(w2.boot,4);
w2.bootstd  = std(w2.boot,[],4);
% w2 = rmfield(w2,'boot');
w2.zscore(isinf(w2.zscore)) = norminv(1-1/st.n_boots/2,0,1);

% Multiple Comparisons Correction within Channel
w2.qval = nan(size(w2.pval));
for ch_ix = 1:numel(w2.label)
    [~, ~, ~, w2.qval(:,ch_ix,:)] = fdr_bh(squeeze(w2.pval(:,ch_ix,:)));
end

%% Print results
% Prep report
sig_report_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id '_' an_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
if st.rt_corr; cond_lab = [st.groups {'RT'}]; else cond_lab = st.groups; end;
result_str = ['%-10s' repmat('%-10i',[1 numel(cond_lab)]) '\n'];

% Print header
fprintf(sig_report,'%s (n = %i)\n',SBJ,numel(w2.label));
fprintf(sig_report,[repmat('%-10s',[1 1+numel(cond_lab)]) '\n'],'label',cond_lab{:});

% Print summary lines (absolute)
if st.rt_corr
    fprintf(sig_report,result_str, 'count', sum(any(w2.qval(:,:,:)<st.alpha,3),2), sum(any(stat.mask(:,1,:),3)));
    fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
        [sum(any(w2.qval(:,:,:)<st.alpha,3),2)', sum(any(stat.mask(:,1,:),3))]./numel(w2.label));
else
    fprintf(sig_report,result_str, 'count', sum(any(w2.qval(:,:,:)<st.alpha,3),2));
    fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
        sum(any(w2.qval(:,:,:)<st.alpha,3),2)'./numel(w2.label));
end

% Print Channel Lines
sig_mat = zeros([numel(w2.label) numel(cond_lab)]);
for ch_ix = 1:numel(w2.label)
    % Consolidate to binary sig/non-sig
    for grp_ix = 1:numel(st.groups)
        if any(squeeze(w2.qval(grp_ix,ch_ix,:))<st.alpha)
            sig_mat(ch_ix,grp_ix) = 1;
        end
    end
    if st.rt_corr && any(stat.mask(ch_ix,1,:))
        sig_mat(ch_ix,grp_ix) = 1;
    end
    
    % Report on significant electrodes for this SBJ
    fprintf(sig_report,result_str,w2.label{ch_ix},sig_mat(ch_ix,:));
end

fclose(sig_report);

%% Save Results
out_fname = [SBJ_vars.dirs.proc SBJ '_smANOVA_ROI_' stat_id '_' an_id '.mat'];
if st.rt_corr
    save(out_fname,'-v7.3','w2','stat','st');
else
    save(out_fname,'-v7.3','w2','st');
end

end
