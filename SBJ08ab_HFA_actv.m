function SBJ08ab_HFA_actv(SBJ,an_id,stat_id)
%% function SBJ08ab_HFA_actv(SBJ,an_id,stat_id)
%   Calculates activation relative to baseline (0) via point wise t-test
%   t-test pval converted to qval via FDR
%   Smart windows based on SBJ-specific RTs
%   Writes text report of significance
% INPUTS:
%   SBJ [str] - subject ID
%   an_id [str] - HFA analysis to run stats
%   stat_id [str] - ID of the statistical parameters
% OUTPUTS:
%   actv [struct] - pseudo-FT structure with main outputs
%   st [struct] - stat params loaded via stat_id

if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end

%% Set up paths
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Data Preparation
SBJ_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m'];
eval(SBJ_vars_cmd);
an_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id '_vars.m'];
eval(an_vars_cmd);
stat_vars_cmd = ['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m'];
eval(stat_vars_cmd);

% Load Data
hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'.mat');
load(hfa_fname);
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));

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

%% Select data in stat window and remove excluded trials
cfg_trim = [];
if st.min_rt
    cfg_trim.trials = ~bad_rt_idx;
else
    cfg_trim.trials = 'all';
end
cfg_trim.latency = st.stat_lim;
% Adjust stat_lim based on events
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
hfa_stat = ft_selectdata(cfg_trim,hfa);

%% Average HFA for D window
if st.cust_win
    fprintf('================== Averaging HFA within D Window =======================\n');
    if ~strcmp(st.evnt_lab,'S') || ~strcmp(st.lim_adj{2},'RT')
        error('This is only meant for D analyses, check your options!');
    end
    % Create template hfa struct
    cfg = []; cfg.latency = hfa_stat.time(1);
    hfa_tmp = ft_selectdata(cfg,hfa_stat);
    
    % Define custom window per trial based on RT
    win_lim = zeros([numel(trial_info.trial_n) 2]);
    for trl_ix = 1:numel(trial_info.trial_n)
        pre_rt_ix = find(hfa_stat.time<trial_info.response_time(trl_ix));
        win_lim(trl_ix,:) = [1 pre_rt_ix(end)];
        hfa_tmp.powspctrm(trl_ix,:) = ...
            squeeze(nanmean(hfa_stat.powspctrm(trl_ix,:,1,win_lim(trl_ix,1):win_lim(trl_ix,2)),4));
    end
    hfa_tmp.time = hfa_stat.time(round(mean(mean(win_lim,2))));
    hfa_stat = hfa_tmp;
    clear hfa_tmp
end

%% Run Statistics
fprintf('===================================================\n');
fprintf('--------------------- Statistics ------------------\n');
fprintf('===================================================\n');

% Create structure for actv in fieldtrip style
actv.cond    = st.groups;
actv.time    = hfa_stat.time;
actv.label   = hfa_stat.label;
actv.dimord  = 'chan_time';
actv.avg     = squeeze(nanmean(hfa_stat.powspctrm(:,:,1,:),1));
if st.cust_win; actv.avg = actv.avg'; end % ensure dimord
actv.pval    = zeros(size(actv.avg));
actv.qval    = zeros(size(actv.avg));
actv.mask    = zeros(size(actv.avg));
if st.cust_win
    actv.win_lim = win_lim;
end

% Compute statistics
for ch_ix = 1:numel(hfa_stat.label)
    % Compute t-test per time point
    for time_ix = 1:numel(hfa_stat.time)
        [~, actv.pval(ch_ix,time_ix)] = ttest(squeeze(hfa_stat.powspctrm(:,ch_ix,1,time_ix)));
    end
    
    % Find epochs with significant task activations
%     [~, qvals] = mafdr(pvals); % Errors on some random ch (e.g., ROF8-9
%     in IR32), so I'm trying the below function
%     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [~, ~, ~, actv.qval(ch_ix,:)] = fdr_bh(actv.pval(ch_ix,:));
    actv.mask(ch_ix,:) = actv.qval(ch_ix,:)<=st.alpha;
    
    % Remove epochs less than actv_win
    actv_chunks = fn_find_chunks(actv.mask(ch_ix,:));
    actv_chunks(actv.mask(ch_ix,actv_chunks(:,1))==0,:) = [];
    actv_chunk_sz = diff(actv_chunks,1,2)+1;
    bad_chunks = actv_chunks(actv_chunk_sz < trial_info.sample_rate*st.actv_win,:);
    for chunk_ix = 1:size(bad_chunks,1)
        actv.mask(ch_ix,bad_chunks(chunk_ix,1):bad_chunks(chunk_ix,2)) = 0;
    end
end

%% Print results
% Compile positive and negative stats
cond_lab = {'actv' 'pos' 'neg'};
sig_mat = zeros([numel(actv.label) numel(cond_lab)]);
for ch_ix = 1:numel(actv.label)
    % Consolidate to binary sig/non-sig
    if any(squeeze(actv.qval(ch_ix,:))<st.alpha)
        sig_mat(ch_ix,1) = 1;
        % Flag whether positive or negative
        sig_idx = squeeze(actv.qval(ch_ix,:))<st.alpha;
        if any(squeeze(actv.avg(ch_ix,sig_idx))>0)
            sig_mat(ch_ix,2) = 1;
        end
        if any(squeeze(actv.avg(ch_ix,sig_idx))<0)
            sig_mat(ch_ix,3) = 1;
        end
    end
end

% Prep report
sig_report_fname = [hfa_fname(1:end-4) '_' stat_id '_sig_report.txt'];
if exist(sig_report_fname)
    system(['mv ' sig_report_fname ' ' sig_report_fname(1:end-4) '_bck.txt']);
end
sig_report = fopen(sig_report_fname,'a');
result_str = ['%-10s' repmat('%-10i',[1 numel(cond_lab)]) '\n'];

% Print header
fprintf(sig_report,'%s (n = %i)\n',SBJ,numel(actv.label));
fprintf(sig_report,[repmat('%-10s',[1 1+numel(cond_lab)]) '\n'],'label',cond_lab{:});

% Print summary lines (absolute)
fprintf(sig_report,result_str, 'count', sum(sig_mat,1));
fprintf(sig_report,strrep(result_str,'i','.3f'), 'percent',...
    sum(sig_mat,1)./numel(actv.label));

% Print Channel Lines
for ch_ix = 1:numel(actv.label)
    % Report on significant electrodes for this SBJ
    fprintf(sig_report,result_str,actv.label{ch_ix},sig_mat(ch_ix,:));
end

fclose(sig_report);

%% Save Results
out_fname = strcat(hfa_fname(1:end-4),'_',stat_id,'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',out_fname);
fprintf('===================================================\n');
save(out_fname,'-v7.3','actv','st');

end
