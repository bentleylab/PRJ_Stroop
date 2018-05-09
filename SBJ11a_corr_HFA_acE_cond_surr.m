function SBJ11a_corr_HFA_acE_cond_surr(SBJ,pipeline_id,stat_id,an_id_main,plt_id)
% Build connectivity matrix based on HFA correlations
%   non-parametric stats via circular shift of trial time series

%% Data Preparation
% Set up paths
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/');
addpath('/home/knight/hoycw/PRJ_Stroop/scripts/utils/');
addpath('/home/knight/hoycw/Apps/fieldtrip/');
ft_defaults

%% Load Results
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/an_vars/' an_id_main '_vars.m']);
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

% Load data
hfa_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id_main,'.mat');
load(hfa_filename1,'hfa');

% Load condition info
load(strcat(SBJ_vars.dirs.events,SBJ,'_trial_info_final.mat'));
% Select Conditions of Interest
if ~strcmp(stat_vars.conditions,'CI')
    error('fix the corr_mat difference!');
end
[cond_lab, ~, ~] = fn_condition_label_styles(stat_vars.conditions);
cond_idx = false([length(cond_lab) length(trial_info.trial_n)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n));
end

% Load ROI and GM/WM info
einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
load(einfo_filename);
% Electrode Info Table:
%   label- name of electrode
%   ROI- specific region
%   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%   ROI2- specific region of second electrode
%   tissue- primary tissue type
%   GM weight- percentage of electrode pair in GM
%   Out- 0/1 flag for whether this is partially out of the brain

%% Prep Data
% Channel Selection
if stat_vars.only_actv_ch
    actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',stat_vars.an_id_s,'_mn',num2str(stat_vars.actv_win),'.mat');
    tmp = load(actv_filename,'actv_ch','actv_ch_epochs');
    actv_ch{1} = tmp.actv_ch; actv_ch_epochs{1} = tmp.actv_ch_epochs;
    actv_filename = strcat(SBJ_vars.dirs.proc,SBJ,'_actv_ROI_',stat_vars.an_id_r,'_mn',num2str(stat_vars.actv_win),'.mat');
    tmp = load(actv_filename,'actv_ch','actv_ch_epochs');
    actv_ch{2} = tmp.actv_ch; actv_ch_epochs{2} = tmp.actv_ch_epochs;
    good_e = union(actv_ch{1},actv_ch{2});
else
    good_e = hfa.label;
end
if stat_vars.actv_epochs_overlap
    error('good idea, but not yet implemented...');
end
if plt_vars.exclude_FWM
    good_e = intersect(good_e,hfa.label(~strcmp(einfo(:,3),'FWM')));
end
if plt_vars.exclude_OUT
    good_e = intersect(good_e,hfa.label(~strcmp(einfo(:,3),'OUT')));
end

% Select Channels, Trials, and Epoch of Interest
cfgs = [];
cfgs.channel = good_e;
cfgs.latency = stat_vars.stat_lim;
cfgs.trials = [];
for cond_ix = 1:numel(cond_lab)
    cfgs.trials = union(cfgs.trials, find(cond_idx(cond_ix,:)==1));
end
hfa = ft_selectdata(cfgs,hfa);

% Clean up einfo after removing channels
good_ix = NaN(size(good_e));
for e_ix = 1:numel(good_e)
    good_ix(e_ix) = find(strcmp(einfo(:,1),good_e{e_ix}));
end
einfo = einfo(good_ix,:);

% Clean up condition index after removing unused trials
cond_idx = false([length(cond_lab) length(cfgs.trials)]);
for cond_ix = 1:length(cond_lab)
    % Get binary condition index
    cond_idx(cond_ix,:) = logical(fn_condition_index(cond_lab{cond_ix},...
        trial_info.condition_n(cfgs.trials)));
end

% Sort by ROI
einfo = sortrows(einfo,plt_vars.sort_vec);
ch_roi_idx = zeros([numel(hfa.label) 1]);
for ch_ix = 1:numel(hfa.label)
    ch_roi_idx(ch_ix) = find(strcmp(hfa.label,einfo{ch_ix,1}));
end
hfa.label     = hfa.label(ch_roi_idx);
hfa.powspctrm = hfa.powspctrm(:,ch_roi_idx,:,:);

%% Smooth HFA time series
for ch_ix = 1:numel(hfa.label)
    if stat_vars.lp_flag && stat_vars.hp_flag
        hfa.powspctrm(:,ch_ix,1,:) = fn_EEGlab_bandpass(...
            hfa.powspctrm(:,ch_ix,1,:), proc_vars.resample_freq, stat_vars.hp_freq, stat_vars.lp_freq);
    elseif stat_vars.lp_flag
        hfa.powspctrm(:,ch_ix,1,:) = fn_EEGlab_lowpass(...
            squeeze(hfa.powspctrm(:,ch_ix,1,:)), proc_vars.resample_freq, stat_vars.lp_freq);
    elseif stat_vars.hp_flag
        error('Why are you only high passing?');
    else
        error('Why did you say yes smooth but no to both low and high pass?');
    end
end

%% Correlation
% Compute correlation
corr_mat = {};
hfa_cond = {};
for cond_ix = 1:numel(cond_lab)
    % Divide by condition
    cfgc = [];
    cfgc.trials = find(cond_idx(cond_ix,:)==1);
    hfa_cond{cond_ix} = ft_selectdata(cfgc,hfa);
    % NOTE: reshape goes down columns, so need to move time to first dimension
    corr_mat{cond_ix} = corrcoef(reshape(permute(hfa_cond{cond_ix}.powspctrm,[4 1 2 3]),...    % make dim_ord: time, trials, ch, empty
        size(hfa_cond{cond_ix}.powspctrm,1)*size(hfa_cond{cond_ix}.powspctrm,4),...            % n_trials*n_samples
        size(hfa_cond{cond_ix}.powspctrm,2)));                                                 % channels
end
% ok next part is only for CI...
diff_mat = corr_mat{2}-corr_mat{1}; % inc-con, so greater inc is positive

% Null Distribution
surr_diff = NaN(stat_vars.n_boots, numel(hfa.label), numel(hfa.label));
fprintf('============================= Bootstrap Correlations: %s =============================\n',stat_vars.event_lab);
for boot_ix = 1:stat_vars.n_boots
    fprintf('%i..',boot_ix);
    if mod(boot_ix,25)==0
        fprintf('\n');
    end
    % Randomly shuffle trials for each channel
    shuffled_data = hfa.powspctrm;
    for ch_ix = 1:size(shuffled_data,2)
        shuffled_data(:,ch_ix,1,:) = shuffled_data(randi(size(shuffled_data,1),size(shuffled_data,1),1),ch_ix,1,:);
    end
    
    surr_mat = {};
    for cond_ix = 1:numel(cond_lab)
        % Divide by condition
        shuffled_cond = shuffled_data(cond_idx(cond_ix,:),:,1,:);
        surr_mat{cond_ix} = corrcoef(reshape(permute(shuffled_cond,[4 1 2 3]),...   % dim_ord: time, trials, ch, empty
            size(shuffled_cond,1)*size(shuffled_cond,4),...                             % n_trials*n_samples
            size(shuffled_cond,2)));                                                    % channels
    end
    surr_diff(boot_ix,:,:) = surr_mat{2}-surr_mat{1};
end
fprintf('\n');  % for a nice reset of the command line

% Compute Statistics
pairs = nchoosek(1:numel(hfa.label),2);
n_pairs = size(pairs,1);
pvals = NaN([1 n_pairs]);
diff_vals = NaN([1 n_pairs]);
for pair_ix = 1:n_pairs
    null_corr_dist = sort(abs(surr_diff(:,pairs(pair_ix,1),pairs(pair_ix,2))));                   % Sort absolute values
    corr_logical = logical(abs(diff_mat(pairs(pair_ix,1),pairs(pair_ix,2))) > null_corr_dist);   % Find n_boots less than real value
    pvals(pair_ix) = (numel(null_corr_dist)-numel(find(corr_logical)))/stat_vars.n_boots;                         % convert to proportion of n_boots
    
    diff_vals(pair_ix) = diff_mat(pairs(pair_ix,1),pairs(pair_ix,2));    % Grab to set plotting limits
end

% Find epochs with significant task activations
%     [~, qvals] = mafdr(pvals); % Errors on some random ch during SBJ08a_HFA_actv (e.g., ROF8-9 in IR32), so try this:
%     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
[~, ~, ~, qvals] = fdr_bh(pvals);

%% Save Data
data_out_filename = [SBJ_vars.dirs.proc SBJ '_' stat_id '_' an_id_main '.mat'];
save(data_out_filename, '-v7.3', 'hfa_cond','pairs','diff_mat','diff_vals','qvals');

end