function SBJ11a_corr_HFA_acE_surr(SBJ,pipeline_id,stat_id,an_id_main,atlas_id,roi_id,plt_id)
% Build connectivity matrix based on HFA correlations
%   non-parametric stats via circular shift of trial time series

%% Data Preparation
% Set up paths
if exist('/home/knight/hoycw/','dir');root_dir='/home/knight/hoycw/';ft_dir=[root_dir 'Apps/fieldtrip/'];
else root_dir='/Volumes/hoycw_clust/';ft_dir='/Users/colinhoy/Code/Apps/fieldtrip/';end
addpath([root_dir 'PRJ_Stroop/scripts/']);
addpath([root_dir 'PRJ_Stroop/scripts/utils/']);
addpath(ft_dir);
ft_defaults

%% Load Results
eval(['run ' root_dir 'PRJ_Stroop/scripts/SBJ_vars/' SBJ '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/stat_vars/' stat_id '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/proc_vars/' pipeline_id '_proc_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/an_vars/' an_id_main '_vars.m']);
eval(['run ' root_dir 'PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

% Load data
hfa_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id_main,'.mat');
load(hfa_fname,'hfa');

%% Load ROI and GM/WM info
elec_fname = [SBJ_vars.dirs.recon SBJ '_elec_' pipeline_id '_pat_' atlas_id '.mat'];
load(elec_fname);

% Sort elecs by stat labels
cfgs = []; cfgs.channel = hfa.label;
elec = fn_select_elec(cfgs,elec);
elec.roi = fn_atlas2roi_labels(elec.atlas_label,atlas_id,roi_id);
elec.roi_id = roi_id;

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
    good_e = intersect(good_e,hfa.label(~strcmp(roi_lab,'FWM')));
end
if plt_vars.exclude_OUT
    good_e = intersect(good_e,hfa.label(~strcmp(roi_lab,'OUT')));
end

% Select Channels and Epoch of Interest
cfgs = [];
cfgs.channel = good_e;
elec = fn_select_elec(cfgs,elec);
cfgs.latency = stat_vars.stat_lim;
hfa = ft_selectdata(cfgs,hfa);

% Sort by ROI
[~,roi_sort_idx] = sort(elec.roi);
hfa.label     = hfa.label(roi_sort_idx);
hfa.powspctrm = hfa.powspctrm(:,roi_sort_idx,:,:);
elec_fields = fieldnames(elec);
for field_ix = 1:numel(elec_fields)
    if any(size(eval(['elec.' elec_fields{field_ix}]))==numel(elec.label))
        eval(['elec.' elec_fields{field_ix} ' = elec.' elec_fields{field_ix} '(roi_sort_idx,:);']);
    end
end

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
% NOTE: reshape goes down columns, so need to move time to first dimension
corr_mat = corrcoef(reshape(permute(hfa.powspctrm,[4 1 2 3]),...    % make dim_ord: time, trials, ch, empty
    size(hfa.powspctrm,1)*size(hfa.powspctrm,4),...                 % n_trials*n_samples
    size(hfa.powspctrm,2)));                                        % channels

% Null Distribution
surr_mat = NaN(stat_vars.n_boots, numel(hfa.label), numel(hfa.label));
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
    
    surr_mat(boot_ix,:,:) = corrcoef(reshape(permute(shuffled_data,[4 1 2 3]),...   % dim_ord: time, trials, ch, empty
        size(shuffled_data,1)*size(shuffled_data,4),...                             % n_trials*n_samples
        size(shuffled_data,2)));                                                    % channels
end
fprintf('\n');  % for a nice reset of the command line

% Compute Statistics
pairs = nchoosek(1:numel(hfa.label),2);
n_pairs = size(pairs,1);
pvals = NaN([1 n_pairs]);
corr_vals = NaN([1 n_pairs]);
for pair_ix = 1:n_pairs
    null_corr_dist = sort(abs(surr_mat(:,pairs(pair_ix,1),pairs(pair_ix,2))));                   % Sort absolute values
    corr_logical = logical(abs(corr_mat(pairs(pair_ix,1),pairs(pair_ix,2))) > null_corr_dist);   % Find n_boots less than real value
    pvals(pair_ix) = (numel(null_corr_dist)-numel(find(corr_logical)))/stat_vars.n_boots;                         % convert to proportion of n_boots
    
    corr_vals(pair_ix) = corr_mat(pairs(pair_ix,1),pairs(pair_ix,2));    % Grab to set plotting limits
end
% Find epochs with significant task activations
%     [~, qvals] = mafdr(pvals); % Errors on some random ch during SBJ08a_HFA_actv (e.g., ROF8-9 in IR32), so try this:
%     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
[~, ~, ~, qvals] = fdr_bh(pvals);

%% Save Data
data_out_filename = [SBJ_vars.dirs.proc SBJ '_' stat_id '_' an_id_main '_' roi_id '_' atlas_id '.mat'];
save(data_out_filename, '-v7.3', 'hfa','pairs','corr_mat','corr_vals','qvals','elec');

end
