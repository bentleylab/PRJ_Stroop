function SBJ11a_corr_HFA_acE_surr(SBJ,pipeline_id,stat_id,an_id,plt_id)
% Build connectivity matrix based on HFA correlations
%   non-parametric stats via circular shift of trial time series
error('stop running response-locked, it doesnt make sense');

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
eval(['run /home/knight/hoycw/PRJ_Stroop/scripts/plt_vars/' plt_id '_vars.m']);

event_lab = {'stim', 'resp'};

% Load data
hfa_filename1 = strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id_s,'.mat');
hfa_filename2 = strcat(SBJ_vars.dirs.proc,SBJ,'_HFA_ROI_',an_id_r,'.mat');
tmp = load(hfa_filename1,'hfa'); hfa{1} = tmp.hfa;
tmp = load(hfa_filename2,'hfa'); hfa{2} = tmp.hfa;
clear tmp

% Load ROI and GM/WM info
einfo_filename = [SBJ_vars.dirs.preproc SBJ '_einfo_' pipeline_id '.mat'];
load(einfo_filename);
einfo = sortrows(einfo,plt_vars.sort_vec);
% Electrode Info Table:
%   label- name of electrode
%   ROI- specific region
%   gROI- general region (LPFC, MPFC, OFC, FWM=frontal white matter)
%   ROI2- specific region of second electrode
%   tissue- primary tissue type
%   GM weight- percentage of electrode pair in GM
%   Out- 0/1 flag for whether this is partially out of the brain

%% Prep Data
% Trim data to stats epoch
cfg_trim = [];
cfg_trim.latency = stat_vars.stat_lim_S;
hfa{1} = ft_selectdata(cfg_trim,hfa{1});
cfg_trim.latency = stat_vars.stat_lim_R;
hfa{2} = ft_selectdata(cfg_trim,hfa{2});

for sr_ix = 1:2
    % Sort by ROI
    ch_roi_idx = zeros([numel(hfa{sr_ix}.label) 1]);
    for ch_ix = 1:numel(hfa{sr_ix}.label)
        ch_roi_idx(ch_ix) = find(strcmp(hfa{sr_ix}.label,einfo{ch_ix,1}));
    end
    hfa{sr_ix}.label     = hfa{sr_ix}.label(ch_roi_idx);
    hfa{sr_ix}.powspctrm = hfa{sr_ix}.powspctrm(:,ch_roi_idx,:,:);
    
    % Smooth HFA time series
    for ch_ix = 1:numel(hfa{sr_ix}.label)
        if stat_vars.lp_flag && stat_vars.hp_flag
            hfa{sr_ix}.powspctrm(:,ch_ix,1,:) = fn_EEGlab_bandpass(...
                hfa{sr_ix}.powspctrm(:,ch_ix,1,:), proc_vars.resample_freq, stat_vars.hp_freq, stat_vars.lp_freq);
        elseif stat_vars.lp_flag
            hfa{sr_ix}.powspctrm(:,ch_ix,1,:) = fn_EEGlab_lowpass(...
                squeeze(hfa{sr_ix}.powspctrm(:,ch_ix,1,:)), proc_vars.resample_freq, stat_vars.lp_freq);
        elseif stat_vars.hp_flag
            error('Why are you only high passing?');
        else
            error('Why did you say yes smooth but no to both low and high pass?');
        end
    end
end

%% Correlation
for sr_ix = 1:2
    % Compute correlation
    corr_mat{sr_ix} = corrcoef(reshape(hfa{sr_ix}.powspctrm,...             % full data
        size(hfa{sr_ix}.powspctrm,1)*size(hfa{sr_ix}.powspctrm,4),...       % n_trials*n_samples
        size(hfa{sr_ix}.powspctrm,2)));                                     % channels
    
    % Null Distribution
    surr_mat{sr_ix} = NaN(stat_vars.n_boots, numel(hfa{sr_ix}.label), numel(hfa{sr_ix}.label));
    fprintf('============================= Bootstrap Correlations: %s =============================\n',event_lab{sr_ix});
    for boot_ix = 1:stat_vars.n_boots
        fprintf('%i..',boot_ix);
        if mod(boot_ix,25)==0
            fprintf('\n');
        end
        % Randomly shuffle trials for each channel
        shuffled_data = hfa{sr_ix}.powspctrm;
        for ch_ix = 1:size(shuffled_data,2)
            shuffled_data(:,ch_ix,1,:) = shuffled_data(randi(size(shuffled_data,1),size(shuffled_data,1),1),ch_ix,1,:);
        end
        
        surr_mat{sr_ix}(boot_ix,:,:) = corrcoef(reshape(shuffled_data,...             % full data
            size(shuffled_data,1)*size(shuffled_data,4),...       % n_trials*n_samples
            size(shuffled_data,2)));                                     % channels
        fprintf('\n');
    end
    
    % Compute Statistics
    pairs = nchoosek(1:numel(hfa{sr_ix}.label),2);
    n_pairs = size(pairs,1);
    pvals = NaN([1 n_pairs]);
    corr_vals{sr_ix} = NaN([1 n_pairs]);
    for pair_ix = 1:n_pairs
        null_corr_dist = sort(abs(surr_mat{sr_ix}(:,pairs(pair_ix,1),pairs(pair_ix,2))));                   % Sort absolute values
        corr_logical = logical(abs(corr_mat{sr_ix}(pairs(pair_ix,1),pairs(pair_ix,2))) > null_corr_dist);   % Find n_boots less than real value
        pvals(pair_ix) = (numel(null_corr_dist)-numel(find(corr_logical)))/stat_vars.n_boots;                         % convert to proportion of n_boots
        
        corr_vals{sr_ix}(pair_ix) = corr_mat{sr_ix}(pairs(pair_ix,1),pairs(pair_ix,2));    % Grab to set plotting limits
    end
    % Find epochs with significant task activations
    %     [~, qvals] = mafdr(pvals); % Errors on some random ch during SBJ08a_HFA_actv (e.g., ROF8-9 in IR32), so try this:
    %     [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
    [~, ~, ~, qvals{sr_ix}] = fdr_bh(pvals);
end

%% Save Data
data_out_filename = [SBJ_vars.dirs.proc SBJ '_' stat_id '_' an_id_s '-' an_id_r '.mat'];
save(data_out_filename, '-v7.3', 'hfa','pairs','corr_mat','corr_vals','qvals');

end
