function SBJ08ab_HFA_stats(SBJ,an_id,stat_id)
% Computes cluster-based statistics for given conditions based on saved HFA
% clear all; %close all;
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

% Select Conditions of Interest
[cond_lab, ~, ~] = fn_condition_label_styles(st.conditions);

%% Separate Trials by Condition
hfa_all = hfa;
hfa     = {};
n_trials = zeros([1 numel(cond_lab)]);
for cond_ix = 1:numel(cond_lab)
    cfgs = [];
    cfgs.trials = find(fn_condition_index(cond_lab{cond_ix}, trial_info.condition_n));
    hfa{cond_ix} = ft_selectdata(cfgs,hfa_all);
    
    % Grab n_trials for design matrix
    n_trials(cond_ix) = size(hfa{cond_ix}.trialinfo,1);
end

%% Run Statistics
fprintf('===================================================\n');
fprintf('--------------------- Statistics ------------------\n');
fprintf('===================================================\n');
% Create design matrix
% design = zeros(2,size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1));
% % Conditions (Independent Variable)
% design(1,1:size(roi_erp_con.trial,1)) = 1;
% design(1,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
% % Trial Numbers
% design(2,1:size(roi_erp_con.trial,1)) = 1;
% design(2,(size(roi_erp_con.trial,1)+1):(size(roi_erp_con.trial,1) + size(roi_erp_inc.trial,1)))= 2;
design = zeros(2,sum(n_trials));
for cond_ix = 1:numel(cond_lab)
    if cond_ix==1
        design(1,1:n_trials(cond_ix)) = cond_ix;                                % Conditions (Independent Variable)
        design(2,1:n_trials(cond_ix)) = 1:n_trials(cond_ix);                    % Trial Numbers
    else
        design(1,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= cond_ix; % Conditions (Independent Variable)
        design(2,sum(n_trials(1:cond_ix-1))+1:sum(n_trials(1:cond_ix)))= 1:n_trials(cond_ix);
    end
end

% Prepare neighbors layout
% cfgn = [];
% cfgn.method  = 'distance';
% cfgn.layout  = 'ordered';
% cfgn.channel = elecs;
% neighbors    = ft_prepare_neighbours(cfgn,roi_erp_allch{1});
% for ch_ix = 1:numel(roi_erp{1}.label)
%     neighbors(ch_ix).label = roi_erp{1}.label{ch_ix};
%     neighbors(ch_ix).neighblabel = {};
% end

% Calculate statistics
cfg_stat.design           = design;
[stat] = ft_freqstatistics(cfg_stat, hfa{:});

%% Save Results
out_fname = strcat(SBJ_vars.dirs.proc,SBJ,'_ROI_',an_id,'_',stat_id,'.mat');
fprintf('===================================================\n');
fprintf('--- Saving %s ------------------\n',out_fname);
fprintf('===================================================\n');
save(out_fname,'-v7.3','stat','st');

end
